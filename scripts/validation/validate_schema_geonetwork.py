#!/usr/bin/env python3
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
import yaml


def utc_now_iso():
    return datetime.now(timezone.utc).isoformat()


def load_config(config_file):
    return yaml.safe_load(Path(config_file).read_text(encoding="utf-8")) or {}


def load_session(session_file):
    payload = json.loads(Path(session_file).read_text(encoding="utf-8"))
    xsrf_token = payload.get("XSRF-TOKEN")
    cookies = payload.get("cookies", {})
    api_user = payload.get("authenticated_user")
    if not xsrf_token:
        raise ValueError("No se encontró XSRF-TOKEN en la sesión de GeoNetwork.")
    return xsrf_token, cookies, api_user


def load_geonetwork_url(config):
    geonetwork_url = config.get("services", {}).get("geonetwork", {}).get("url", "")
    if not geonetwork_url:
        raise ValueError("No se encontró services.geonetwork.url en config.yaml.")
    return geonetwork_url.rstrip("/")


def load_schema_validation_settings(config):
    schema_cfg = (
        config.get("validation", {}).get("geonetwork_schema", {})
        if isinstance(config.get("validation"), dict)
        else {}
    )
    if not isinstance(schema_cfg, dict):
        schema_cfg = {}

    retry_cfg = schema_cfg.get("retry", {})
    if not isinstance(retry_cfg, dict):
        retry_cfg = {}

    enabled = schema_cfg.get("enabled", True) is True
    max_attempts = int(retry_cfg.get("max_attempts", 6))
    sleep_seconds = float(retry_cfg.get("sleep_seconds", 10))

    if max_attempts < 1:
        max_attempts = 1
    if sleep_seconds < 0:
        sleep_seconds = 0

    return {
        "enabled": enabled,
        "max_attempts": max_attempts,
        "sleep_seconds": sleep_seconds,
    }


def build_record_candidates(upload_response_file, record_id_file):
    candidates = []
    seen = set()

    def add_candidate(value):
        if value is None:
            return
        candidate = str(value).strip()
        if not candidate:
            return
        if candidate in seen:
            return
        seen.add(candidate)
        candidates.append(candidate)

    # Canonical workflow UUID.
    if Path(record_id_file).exists():
        add_candidate(Path(record_id_file).read_text(encoding="utf-8").strip())

    upload_response = {}
    if Path(upload_response_file).exists():
        upload_response = json.loads(Path(upload_response_file).read_text(encoding="utf-8"))

    png = upload_response.get("png_response", {}) if isinstance(upload_response, dict) else {}
    if isinstance(png, dict):
        add_candidate(png.get("metadataUuid"))
        add_candidate(png.get("metadataId"))

    if isinstance(upload_response, dict):
        add_candidate(upload_response.get("metadataUuid"))
        add_candidate(upload_response.get("metadataId"))

    if not candidates:
        raise ValueError("No se pudo resolver metadataUuid/metadataId para validar esquema.")
    return candidates


def parse_json_response(response):
    try:
        return response.json()
    except ValueError:
        return {"raw_response": response.text}


def extract_xsd_messages(payload, max_items=5):
    messages = []
    if not isinstance(payload, dict):
        return messages
    reports = payload.get("report")
    if not isinstance(reports, list):
        return messages
    for report in reports:
        if not isinstance(report, dict):
            continue
        if str(report.get("id", "")).lower() != "xsd":
            continue
        patterns = report.get("patterns", {}).get("pattern", [])
        if not isinstance(patterns, list):
            continue
        for pattern in patterns:
            if not isinstance(pattern, dict):
                continue
            title = str(pattern.get("title", "")).strip()
            if title:
                messages.append(title)
            rules = pattern.get("rules", {}).get("rule", [])
            if isinstance(rules, list):
                for rule in rules:
                    if not isinstance(rule, dict):
                        continue
                    details = str(rule.get("details", "")).strip()
                    if details:
                        messages.append(details)
            if len(messages) >= max_items:
                return messages[:max_items]
    return messages[:max_items]


def evaluate_schema_payload(payload):
    """
    Expected successful GeoNetwork payload (JSON):
    {
      "report": [{"id":"xsd","error":0,...}, ...]
    }
    """
    failures = []

    if not isinstance(payload, dict):
        failures.append({"path": "payload", "reason": "Formato de respuesta inesperado."})
        return False, failures

    reports = payload.get("report")
    if not isinstance(reports, list):
        failures.append({"path": "report", "reason": "No se encontró lista 'report'."})
        return False, failures

    xsd_report = None
    for report in reports:
        if isinstance(report, dict) and str(report.get("id", "")).lower() == "xsd":
            xsd_report = report
            break

    if xsd_report is None:
        failures.append({"path": "report", "reason": "No se encontró entrada de validación 'xsd'."})
        return False, failures

    try:
        xsd_errors = int(xsd_report.get("error", 0))
    except (TypeError, ValueError):
        failures.append({"path": "report.xsd.error", "reason": "Valor 'error' no numérico."})
        return False, failures

    if xsd_errors > 0:
        failures.append(
            {
                "path": "report.xsd.error",
                "reason": f"Se detectaron {xsd_errors} errores de esquema XSD.",
            }
        )
        return False, failures

    return True, failures


def validate_schema(config_file, session_file, record_id_file, upload_response_file, output_file):
    config = load_config(config_file)
    geonetwork_url = load_geonetwork_url(config)
    settings = load_schema_validation_settings(config)
    xsrf_token, cookies, api_user = load_session(session_file)
    record_candidates = build_record_candidates(upload_response_file, record_id_file)

    if not settings["enabled"]:
        report = {
            "timestamp": utc_now_iso(),
            "status": "skipped",
            "message": "validation.geonetwork_schema.enabled = false",
            "record_candidates": record_candidates,
            "schema_valid": True,
            "settings": settings,
            "api_user": api_user,
        }
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        Path(output_file).write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        return report

    session = requests.Session()
    session.cookies.update(cookies)
    api_base = f"{geonetwork_url}/srv/api"

    attempts_log = []
    last_report = None
    companion_path = Path(output_file).with_name(Path(output_file).stem + ".last.json")

    for attempt in range(1, settings["max_attempts"] + 1):
        attempt_entry = {"attempt": attempt, "candidates": []}
        selected_identifier = None

        for candidate in record_candidates:
            record_url = f"{api_base}/records/{candidate}"
            record_response = session.get(
                record_url,
                headers={"X-XSRF-TOKEN": xsrf_token, "Accept": "application/json"},
                timeout=30,
            )

            candidate_entry = {
                "record_identifier": candidate,
                "record_check_status": record_response.status_code,
            }

            if record_response.status_code == 200 and selected_identifier is None:
                selected_identifier = candidate

            attempt_entry["candidates"].append(candidate_entry)

        if selected_identifier is None:
            attempts_log.append(attempt_entry)
            if attempt < settings["max_attempts"]:
                time.sleep(settings["sleep_seconds"])
                continue

            last_report = {
                "timestamp": utc_now_iso(),
                "status": "failed",
                "error": "Ningún identificador de registro estuvo disponible en la API de GeoNetwork.",
                "record_candidates": record_candidates,
                "attempts": attempts_log,
                "settings": settings,
                "schema_valid": False,
                "api_user": api_user,
            }
            break

        editor_url = f"{api_base}/records/{selected_identifier}/editor"
        validate_url = f"{api_base}/records/{selected_identifier}/validate/internal"

        editor_response = session.get(
            editor_url,
            headers={"X-XSRF-TOKEN": xsrf_token},
            timeout=30,
        )
        attempt_entry["editor_lock_status"] = editor_response.status_code

        validate_response = session.put(
            validate_url,
            headers={"X-XSRF-TOKEN": xsrf_token, "Accept": "application/json"},
            timeout=90,
        )
        payload = parse_json_response(validate_response)

        attempt_entry["selected_identifier"] = selected_identifier
        attempt_entry["validate_status"] = validate_response.status_code

        unlock_status = None
        if editor_response.status_code in (200, 201, 204):
            unlock_response = session.delete(
                editor_url,
                headers={"X-XSRF-TOKEN": xsrf_token},
                timeout=30,
            )
            unlock_status = unlock_response.status_code
        attempt_entry["editor_unlock_status"] = unlock_status

        if validate_response.status_code in (200, 201, 204):
            schema_valid, schema_failures = evaluate_schema_payload(payload)
            xsd_messages = extract_xsd_messages(payload)
            last_report = {
                "timestamp": utc_now_iso(),
                "status": "ok" if schema_valid else "failed",
                "record_identifier": selected_identifier,
                "validate_url": validate_url,
                "http_status": validate_response.status_code,
                "schema_valid": schema_valid,
                "schema_failures": schema_failures,
                "xsd_messages": xsd_messages,
                "payload": payload,
                "attempts": attempts_log + [attempt_entry],
                "settings": settings,
                "api_user": api_user,
                "record_candidates": record_candidates,
            }
            attempts_log.append(attempt_entry)
            break

        # Retry only for transient visibility issues.
        if validate_response.status_code == 404 and attempt < settings["max_attempts"]:
            attempts_log.append(attempt_entry)
            time.sleep(settings["sleep_seconds"])
            continue

        last_report = {
            "timestamp": utc_now_iso(),
            "status": "failed",
            "record_identifier": selected_identifier,
            "validate_url": validate_url,
            "http_status": validate_response.status_code,
            "schema_valid": False,
            "schema_failures": [],
            "payload": payload,
            "attempts": attempts_log + [attempt_entry],
            "settings": settings,
            "api_user": api_user,
            "record_candidates": record_candidates,
            "error": f"GeoNetwork devolvió HTTP {validate_response.status_code} en validate/internal.",
        }
        attempts_log.append(attempt_entry)
        break

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(output_file).write_text(json.dumps(last_report, indent=2, ensure_ascii=False), encoding="utf-8")
    companion_path.parent.mkdir(parents=True, exist_ok=True)
    companion_path.write_text(json.dumps(last_report, indent=2, ensure_ascii=False), encoding="utf-8")

    if last_report.get("status") != "ok":
        xsd_messages = last_report.get("xsd_messages", [])
        detail = f" Primer error: {xsd_messages[0]}" if xsd_messages else ""
        raise RuntimeError(
            "La validación de esquema en GeoNetwork falló. "
            "Consulta logs/metadata_schema_validation.last.json."
            + detail
        )

    return last_report


if __name__ == "__main__":
    try:
        config_file = snakemake.input.config
        session_file = snakemake.input.session
        record_id_file = snakemake.input.record_id
        upload_response_file = snakemake.input.upload_response
        output_file = snakemake.output.report
        gate_file = snakemake.output.ok
    except NameError:
        if len(sys.argv) < 6:
            print(
                "Uso: validate_schema_geonetwork.py "
                "<config.yaml> <session.json> <record_id.txt> <upload_response.json> <output_report.json> [gate.txt]"
            )
            sys.exit(1)
        config_file = sys.argv[1]
        session_file = sys.argv[2]
        record_id_file = sys.argv[3]
        upload_response_file = sys.argv[4]
        output_file = sys.argv[5]
        gate_file = sys.argv[6] if len(sys.argv) > 6 else "temp_files/metadata_schema_ok.txt"

    try:
        report = validate_schema(
            config_file=config_file,
            session_file=session_file,
            record_id_file=record_id_file,
            upload_response_file=upload_response_file,
            output_file=output_file,
        )
        marker = "SKIPPED\n" if report.get("status") == "skipped" else "OK\n"
        Path(gate_file).parent.mkdir(parents=True, exist_ok=True)
        Path(gate_file).write_text(marker, encoding="utf-8")
        print("Validación de esquema GeoNetwork completada.")
    except Exception as exc:
        report_path = Path(output_file)
        if not report_path.exists():
            report_path.parent.mkdir(parents=True, exist_ok=True)
            fallback_report = {
                "timestamp": utc_now_iso(),
                "status": "error",
                "schema_valid": False,
                "error": str(exc),
            }
            report_path.write_text(
                json.dumps(fallback_report, indent=2, ensure_ascii=False), encoding="utf-8"
            )
            report_path.with_name(report_path.stem + ".last.json").write_text(
                json.dumps(fallback_report, indent=2, ensure_ascii=False), encoding="utf-8"
            )
        Path(gate_file).parent.mkdir(parents=True, exist_ok=True)
        Path(gate_file).write_text("FAIL\n", encoding="utf-8")
        print(f"Error en validación de esquema GeoNetwork: {exc}", file=sys.stderr)
        sys.exit(1)
