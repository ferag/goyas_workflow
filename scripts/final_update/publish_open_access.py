#!/usr/bin/env python3
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests
import yaml


def utc_now_iso():
    return datetime.now(timezone.utc).isoformat()


def write_text_atomic(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(content, encoding="utf-8")
    tmp.replace(path)


def write_json_atomic(path, payload):
    write_text_atomic(path, json.dumps(payload, indent=2, ensure_ascii=False))


def _first_non_empty(*values):
    for value in values:
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def read_gate_value(path):
    if not Path(path).exists():
        return ""
    return Path(path).read_text(encoding="utf-8").strip().upper()


def load_config(path):
    return yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}


def load_session(path):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    xsrf_token = payload.get("XSRF-TOKEN")
    cookies = payload.get("cookies", {})
    if not xsrf_token:
        raise ValueError("No se encontró XSRF-TOKEN en la sesión de GeoNetwork.")
    return xsrf_token, cookies


def ensure_geonetwork_url(config):
    geonetwork_url = config.get("services", {}).get("geonetwork", {}).get("url", "")
    geonetwork_url = geonetwork_url.strip().rstrip("/")
    if not geonetwork_url:
        raise ValueError("No se encontró services.geonetwork.url en config.yaml.")
    return geonetwork_url


def resolve_metadata_uuid(upload_response_file, record_id_file):
    record_id = ""
    if Path(record_id_file).exists():
        record_id = Path(record_id_file).read_text(encoding="utf-8").strip()

    upload_data = {}
    if Path(upload_response_file).exists():
        upload_data = json.loads(Path(upload_response_file).read_text(encoding="utf-8"))

    png = upload_data.get("png_response", {}) if isinstance(upload_data, dict) else {}
    tif = upload_data.get("tif_response", {}) if isinstance(upload_data, dict) else {}
    return _first_non_empty(
        png.get("metadataUuid"),
        png.get("metadataId"),
        tif.get("metadataUuid"),
        tif.get("metadataId"),
        upload_data.get("metadataUuid") if isinstance(upload_data, dict) else "",
        upload_data.get("metadataId") if isinstance(upload_data, dict) else "",
        record_id,
    )


def publish_record(geonetwork_url, session_file, metadata_uuid):
    xsrf_token, cookies = load_session(session_file)
    session = requests.Session()
    session.cookies.update(cookies)

    publish_url = f"{geonetwork_url}/srv/api/records/{metadata_uuid}/publish"
    headers = {"X-XSRF-TOKEN": xsrf_token, "Accept": "application/json"}
    response = session.put(publish_url, headers=headers, timeout=60)

    payload = None
    try:
        payload = response.json()
    except ValueError:
        payload = {"raw_response": response.text}

    if response.status_code not in (200, 201, 204):
        raise RuntimeError(
            f"Error publicando registro ({response.status_code}) en {publish_url}: {payload}"
        )

    return {
        "publish_url": publish_url,
        "http_status": response.status_code,
        "response": payload,
    }


def run(config_file, session_file, record_id_file, upload_response_file, schema_gate_file):
    config = load_config(config_file)
    enabled = config.get("open_access", False) is True
    schema_gate = read_gate_value(schema_gate_file)

    report = {
        "timestamp": utc_now_iso(),
        "enabled": enabled,
        "schema_gate": schema_gate,
        "status": "skipped",
        "metadata_uuid": None,
        "errors": [],
    }

    if not enabled:
        report["reason"] = "open_access = false"
        return report

    if schema_gate != "OK":
        report["reason"] = "Requiere metadata_schema_ok=OK"
        return report

    metadata_uuid = resolve_metadata_uuid(upload_response_file, record_id_file)
    if not metadata_uuid:
        raise ValueError("No se pudo resolver metadata UUID para publicar el registro.")
    report["metadata_uuid"] = metadata_uuid

    geonetwork_url = ensure_geonetwork_url(config)
    publish_result = publish_record(geonetwork_url, session_file, metadata_uuid)
    report["publish"] = publish_result
    report["status"] = "ok"
    return report


if __name__ == "__main__":
    try:
        if "snakemake" in globals():
            config_file = str(snakemake.input.config)
            session_file = str(snakemake.input.session)
            record_id_file = str(snakemake.input.record_id)
            upload_response_file = str(snakemake.input.upload_response)
            schema_gate_file = str(snakemake.input.schema_gate)
            report_file = Path(str(snakemake.output.report))
            gate_file = Path(str(snakemake.output.ok))
        else:
            if len(sys.argv) < 8:
                print(
                    "Uso: publish_open_access.py "
                    "<config.yaml> <session.json> <record_id.txt> <upload_response.json> "
                    "<schema_gate.txt> <report.json> <gate.txt>"
                )
                sys.exit(1)
            config_file = sys.argv[1]
            session_file = sys.argv[2]
            record_id_file = sys.argv[3]
            upload_response_file = sys.argv[4]
            schema_gate_file = sys.argv[5]
            report_file = Path(sys.argv[6])
            gate_file = Path(sys.argv[7])

        result = run(
            config_file=config_file,
            session_file=session_file,
            record_id_file=record_id_file,
            upload_response_file=upload_response_file,
            schema_gate_file=schema_gate_file,
        )
        write_json_atomic(report_file, result)
        marker = "OK\n" if result.get("status") == "ok" else "SKIPPED\n"
        write_text_atomic(gate_file, marker)
        print(f"Publicación open_access completada con estado: {result.get('status')}")
        if result.get("status") == "failed":
            sys.exit(1)
    except Exception as exc:
        error_report = {
            "timestamp": utc_now_iso(),
            "status": "failed",
            "errors": [str(exc)],
        }
        if "report_file" in locals():
            write_json_atomic(report_file, error_report)
        if "gate_file" in locals():
            write_text_atomic(gate_file, "FAIL\n")
        print(f"Error en publish_open_access.py: {exc}", file=sys.stderr)
        sys.exit(1)
