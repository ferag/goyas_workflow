#!/usr/bin/env python3
import importlib
import json
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


DEFAULT_TESTS = [
    "rda_a1_01m",
    "rda_a1_02d",
    "rda_a1_04m",
    "rda_a1_1_01m",
    "rda_a1_2_01d",
    "rda_f1_01d",
    "rda_f1_01m",
    "rda_f1_02d",
    "rda_f1_02m",
    "rda_f2_01m",
    "rda_f2_01m_disciplinar",
    "rda_f2_01m_generic",
    "rda_f3_01m",
    "rda_f4_01m",
    "rda_i1_01m",
    "rda_i3_01d",
    "rda_i3_02d",
    "rda_i3_02m",
    "rda_i3_03m",
    "rda_i3_04m",
    "rda_r1_01m",
    "rda_r1_1_01m",
    "rda_r1_1_02m",
]

FAIR_EVA_PLUGIN = "goyas"
FAIR_EVA_LANG = "en"
FAIR_EVA_MIN_SCORE = 75.0


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


def extract_item_id(upload_response_file, record_id_file):
    record_id = ""
    if record_id_file.exists():
        record_id = record_id_file.read_text(encoding="utf-8").strip()

    upload_data = {}
    if upload_response_file.exists():
        upload_data = json.loads(upload_response_file.read_text(encoding="utf-8"))

    upload_id = _first_non_empty(
        upload_data.get("png_response", {}).get("metadataId"),
        upload_data.get("tif_response", {}).get("metadataId"),
        upload_data.get("metadataId"),
    )

    item_id = _first_non_empty(upload_id, record_id)
    if not item_id:
        raise ValueError(
            "No se pudo obtener el identificador del registro (metadataId/record_id)."
        )

    return item_id


def resolve_api_endpoint(config, fair_cfg):
    configured = (fair_cfg.get("api_endpoint") or "").strip()
    if configured:
        return configured.rstrip("/")

    geonetwork_url = (
        config.get("services", {}).get("geonetwork", {}).get("url", "").strip()
    )
    if not geonetwork_url:
        raise ValueError(
            "No se encontró services.geonetwork.url y validation.fair_eva.api_endpoint no está definido."
        )
    return f"{geonetwork_url.rstrip('/')}/srv/api"


def _install_local_metadata_loader(plugin_module, metadata_xml_path):
    """
    Sustituye temporalmente get_metadata para que el plugin lea el XML local
    generado por el workflow, evitando dependencias de visibilidad pública del registro.
    """
    xml_path = Path(metadata_xml_path)
    if not xml_path.exists():
        return None

    original = plugin_module.Plugin.get_metadata

    def _get_metadata_local(self):
        xml_payload = xml_path.read_bytes()
        namespaces = self._collect_namespaces(xml_payload)
        root = ET.fromstring(xml_payload)
        rows = self._parse_xml(root, namespaces)
        if not rows:
            return self._empty_metadata_df()
        return pd.DataFrame(
            rows,
            columns=["metadata_schema", "element", "text_value", "qualifier"],
        ).drop_duplicates(ignore_index=True)

    plugin_module.Plugin.get_metadata = _get_metadata_local
    return original


def run_fair_eva_validation(config_file, upload_response_file, record_id_file, metadata_xml_file):
    cfg = yaml.safe_load(Path(config_file).read_text(encoding="utf-8")) or {}
    fair_cfg = cfg.get("validation", {}).get("fair_eva", {}) or {}

    enabled = fair_cfg.get("enabled", False) is True
    plugin_name = FAIR_EVA_PLUGIN
    lang = FAIR_EVA_LANG
    min_score = FAIR_EVA_MIN_SCORE
    tests = list(DEFAULT_TESTS)

    base_report = {
        "timestamp": utc_now_iso(),
        "enabled": enabled,
        "plugin": plugin_name,
        "lang": lang,
        "min_score": min_score,
        "tests_requested": tests,
        "api_endpoint": None,
        "status": "skipped",
        "item_id": None,
        "all_passed": True,
        "summary": {"total": 0, "passed": 0, "failed": 0, "errors": 0},
        "results": {},
        "errors": [],
    }

    if not enabled:
        base_report["status"] = "skipped"
        base_report["reason"] = "validation.fair_eva.enabled is false"
        return base_report

    api_endpoint = resolve_api_endpoint(cfg, fair_cfg)
    base_report["api_endpoint"] = api_endpoint

    item_id = extract_item_id(Path(upload_response_file), Path(record_id_file))
    base_report["item_id"] = item_id

    plugin_path = f"fair_eva.plugin.{plugin_name}"
    try:
        plugin_module = importlib.import_module(f"{plugin_path}.plugin")
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "No se pudo importar el plugin FAIR-EVA de GOYAS. "
            "Instala la dependencia en el entorno del workflow, por ejemplo: "
            "'pip install git+https://github.com/IFCA-Advanced-Computing/fair_eva_plugin_goyas'"
        ) from exc

    config_data = plugin_module.Plugin.load_config(plugin_path)
    original_get_metadata = _install_local_metadata_loader(plugin_module, metadata_xml_file)
    eva = plugin_module.Plugin(
        item_id=item_id,
        api_endpoint=api_endpoint,
        lang=lang,
        config=config_data,
        name=plugin_name,
    )
    if original_get_metadata is not None:
        plugin_module.Plugin.get_metadata = original_get_metadata

    passed = 0
    failed = 0
    errors = 0
    all_passed = True
    results = {}

    for test_id in tests:
        test_payload = {
            "test_id": test_id,
            "status": "error",
            "points": 0,
            "min_score": min_score,
            "passed": False,
            "messages": [],
        }
        method = getattr(eva, test_id, None)
        if not callable(method):
            test_payload["messages"] = [f"El test '{test_id}' no existe en el plugin '{plugin_name}'."]
            all_passed = False
            errors += 1
            results[test_id] = test_payload
            continue

        try:
            points, messages = method()
            points_value = float(points)
            test_passed = points_value >= min_score

            test_payload["status"] = "pass" if test_passed else "fail"
            test_payload["points"] = points_value
            test_payload["passed"] = test_passed
            test_payload["messages"] = messages

            if test_passed:
                passed += 1
            else:
                failed += 1
                all_passed = False
        except Exception as exc:
            test_payload["messages"] = [f"Excepción al ejecutar '{test_id}': {exc}"]
            all_passed = False
            errors += 1

        results[test_id] = test_payload

    base_report["status"] = "ok" if all_passed else "failed"
    base_report["all_passed"] = all_passed
    base_report["results"] = results
    base_report["summary"] = {
        "total": len(tests),
        "passed": passed,
        "failed": failed,
        "errors": errors,
    }
    if not all_passed:
        base_report["errors"].append(
            f"Al menos un test FAIR-EVA no alcanzó el umbral mínimo de {min_score}."
        )

    return base_report


if __name__ == "__main__":
    try:
        if "snakemake" in globals():
            config_file = str(snakemake.input.config)
            upload_response_file = str(snakemake.input.upload_response)
            record_id_file = str(snakemake.input.record_id)
            metadata_xml_file = str(snakemake.input.metadata)
            report_output = Path(str(snakemake.output.report))
            ok_output = Path(str(snakemake.output.ok))
        else:
            if len(sys.argv) < 6:
                print(
                    "Uso: validate_fair_eva.py <config.yaml> <upload_response.json> <record_id.txt> <metadata.xml> <report_output.json> [ok_output.txt]"
                )
                sys.exit(1)
            config_file = sys.argv[1]
            upload_response_file = sys.argv[2]
            record_id_file = sys.argv[3]
            metadata_xml_file = sys.argv[4]
            report_output = Path(sys.argv[5])
            ok_output = Path(sys.argv[6]) if len(sys.argv) > 6 else Path(
                "temp_files/fair_eva_ok.txt"
            )

        report = run_fair_eva_validation(
            config_file,
            upload_response_file,
            record_id_file,
            metadata_xml_file,
        )
        write_json_atomic(report_output, report)

        marker = "OK\n" if report.get("all_passed", False) else "FAIL\n"
        if report.get("status") == "skipped":
            marker = "SKIPPED\n"
        write_text_atomic(ok_output, marker)

        print(f"Reporte FAIR-EVA guardado en: {report_output}")
        print(f"Marcador FAIR-EVA guardado en: {ok_output} ({marker.strip()})")

        if report.get("status") == "failed":
            sys.exit(1)

    except Exception as exc:
        error_report = {
            "timestamp": utc_now_iso(),
            "status": "error",
            "all_passed": False,
            "errors": [str(exc)],
        }
        if "snakemake" in globals():
            report_output = Path(str(snakemake.output.report))
            ok_output = Path(str(snakemake.output.ok))
            write_json_atomic(report_output, error_report)
            write_text_atomic(ok_output, "FAIL\n")
        print(f"Error en validación FAIR-EVA: {exc}", file=sys.stderr)
        sys.exit(1)
