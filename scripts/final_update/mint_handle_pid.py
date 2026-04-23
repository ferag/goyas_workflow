#!/usr/bin/env python3
import base64
import json
import os
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import requests
import yaml


NAMESPACES = {
    "gmd": "http://www.isotc211.org/2005/gmd",
    "gco": "http://www.isotc211.org/2005/gco",
}

DEFAULT_HANDLE_ENDPOINT = "https://hdl.grnet.gr:8001/api/handles"
DEFAULT_HANDLE_PREFIX = "21.T15999"
DEFAULT_PERMISSIONS = "011111110011"
DEFAULT_ADMIN_INDEX = 200
DEFAULT_TARGET_URL_TEMPLATE = "{geonetwork_url}/srv/spa/catalog.search#/metadata/{metadata_uuid}"


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


def load_config(config_file):
    return yaml.safe_load(Path(config_file).read_text(encoding="utf-8")) or {}


def load_session(session_file):
    payload = json.loads(Path(session_file).read_text(encoding="utf-8"))
    xsrf_token = payload.get("XSRF-TOKEN")
    cookies = payload.get("cookies", {})
    if not xsrf_token:
        raise ValueError("No se encontró XSRF-TOKEN en la sesión de GeoNetwork.")
    return xsrf_token, cookies


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


def ensure_geonetwork_url(config):
    geonetwork_url = config.get("services", {}).get("geonetwork", {}).get("url", "")
    geonetwork_url = geonetwork_url.strip().rstrip("/")
    if not geonetwork_url:
        raise ValueError("No se encontró services.geonetwork.url en config.yaml.")
    return geonetwork_url


def _derive_admin_handle(auth_user):
    if not auth_user:
        return ""
    if ":" in auth_user:
        return auth_user.split(":", 1)[1]
    return auth_user


def load_handle_settings(config):
    handle_cfg = config.get("services", {}).get("handle", {})
    if not isinstance(handle_cfg, dict):
        handle_cfg = {}

    endpoint = (
        str(handle_cfg.get("endpoint") or os.getenv("HANDLE_ENDPOINT") or DEFAULT_HANDLE_ENDPOINT)
        .strip()
        .rstrip("/")
    )
    prefix = str(handle_cfg.get("prefix") or os.getenv("HANDLE_PREFIX") or DEFAULT_HANDLE_PREFIX).strip()
    auth_user = str(
        handle_cfg.get("user")
        or os.getenv("HANDLE_USER_ENC")
        or os.getenv("HANDLE_USER_RAW")
        or ""
    ).strip()
    auth_pass = str(handle_cfg.get("password") or os.getenv("HANDLE_PASS") or "").strip()

    admin_handle = str(
        handle_cfg.get("admin_handle")
        or os.getenv("HANDLE_ADMIN_HANDLE")
        or _derive_admin_handle(auth_user)
    ).strip()
    admin_index = int(handle_cfg.get("admin_index", DEFAULT_ADMIN_INDEX))
    permissions = str(handle_cfg.get("permissions", DEFAULT_PERMISSIONS)).strip() or DEFAULT_PERMISSIONS
    target_url_template = str(
        handle_cfg.get("target_url_template", DEFAULT_TARGET_URL_TEMPLATE)
    ).strip()

    if not endpoint:
        raise ValueError("services.handle.endpoint o HANDLE_ENDPOINT es obligatorio.")
    if not prefix:
        raise ValueError("services.handle.prefix o HANDLE_PREFIX es obligatorio.")
    if not auth_user:
        raise ValueError("services.handle.user o HANDLE_USER_ENC/HANDLE_USER_RAW es obligatorio.")
    if not auth_pass:
        raise ValueError("services.handle.password o HANDLE_PASS es obligatorio.")
    if not admin_handle:
        raise ValueError("No se pudo resolver services.handle.admin_handle/HANDLE_ADMIN_HANDLE.")

    return {
        "endpoint": endpoint,
        "prefix": prefix,
        "auth_user": auth_user,
        "auth_pass": auth_pass,
        "admin_handle": admin_handle,
        "admin_index": admin_index,
        "permissions": permissions,
        "target_url_template": target_url_template,
    }


def _basic_auth_header(username, password):
    payload = f"{username}:{password}".encode("utf-8")
    token = base64.b64encode(payload).decode("ascii")
    return f"Basic {token}"


def mint_handle_pid(settings, metadata_uuid, target_url):
    pid = f"{settings['prefix']}/{metadata_uuid}"
    mint_url = f"{settings['endpoint']}/{pid}"

    body = {
        "values": [
            {
                "index": 1,
                "type": "URL",
                "data": {"format": "string", "value": target_url},
            },
            {
                "index": 100,
                "type": "HS_ADMIN",
                "data": {
                    "format": "admin",
                    "value": {
                        "handle": settings["admin_handle"],
                        "index": settings["admin_index"],
                        "permissions": settings["permissions"],
                    },
                },
            },
        ]
    }

    headers = {
        "Authorization": _basic_auth_header(settings["auth_user"], settings["auth_pass"]),
        "Accept": "application/json",
        "Content-Type": "application/json",
    }
    response = requests.put(mint_url, headers=headers, json=body, timeout=60)
    payload = None
    try:
        payload = response.json()
    except ValueError:
        payload = {"raw_response": response.text}

    if response.status_code not in (200, 201):
        raise RuntimeError(
            f"Error minting PID en Handle API ({response.status_code}): {payload}"
        )

    return {
        "pid": pid,
        "mint_url": mint_url,
        "http_status": response.status_code,
        "response": payload,
    }


def register_namespaces(xml_path):
    for _, (prefix, uri) in ET.iterparse(xml_path, events=["start-ns"]):
        try:
            ET.register_namespace(prefix, uri)
        except ValueError:
            pass


def _qn(ns, local_name):
    return f"{{{NAMESPACES[ns]}}}{local_name}"


def set_pid_identifier_in_xml(input_xml, output_xml, pid_value):
    register_namespaces(input_xml)
    tree = ET.parse(input_xml)
    root = tree.getroot()

    md_data_ident = root.find(
        ".//gmd:identificationInfo/gmd:MD_DataIdentification", NAMESPACES
    )
    if md_data_ident is None:
        raise ValueError("No se encontró gmd:MD_DataIdentification en el XML.")

    citation = md_data_ident.find("gmd:citation/gmd:CI_Citation", NAMESPACES)
    if citation is None:
        citation_parent = md_data_ident.find("gmd:citation", NAMESPACES)
        if citation_parent is None:
            citation_parent = ET.SubElement(md_data_ident, _qn("gmd", "citation"))
        citation = ET.SubElement(citation_parent, _qn("gmd", "CI_Citation"))

    for identifier in citation.findall("gmd:identifier", NAMESPACES):
        citation.remove(identifier)

    identifier = ET.Element(_qn("gmd", "identifier"))
    md_identifier = ET.SubElement(identifier, _qn("gmd", "MD_Identifier"))
    code = ET.SubElement(md_identifier, _qn("gmd", "code"))
    code_value = ET.SubElement(code, _qn("gco", "CharacterString"))
    code_value.text = pid_value

    children = list(citation)
    insert_idx = len(children)
    for idx, child in enumerate(children):
        localname = child.tag.split("}", 1)[-1]
        if localname == "edition":
            insert_idx = idx + 1
            break
    citation.insert(insert_idx, identifier)

    tree.write(str(output_xml), encoding="utf-8", xml_declaration=True)


def update_remote_metadata(geonetwork_url, session_file, metadata_xml_path):
    xsrf_token, cookies = load_session(session_file)
    session = requests.Session()
    session.cookies.update(cookies)

    update_url = (
        f"{geonetwork_url}/srv/api/records?"
        "metadataType=METADATA&uuidProcessing=OVERWRITE&transformWith=_none_"
    )
    xml_bytes = Path(metadata_xml_path).read_bytes()
    headers = {
        "Content-Type": "application/xml",
        "Accept": "application/json",
        "X-XSRF-TOKEN": xsrf_token,
    }
    response = session.put(update_url, data=xml_bytes, headers=headers, timeout=90)
    payload = None
    try:
        payload = response.json()
    except ValueError:
        payload = {"raw_response": response.text}

    if response.status_code not in (200, 201):
        raise RuntimeError(
            f"Error actualizando XML en GeoNetwork ({response.status_code}): {payload}"
        )

    return {"update_url": update_url, "http_status": response.status_code, "response": payload}


def run(
    config_file,
    session_file,
    record_id_file,
    upload_response_file,
    metadata_input_file,
    fair_gate_file,
    schema_gate_file,
    metadata_output_file,
):
    config = load_config(config_file)
    enabled = config.get("PID_minting", False) is True
    fair_gate = read_gate_value(fair_gate_file)
    schema_gate = read_gate_value(schema_gate_file)

    report = {
        "timestamp": utc_now_iso(),
        "enabled": enabled,
        "fair_eva_gate": fair_gate,
        "schema_gate": schema_gate,
        "status": "skipped",
        "pid": None,
        "metadata_uuid": None,
        "errors": [],
    }

    input_xml = Path(metadata_input_file)
    output_xml = Path(metadata_output_file)
    output_xml.parent.mkdir(parents=True, exist_ok=True)

    if not enabled:
        output_xml.write_bytes(input_xml.read_bytes())
        report["reason"] = "PID_minting = false"
        return report

    if fair_gate != "OK" or schema_gate != "OK":
        output_xml.write_bytes(input_xml.read_bytes())
        report["reason"] = "Requiere fair_eva_ok=OK y metadata_schema_ok=OK"
        return report

    metadata_uuid = resolve_metadata_uuid(upload_response_file, record_id_file)
    if not metadata_uuid:
        raise ValueError("No se pudo resolver metadata UUID para minting de PID.")
    report["metadata_uuid"] = metadata_uuid

    geonetwork_url = ensure_geonetwork_url(config)
    settings = load_handle_settings(config)
    target_url = settings["target_url_template"].format(
        metadata_uuid=metadata_uuid,
        geonetwork_url=geonetwork_url,
    )

    handle_result = mint_handle_pid(settings, metadata_uuid, target_url)
    report["pid"] = handle_result["pid"]
    report["handle"] = handle_result

    set_pid_identifier_in_xml(input_xml, output_xml, handle_result["pid"])
    remote_update = update_remote_metadata(geonetwork_url, session_file, output_xml)
    report["remote_update"] = remote_update
    report["status"] = "ok"
    return report


if __name__ == "__main__":
    try:
        if "snakemake" in globals():
            config_file = str(snakemake.input.config)
            session_file = str(snakemake.input.session)
            record_id_file = str(snakemake.input.record_id)
            upload_response_file = str(snakemake.input.upload_response)
            metadata_input_file = str(snakemake.input.metadata)
            fair_gate_file = str(snakemake.input.fair_gate)
            schema_gate_file = str(snakemake.input.schema_gate)
            metadata_output_file = str(snakemake.output.metadata)
            report_file = Path(str(snakemake.output.report))
            gate_file = Path(str(snakemake.output.ok))
        else:
            if len(sys.argv) < 11:
                print(
                    "Uso: mint_handle_pid.py "
                    "<config.yaml> <session.json> <record_id.txt> <upload_response.json> "
                    "<metadata_in.xml> <fair_gate.txt> <schema_gate.txt> "
                    "<metadata_out.xml> <report.json> <gate.txt>"
                )
                sys.exit(1)
            config_file = sys.argv[1]
            session_file = sys.argv[2]
            record_id_file = sys.argv[3]
            upload_response_file = sys.argv[4]
            metadata_input_file = sys.argv[5]
            fair_gate_file = sys.argv[6]
            schema_gate_file = sys.argv[7]
            metadata_output_file = sys.argv[8]
            report_file = Path(sys.argv[9])
            gate_file = Path(sys.argv[10])

        result = run(
            config_file=config_file,
            session_file=session_file,
            record_id_file=record_id_file,
            upload_response_file=upload_response_file,
            metadata_input_file=metadata_input_file,
            fair_gate_file=fair_gate_file,
            schema_gate_file=schema_gate_file,
            metadata_output_file=metadata_output_file,
        )
        write_json_atomic(report_file, result)
        marker = "OK\n" if result.get("status") == "ok" else "SKIPPED\n"
        write_text_atomic(gate_file, marker)
        print(f"Mint PID completado con estado: {result.get('status')}")
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
        print(f"Error en mint_handle_pid.py: {exc}", file=sys.stderr)
        sys.exit(1)
