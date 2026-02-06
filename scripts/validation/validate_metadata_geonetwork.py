#!/usr/bin/env python3
import json
import sys
from pathlib import Path

import requests
import yaml


def load_session(session_file):
    with open(session_file, 'r', encoding='utf-8') as f:
        session_info = json.load(f)
    xsrf_token = session_info.get("XSRF-TOKEN")
    cookies = session_info.get("cookies", {})
    if not xsrf_token:
        raise ValueError("No se encontró el token XSRF en el archivo de sesión.")
    return xsrf_token, cookies


def read_record_id(record_id_file):
    with open(record_id_file, 'r', encoding='utf-8') as f:
        record_id = f.read().strip()
    if not record_id:
        raise ValueError("El archivo de record_id está vacío.")
    return record_id


def load_geonetwork_url(config_file):
    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    geonetwork_url = config.get("services", {}).get("geonetwork", {}).get("url")
    if not geonetwork_url:
        raise ValueError("No se encontró services.geonetwork.url en el config.")
    return geonetwork_url.rstrip('/')


def validate_metadata(config_file, session_file, record_id_file, output_file):
    geonetwork_url = load_geonetwork_url(config_file)
    xsrf_token, cookies = load_session(session_file)
    record_id = read_record_id(record_id_file)

    session = requests.Session()
    session.cookies.update(cookies)
    headers = {
        'Accept': 'application/json',
        'X-XSRF-TOKEN': xsrf_token,
    }

    validate_url = f"{geonetwork_url}/srv/api/records/{record_id}/validate"
    response = session.post(validate_url, headers=headers)

    if response.status_code not in [200, 201]:
        print(f"Error al validar metadatos: {response.status_code} - {response.text}")
        sys.exit(1)

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        payload = response.json()
    except ValueError:
        payload = {"raw_response": response.text}

    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding='utf-8')
    print(f"Validación completada. Resultado guardado en: {output_file}")


if __name__ == "__main__":
    try:
        config_file = snakemake.input.config
        session_file = snakemake.input.session
        record_id_file = snakemake.input.record_id
        output_file = snakemake.output[0]
    except NameError:
        if len(sys.argv) < 5:
            print("Uso: validate_metadata_geonetwork.py <config_file> <session_file> <record_id_file> <output_file>")
            sys.exit(1)
        config_file = sys.argv[1]
        session_file = sys.argv[2]
        record_id_file = sys.argv[3]
        output_file = sys.argv[4]

    validate_metadata(config_file, session_file, record_id_file, output_file)
