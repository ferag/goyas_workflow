#!/usr/bin/env python3
import os
import requests
import sys
import json
import yaml
from datetime import datetime, timezone

# Si se ejecuta dentro de Snakemake, la variable 'snakemake' estará definida.
try:
    yaml_file = snakemake.input.config
    session_file = snakemake.output.session
except NameError:
    if len(sys.argv) < 2:
        print("Uso: geonetwork_login.py <output_file>")
        sys.exit(1)
    yaml_file = sys.argv[1]
    session_file = sys.argv[2]

# Leer el YAML de configuración
with open(yaml_file, 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

def login_geonetwork(config, output_file):
    # URL base y credenciales
    geonetwork_url = config.get("services")['geonetwork']['url']
    username = config.get("services")['geonetwork']['username']
    password = config.get("services")['geonetwork']['password']

    # Crear sesión
    session = requests.Session()

    # Autenticación inicial para capturar cookies y obtener el token XSRF
    auth_response = session.post(f"{geonetwork_url}/srv/me")
    xsrf_token = session.cookies.get("XSRF-TOKEN")
    if not xsrf_token:
        # Fallback por compatibilidad con instancias que solo exponen la cabecera.
        if "Set-Cookie" in auth_response.headers:
            cookie_value = auth_response.headers["Set-Cookie"]
            if "XSRF-TOKEN=" in cookie_value:
                xsrf_token = cookie_value.split("XSRF-TOKEN=")[-1].split(";")[0]
    if not xsrf_token:
        print("Error: no se pudo obtener XSRF-TOKEN en la respuesta de /srv/me")
        sys.exit(1)
    print("Token XSRF obtenido.")

    # Preparar cabeceras con el token XSRF
    headers = {'X-XSRF-TOKEN': xsrf_token, 'Accept': 'application/json'}

    # Autenticación con credenciales en la URL de login
    auth_response = session.post(f"{geonetwork_url}/srv/eng/me", auth=(username, password), headers=headers)
    if auth_response.status_code == 200:
        print("Login inicial en /srv/eng/me completado.")
    else:
        print(f"Error de autenticación: {auth_response.text}")
        sys.exit(1)

    # Verifica que la sesión API esté realmente autenticada.
    me_response = session.get(f"{geonetwork_url}/srv/api/me", headers=headers)
    if me_response.status_code != 200:
        print(
            f"Error de autenticación API: /srv/api/me devolvió {me_response.status_code} - {me_response.text}"
        )
        sys.exit(1)

    try:
        me_payload = me_response.json()
    except ValueError:
        print("Error: /srv/api/me no devolvió un JSON válido.")
        sys.exit(1)

    api_username = str(me_payload.get("username", "")).strip()
    if not api_username or api_username.lower() == "anonymous":
        print("Error: la sesión API no está autenticada (usuario anónimo o vacío).")
        sys.exit(1)

    print(f"Sesión API autenticada como: {api_username}")

    # Guardar la información de la sesión (token y cookies) para usar en pasos posteriores
    session_info = {
        "XSRF-TOKEN": xsrf_token,
        "cookies": requests.utils.dict_from_cookiejar(session.cookies),
        "authenticated_user": me_payload,
        "authenticated_at": datetime.now(timezone.utc).isoformat(),
        "geonetwork_url": geonetwork_url,
    }

    # Asegurarse de que el directorio del archivo output_file exista
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(session_info, f)

    print(f"Información de la sesión guardada en: {output_file}")

if __name__ == "__main__":
    login_geonetwork(config, session_file)
