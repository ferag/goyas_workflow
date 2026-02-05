#!/usr/bin/env python3
import os
import requests
import sys
import json
import yaml

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
    if 'Set-Cookie' not in auth_response.headers:
        print("Error: No se recibió la cabecera Set-Cookie en la respuesta de /srv/me")
        sys.exit(1)
    # Extraer el token XSRF de la cabecera. (Ajusta la extracción si es necesario.)
    xsrf_token = auth_response.headers['Set-Cookie'][11:].split(';')[0]
    print("Token obtenido:", xsrf_token)

    # Preparar cabeceras con el token XSRF
    headers = {'X-XSRF-TOKEN': xsrf_token, 'Accept': 'application/json'}

    # Autenticación con credenciales en la URL de login
    auth_response = session.post(f"{geonetwork_url}/srv/eng/me", auth=(username, password), headers=headers)
    if auth_response.status_code == 200:
        print("Autenticación exitosa")
    else:
        print(f"Error de autenticación: {auth_response.text}")
        sys.exit(1)

    # Guardar la información de la sesión (token y cookies) para usar en pasos posteriores
    session_info = {
        "XSRF-TOKEN": xsrf_token,
        "cookies": requests.utils.dict_from_cookiejar(session.cookies)
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
