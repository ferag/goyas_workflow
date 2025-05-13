#!/usr/bin/env python3
import requests
import sys
import json

def upload_data(session_file, record_id_file, tif_file_path, png_file_path, output_file):
    geonetwork_url = "https://goyas.csic.es/geonetwork"

    # Cargar la información de la sesión (token y cookies)
    with open(session_file, 'r', encoding='utf-8') as f:
        session_info = json.load(f)
    xsrf_token = session_info.get("XSRF-TOKEN")
    cookies = session_info.get("cookies", {})

    # Crear sesión y actualizar cookies
    session = requests.Session()
    session.cookies.update(cookies)
    headers = {'X-XSRF-TOKEN': xsrf_token, 'Accept': 'application/json'}

    # Leer el record_id generado en el paso de metadatos
    with open(record_id_file, 'r', encoding='utf-8') as f:
        record_id = f.read().strip()

    # URL para subir archivos asociados al record_id
    upload_url = f"{geonetwork_url}/srv/api/records/{record_id}/attachments"

    # Subir el archivo TIF
    with open(tif_file_path, 'rb') as tif_file:
        files = {'file': tif_file}
        tif_response = session.post(upload_url, files=files, headers=headers)

    if tif_response.status_code in [200, 201]:
        print("Archivo TIF asociado correctamente")
    else:
        print(f"Error al asociar TIF: {tif_response.text}")
        sys.exit(1)

    # Subir el archivo PNG (miniatura)
    with open(png_file_path, 'rb') as png_file:
        files_png = {'file': png_file}
        png_response = session.post(upload_url, files=files_png, headers=headers)

    if png_response.status_code in [200, 201]:
        print("PNG thumbnail subido correctamente")
    else:
        print(f"Error al asociar PNG: {png_response.text}")
        sys.exit(1)

    # Guardar las respuestas combinadas en un archivo JSON
    combined_response = {
        "tif_response": tif_response.json(),
        "png_response": png_response.json()
    }
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(json.dumps(combined_response, indent=4))
    print(f"Respuestas de upload guardadas en: {output_file}")

if __name__ == "__main__":
    try:
        # Acceder a inputs y outputs mediante snakemake
        session_file = snakemake.input.session
        record_id_file = snakemake.input.record_id
        tif_file_path = snakemake.input.data
        png_file_path = snakemake.input.png
        output_file = snakemake.output[0]
    except NameError:
        # Modo manual
        if len(sys.argv) < 6:
            print("Uso: upload_data.py <session_file> <record_id_file> <tif_file_path> <png_file_path> <output_file>")
            sys.exit(1)
        session_file = sys.argv[1]
        record_id_file = sys.argv[2]
        tif_file_path = sys.argv[3]
        png_file_path = sys.argv[4]
        output_file = sys.argv[5]

    upload_data(session_file, record_id_file, tif_file_path, png_file_path, output_file)
