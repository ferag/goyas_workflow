#!/usr/bin/env python3
import requests
import sys
import json
import sys
import boto3
from botocore.client import Config as BotoConfig
from boto3.s3.transfer import TransferConfig, S3Transfer
import yaml
import os

def upload_data(session_file, record_id_file, config, png_file_path, output_file):
    geonetwork_url = config.get("geonetwork")['url']
    local_file = config['file']
    endpoint_url = config.get("storage")['endpoint_url']
    access_key = config.get("storage")['access_key']
    secret_key = config.get("storage")['secret_key']
    bucket_name = config.get("storage")['bucket']

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

    # Subir el archivo PNG (miniatura)
    with open(png_file_path, 'rb') as png_file:
        files_png = {'file': png_file}
        png_response = session.post(upload_url, files=files_png, headers=headers)

    if png_response.status_code in [200, 201]:
        print("PNG thumbnail subido correctamente")
    else:
        print(f"Error al asociar PNG: {png_response.text}")
        sys.exit(1)

    # === CREAR CLIENTE S3 ===
    s3 = boto3.client(
        's3',
        endpoint_url=endpoint_url,
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        config=BotoConfig(signature_version='s3'),
    )

    # Configuración para multipart upload
    tfr_cfg = TransferConfig(
        multipart_threshold=50 * 1024 * 1024,
        multipart_chunksize=25 * 1024 * 1024,
        use_threads=True
    )

    # Callback de progreso
    class ProgressPercentage:
        def __init__(self, filename):
            self._filename = filename
            self._size = float(os.path.getsize(filename))
            self._seen = 0
        def __call__(self, bytes_amount):
            self._seen += bytes_amount
            pct = (self._seen / self._size) * 100
            sys.stdout.write(f"\r{self._filename} — {self._seen:.0f}/{self._size:.0f} bytes ({pct:.2f}%)")
            sys.stdout.flush()

    # === SUBIDA A S3 ===
    transfer = S3Transfer(client=s3, config=tfr_cfg)
    transfer.upload_file(
        local_file,
        bucket_name,
        local_file,
        callback=ProgressPercentage(local_file),
        extra_args={'ACL': 'public-read'}
    )
    print("\n✅ Subida de TIFF completada.")

    # Generar URL firmada (1 hora)
    signed_url = s3.generate_presigned_url(
        'get_object',
        Params={'Bucket': bucket_name, 'Key': local_file},
        ExpiresIn=3600
    )
    print(f"🔐 URL firmada (1 h): {signed_url}")

    # Guardar las respuestas combinadas en un archivo JSON
    combined_response = {
        "tif_response": {"url": signed_url, "filename": local_file},
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
        yaml_file = snakemake.input.config
        png_file_path = snakemake.input.png
        output_file = snakemake.output[0]
    except NameError:
        # Modo manual
        if len(sys.argv) < 6:
            print("Uso: upload_data.py <session_file> <record_id_file> <tif_file_path> <png_file_path> <output_file>")
            sys.exit(1)
        session_file = sys.argv[1]
        record_id_file = sys.argv[2]
        yaml_file = sys.argv[3]
        png_file_path = sys.argv[4]
        output_file = sys.argv[4]
    # Leer el YAML de configuración
    with open(yaml_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    upload_data(session_file, record_id_file, config, png_file_path, output_file)
