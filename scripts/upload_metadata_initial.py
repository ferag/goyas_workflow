#!/usr/bin/env python3
import requests
import sys
import json

def upload_metadata(session_file, xml_file_path, output_file):
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

    # Construir URL para subir los metadatos
    upload_url = f"{geonetwork_url}/srv/api/records?metadataType=METADATA&uuidProcessing=GENERATEUUID&transformWith=_none_"

    # Leer y subir el fichero XML de metadatos
    with open(xml_file_path, 'rb') as xml_file:
        files = {'file': xml_file}
        metadata_response = session.post(upload_url, files=files, headers=headers)

    if metadata_response.status_code in [200, 201]:
        print("Metadatos subidos correctamente")
        response_json = metadata_response.json()
        record_id = None
        # Extraer el record_id de la respuesta (se asume que está en 'metadataInfos')
        for key in response_json.get('metadataInfos', {}):
            record_list = response_json['metadataInfos'][key]
            if record_list and isinstance(record_list, list):
                record_id = record_list[0].get('uuid')
                break
        if record_id is None:
            print("No se pudo extraer el record_id del JSON de respuesta.")
            sys.exit(1)
        else:
            print("UUID: %s" % record_id)
        # Guardar el record_id en el fichero de salida para el siguiente paso
        with open(str(output_file), 'w', encoding='utf-8') as f:
            f.write(record_id)
    else:
        print(f"Error al subir metadatos: {metadata_response.text}")
        sys.exit(1)

if __name__ == "__main__":
    try:
        # Acceso a los inputs y outputs mediante el objeto snakemake
        session_file = snakemake.input.session
        xml_file_path = snakemake.input.metadata
        output_file = snakemake.output
        print(session_file)
    except NameError:
        # Modo de ejecución manual
        if len(sys.argv) < 4:
            print("Uso: upload_metadata_initial.py <session_file> <xml_file_path> <output_file>")
            sys.exit(1)
            session_file = sys.argv[1]
            xml_file_path = sys.argv[2]
            output_file = sys.argv[3]

    upload_metadata(session_file, xml_file_path, output_file)
