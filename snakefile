configfile: "config.yaml"

# Paso 1: Actualización del XML en base al YAML
rule update_xml:
    input:
        base_xml="base_iso19139.xml",
        config="config.yaml"
    output:
        "temp_files/updated_metadata.xml"
    script:
        "scripts/update_xml.py"

# Paso 2: Publicación en Geoserver (si está activada en el YAML)
rule publish_geoserver:
    input:
        data=config["file"],
        config="config.yaml"
    output:
        "temp_files/geoserver_publish_done.txt"
    script:
        "scripts/publish_geoserver.py"

# Paso 3: Login en Geonetwork
rule geonetwork_login:
    input:
        config="config.yaml"
    output:
        session="temp_files/geonetwork_session.txt"
    script:
        "scripts/geonetwork_login.py"

# Paso 4: Subida de la primera versión de metadatos
rule upload_metadata_initial:
    input:
        session="temp_files/geonetwork_session.txt",
        metadata="temp_files/updated_metadata.xml"
    output:
        "temp_files/metadata_uploaded.txt"
    script:
        "scripts/upload_metadata_initial.py"

rule add_coverage:
    input:
        xml="temp_files/updated_metadata.xml",
        config="config.yaml",
        data=config['file']
    output:
        "temp_files/updated_metadata_with_coverage.xml"
    script:
        "scripts/add_coverage.py"

rule add_contentinfo:
    input:
        xml="temp_files/updated_metadata_with_coverage.xml",
        config="config.yaml",
        data=config['file']
    output:
        "temp_files/updated_metadata_with_content.xml"
    script:
        "scripts/add_contentinfo.py"

rule tif_snapshot:
    input:
        data=config['file']
    output:
        "temp_files/snapshot.png"
    script:
        "scripts/tif_to_png.py"

rule upload_data:
    input:
        session="temp_files/geonetwork_session.txt",
        record_id="temp_files/metadata_uploaded.txt",
        config="config.yaml",
        png="temp_files/snapshot.png"

    output:
        "temp_files/upload_response.json"
    script:
        "scripts/upload_data.py"

rule update_metadata:
    input:
        xml="temp_files/updated_metadata_with_content.xml",
        config="config.yaml",
        record_id="temp_files/metadata_uploaded.txt",
        upload_response="temp_files/upload_response.json",
        wms=rules.publish_geoserver.output,
        session="temp_files/geonetwork_session.txt"
    output:
        "temp_files/final_metadata.xml"
    script:
        "scripts/update_metadata.py"



# Regla principal que engloba todos los pasos
rule all:
    input:
        rules.update_xml.output,
