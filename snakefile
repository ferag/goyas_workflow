# GOYAS Workflow - Snakefile
# Fase 1: Validacion 
# Fase 2: Generar xml base
# Fase 3: Carga de datos
# Fase 4: Actualización final XML
#############################

configfile: "config.yaml"

RUN_ALL = config.get("run_all", False)

### Fase 1
##########

# Paso 1.1: Verificación del entorno (Python + dependencias + directorios)
rule check_environment:
    input:
        "requirements.txt"
    output:
        "logs/env_check.json"
    script:
        "scripts/validation/check_environment.py"

# Paso 1.2: Validación del archivo config.yaml
rule validate_config:
    input:
        "config.yaml"
    output:
        "logs/config_validation.json"
    script:
        "scripts/validation/validate_config.py"

# Ejecuta todos los pasos de validación
rule validation:
    input:
        rules.check_environment.output,
        rules.validate_config.output

### Fase 2
##########
        

# Paso 2.1: Actualización del XML en base al YAML
rule update_xml:
    input:
        base_xml="base_iso19139.xml",
        config="config.yaml"
    output:
        "temp_files/updated_metadata.xml"
    script:
        "scripts/generate_xml/update_xml.py"

# Paso 2.2: Publicación en Geoserver (si está activada en el YAML)
rule publish_geoserver:
    input:
        data=config["dataset"]["file"],
        config="config.yaml"
    output:
        "temp_files/geoserver_publish_done.txt"
    script:
        "scripts/generate_xml/publish_geoserver.py"

# Paso 2.3: Login en Geonetwork
rule geonetwork_login:
    input:
        config="config.yaml"
    output:
        session="temp_files/geonetwork_session.txt"
    script:
        "scripts/generate_xml/geonetwork_login.py"

# Paso 2.4: Subida de la primera versión de metadatos
rule upload_metadata_initial:
    input:
        session="temp_files/geonetwork_session.txt",
        metadata="temp_files/updated_metadata.xml"
    output:
        "temp_files/metadata_uploaded.txt"
    script:
        "scripts/generate_xml/upload_metadata_initial.py"

# Paso 2.5: Incluir coverage
rule add_coverage:
    input:
        xml="temp_files/updated_metadata.xml",
        config="config.yaml",
        data=config['dataset']['file']
    output:
        "temp_files/updated_metadata_with_coverage.xml"
    script:
        "scripts/generate_xml/add_coverage.py"

# Paso 2.6: Actualiza content info
rule add_contentinfo:
    input:
        xml="temp_files/updated_metadata_with_coverage.xml",
        config="config.yaml",
        data=config['dataset']['file']
    output:
        "temp_files/updated_metadata_with_content.xml"
    script:
        "scripts/generate_xml/add_contentinfo.py"

# Paso 2.7: Genera miniatura
rule tif_snapshot:
    input:
        data=config['dataset']['file']
    output:
        "temp_files/snapshot.png"
    script:
        "scripts/generate_xml/tif_to_png.py"

# Ejecuta todos los pasos de validación
rule generate_xml:
    input:
        rules.update_xml.output,
        rules.publish_geoserver.output,
        rules.geonetwork_login.output,
        rules.upload_metadata_initial.output,
        rules.add_coverage.output,
        rules.add_contentinfo.output,
        rules.tif_snapshot.output,

### Fase 3
##########

rule upload_data:
    input:
        session="temp_files/geonetwork_session.txt",
        record_id="temp_files/metadata_uploaded.txt",
        config="config.yaml",
        png="temp_files/snapshot.png"

    output:
        "temp_files/upload_response.json"
    script:
        "scripts/upload_data/upload_data.py"

### Fase 4
##########

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
        "scripts/final_update/update_metadata.py"

rule validate_metadata_geonetwork:
    input:
        config="config.yaml",
        session="temp_files/geonetwork_session.txt",
        record_id="temp_files/metadata_uploaded.txt",
        metadata=rules.update_metadata.output
    output:
        "logs/metadata_validation.json"
    script:
        "scripts/validation/validate_metadata_geonetwork.py"



# Regla para ejecutar todos los pasos y permitir reanudación desde el último correcto
rule run_all:
    input:
        rules.validation.output,
        rules.validate_metadata_geonetwork.output
    output:
        "logs/run_all.done"
    shell:
        "touch {output}"

# Regla principal que engloba todos los pasos
rule all:
    input:
        "temp_files/final_metadata.xml",
        "logs/metadata_validation.json",
