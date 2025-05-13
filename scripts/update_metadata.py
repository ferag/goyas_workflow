 
#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import json
import sys
import requests
import yaml

namespaces = {
    'gmd': 'http://www.isotc211.org/2005/gmd',
    'gmx': 'http://www.isotc211.org/2005/gmx',
    'gco': 'http://www.isotc211.org/2005/gco',
    'gmi': "http://www.isotc211.org/2005/gmi",
    'gml': 'http://www.opengis.net/gml/3.2'
}

def registrar_espacios_de_nombres(xml_file):
    # Itera sobre los eventos de inicio de espacio de nombres en el archivo XML
    for event, (prefix, uri) in ET.iterparse(xml_file, events=['start-ns']):
        try:
            ET.register_namespace(prefix, uri)
        except ValueError as e:
            # Si el prefijo está reservado, se salta
            print(f"Skipping reserved namespace prefix: {prefix} -> {uri} ({e})")

def update_metadata(xml_file, record_id, combined_response, output_file, wms_url, yaml_file, session_file):

    # Parámetros de conexión y configuración
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)

    registrar_espacios_de_nombres(xml_file)

    # Parsear el XML de entrada
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Inserción del bloque onLine dentro de transferOptions en MD_Distribution
    if wms_url != "":
        # Buscar el elemento <gmd:distributionInfo/gmd:MD_Distribution>
        md_distribution = root.find(".//gmd:distributionInfo/gmd:MD_Distribution", namespaces)
        if md_distribution is None:
            # Si no existe, se puede crear un bloque nuevo en root (o según la estructura deseada)
            distribution_info = ET.SubElement(root, f"{{{namespaces['gmd']}}}distributionInfo")
            md_distribution = ET.SubElement(distribution_info, f"{{{namespaces['gmd']}}}MD_Distribution")

        # Crear el elemento <gmd:transferOptions> dentro de MD_Distribution
        transfer_options = ET.SubElement(md_distribution, f"{{{namespaces['gmd']}}}transferOptions")

        # Crear el bloque <gmd:MD_DigitalTransferOptions>
        md_digital_transfer = ET.Element(f"{{{namespaces['gmd']}}}MD_DigitalTransferOptions")

        # Crear el bloque <gmd:onLine>
        online = ET.SubElement(md_digital_transfer, f"{{{namespaces['gmd']}}}onLine")
        ci_online_resource = ET.SubElement(online, f"{{{namespaces['gmd']}}}CI_OnlineResource")

        # Elemento <gmd:linkage> con el <gmd:URL>
        linkage = ET.SubElement(ci_online_resource, f"{{{namespaces['gmd']}}}linkage")
        url_elem = ET.SubElement(linkage, f"{{{namespaces['gmd']}}}URL")
        url_elem.text = wms_url

        # Elemento <gmd:protocol>
        protocol = ET.SubElement(ci_online_resource, f"{{{namespaces['gmd']}}}protocol")
        char_protocol = ET.SubElement(protocol, f"{{{namespaces['gco']}}}CharacterString")
        char_protocol.text = "OGC:WMS"

        # Elemento <gmd:name>
        name_elem = ET.SubElement(ci_online_resource, f"{{{namespaces['gmd']}}}name")
        char_name = ET.SubElement(name_elem, f"{{{namespaces['gco']}}}CharacterString")
        char_name.text = config.get("title")

        # Elemento <gmd:function>
        name_elem = ET.SubElement(ci_online_resource, f"{{{namespaces['gmd']}}}function")
        function_code = ET.SubElement(name_elem, f"{{{namespaces['gmd']}}}CI_OnLineFunctionCode")
        function_code.set('codeList', 'http://standards.iso.org/iso/19139/resources/gmxCodelists.xml#CI_OnLineFunctionCode')
        function_code.set('codeListValue', 'download')


        # Añadir el bloque MD_DigitalTransferOptions a transferOptions
        transfer_options.append(md_digital_transfer)
        print("Bloque onLine insertado dentro de transferOptions en MD_Distribution.")

    # Actualizar <gmd:fileIdentifier>/<gco:CharacterString> con el record_id
    file_identifier = root.find('gmd:fileIdentifier', namespaces)
    if file_identifier is not None:
        char_string = file_identifier.find('gco:CharacterString', namespaces)
        if char_string is not None:
            char_string.text = record_id
        else:
            print("No se encontró <gco:CharacterString> en fileIdentifier")
    else:
        print("No se encontró <gmd:fileIdentifier>")

    # Actualizar recursos en línea existentes para FILE y OGC
    online_resources = root.findall('.//gmd:CI_OnlineResource', namespaces)
    for resource in online_resources:
        resource_url = resource.find('.//gmd:URL', namespaces)
        if resource_url is not None:
            if resource_url.text == 'FILE':
                print(combined_response.get('tif_response', {}).get('url', ''))
                resource_url.text = combined_response.get('tif_response', {}).get('url', '')
                name_elem = resource.find('gmd:name', namespaces)
                if name_elem is not None:
                    char_str = name_elem.find('gco:CharacterString', namespaces)
                    if char_str is not None:
                        char_str.text = combined_response.get('tif_response', {}).get('filename', '')

    # Añadir un nuevo bloque para la miniatura (THUMBNAIL)
    ns_gmd = namespaces['gmd']
    ns_gco = namespaces['gco']

    graphicOverview = ET.Element(f"{{{ns_gmd}}}graphicOverview")
    mdBrowseGraphic = ET.SubElement(graphicOverview, f"{{{ns_gmd}}}MD_BrowseGraphic")
    fileName = ET.SubElement(mdBrowseGraphic, f"{{{ns_gmd}}}fileName")
    charString = ET.SubElement(fileName, f"{{{ns_gco}}}CharacterString")
    # Asigna la URL proveniente de png_response
    charString.text = combined_response.get('png_response', {}).get('url', '')

    # Intentar insertar la miniatura en el MD_DataIdentification
    ci_citation = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification", namespaces)
    if ci_citation is not None:
        ci_citation.append(graphicOverview)
        print("Miniatura añadida a CI_Citation.")
    else:
        # Si no se encuentra, añadirla al final del XML
        root.append(graphicOverview)
        print("No se encontró CI_Citation; miniatura añadida al final del XML.")

      # Guardar el XML actualizado localmente
    tree.write(str(output_file), encoding='utf-8', xml_declaration=True)
    print(f"XML actualizado guardado en: {output_file}")

    # -------------------------------------------------------------------------
    # Guardar el XML modificado
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"XML actualizado y guardado en: {output_file}")

    # Actualizar los metadatos en Geonetwork
    with open(session_file, 'r', encoding='utf-8') as f:
        session_info = json.load(f)
    xsrf_token = session_info.get("XSRF-TOKEN")
    cookies = session_info.get("cookies", {})
    sess = requests.Session()
    sess.cookies.update(cookies)
    headers = {
        'Content-Type': 'application/xml',
        'Accept': 'application/json',
        'X-XSRF-TOKEN': xsrf_token
    }
    update_url = ("https://goyas.csic.es/geonetwork/srv/api/records?"
                  "metadataType=METADATA&uuidProcessing=OVERWRITE&transformWith=_none_")
    updated_xml_bytes = ET.tostring(root, encoding='utf-8', method='xml')
    update_response = sess.put(update_url, data=updated_xml_bytes, headers=headers)
    print("Código de respuesta de actualización remota:", update_response.status_code)
    if update_response.status_code in [200, 201]:
        print("Metadatos actualizados correctamente en Geonetwork.")
        print(update_response.json())
    else:
        print("Error al actualizar metadatos:", update_response.text)
        sys.exit(1)

if __name__ == '__main__':
    try:
        xml_file = snakemake.input.xml              # Ej.: "updated_metadata_with_coverage.xml"
        yaml_file = snakemake.input.config
        record_id_file = snakemake.input.record_id    # Ej.: "metadata_uploaded.txt"
        combined_response_file = snakemake.input.upload_response  # Ej.: "upload_response.json"
        session_file = snakemake.input.session        # Ej.: "geonetwork_session.txt"
        wms = snakemake.input.wms
        output_file = snakemake.output[0]             # Ej.: "final_metadata.xml"
    except NameError:
        if len(sys.argv) < 6:
            print("Uso: update_metadata.py <xml_file> <record_id_file> <combined_response_file> <session_file> <output_file>")
            sys.exit(1)
        xml_file = sys.argv[1]
        yaml_file = sys.argv[2]
        record_id_file = sys.argv[3]
        combined_response_file = sys.argv[4]
        session_file = sys.argv[5]
        wms = sys.argv[6]
        output_file = sys.argv[7]

    with open(record_id_file, 'r', encoding='utf-8') as f:
        record_id = f.read().strip()

    with open(combined_response_file, 'r', encoding='utf-8') as f:
        combined_response = json.load(f)

    with open(wms[0], 'r', encoding='utf-8') as f:
        wms_url = f.readline()

    update_metadata(xml_file, record_id, combined_response, output_file, wms_url, yaml_file, session_file)
