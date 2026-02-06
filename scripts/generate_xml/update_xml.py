#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import yaml
import xml.dom.minidom as md
from datetime import datetime
from zoneinfo import ZoneInfo

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

def update_xml_with_processing(xml_file, yaml_file, output_file):
    # Registrar los namespaces para la salida
    registrar_espacios_de_nombres(xml_file)

    # Cargar la configuración YAML
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)

    # Parsear el XML de entrada
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Actualizar el título en la sección identificationInfo, si se define en el YAML
    yaml_title = config.get("metadata").get("title")
    if yaml_title:
        # Buscar el elemento <gco:CharacterString> dentro de identificationInfo > MD_DataIdentification > citation > CI_Citation > title
        title_elem = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString",
                              {"gmd": "http://www.isotc211.org/2005/gmd",
                               "gco": "http://www.isotc211.org/2005/gco"})
        if title_elem is not None:
            title_elem.text = yaml_title
            print(f"Título actualizado a: {yaml_title}")
        else:
            print("No se encontró el elemento title en identificationInfo.")

    # Actualizar la fecha de publicación en la sección identificationInfo, si se define en el YAML

    date_elem = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:DateTime",
                            {"gmd": "http://www.isotc211.org/2005/gmd",
                            "gco": "http://www.isotc211.org/2005/gco"})
    if date_elem is not None:
        now = datetime.now(ZoneInfo("Europe/Madrid"))
        formatted = now.isoformat()
        date_elem.text = formatted
        print(f"Título actualizado a: {formatted}")
    else:
        print("No se encontró el elemento title en identificationInfo.")

    # Actualizar el abstract en la sección identificationInfo, si se define en el YAML
    yaml_abstract = config.get("metadata").get("abstract")
    if yaml_title:
        # Buscar el elemento <gco:CharacterString> dentro de identificationInfo > MD_DataIdentification > citation > CI_Citation > title
        abstract_elem = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:abstract/gco:CharacterString",
                              {"gmd": "http://www.isotc211.org/2005/gmd",
                               "gco": "http://www.isotc211.org/2005/gco"})
        if abstract_elem is not None:
            abstract_elem.text = yaml_abstract
            print(f"Título actualizado a: {yaml_abstract}")
        else:
            print("No se encontró el elemento title en identificationInfo.")

    # -------------------------------------------------------------------------
    # (1) INFORMACIÓN DE CONTACTO
    # YAML:
    # contact:
    #   individualName: "Fernando Aguilar"
    #   organisationName: "IFCA-CSIC"
    #   email: "aguilarf@ifca.unican.es"
    # -------------------------------------------------------------------------
    contact_cfg = config.get("metadata").get("contact", {})
    if contact_cfg:
        # Buscar el bloque <gmd:contact>/<gmd:CI_ResponsibleParty>
        ci_responsible = root.find(".//gmd:contact/gmd:CI_ResponsibleParty", namespaces)
        if ci_responsible is not None:
            # individualName
            indiv_elem = ci_responsible.find("gmd:individualName/gco:CharacterString", namespaces)
            if indiv_elem is not None:
                indiv_elem.text = contact_cfg.get("individualName", "")

            # organisationName
            org_elem = ci_responsible.find("gmd:organisationName/gco:CharacterString", namespaces)
            if org_elem is not None:
                org_elem.text = contact_cfg.get("organisationName", "")

            # email
            email_elem = ci_responsible.find(
                "./gmd:contactInfo/gmd:CI_Contact/gmd:address/gmd:CI_Address/"
                "gmd:electronicMailAddress/gco:CharacterString",
                namespaces
            )
            if email_elem is not None:
                email_elem.text = contact_cfg.get("email", "")
        else:
            print("No se encontró gmd:CI_ResponsibleParty dentro de gmd:contact. No se actualizará la sección de contacto.")

        # Buscar el bloque <gmd:contact>/<gmd:CI_ResponsibleParty>
        poc_responsible = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:pointOfContact/gmd:CI_ResponsibleParty", namespaces)
        if poc_responsible is not None:
            # individualName
            indiv_elem = poc_responsible.find("gmd:individualName/gco:CharacterString", namespaces)
            if indiv_elem is not None:
                indiv_elem.text = contact_cfg.get("individualName", "")

            # organisationName
            org_elem = poc_responsible.find("gmd:organisationName/gco:CharacterString", namespaces)
            if org_elem is not None:
                org_elem.text = contact_cfg.get("organisationName", "")

            # positionName
            indiv_elem = poc_responsible.find("gmd:positionName/gco:CharacterString", namespaces)
            if indiv_elem is not None:
                indiv_elem.text = contact_cfg.get("role", "")

            # email
            email_elem = poc_responsible.find(
                "./gmd:contactInfo/gmd:CI_Contact/gmd:address/gmd:CI_Address/"
                "gmd:electronicMailAddress/gco:CharacterString",
                namespaces
            )
            if email_elem is not None:
                email_elem.text = contact_cfg.get("email", "")
        else:
            print("No se encontró gmd:CI_ResponsibleParty dentro de gmd:contact. No se actualizará la sección de contacto.")



    # -------------------------------------------------------------------------
    # (2) RESOLUCIÓN ESPACIAL
    # YAML:
    # resolution:
    #   x: 10
    #   y: 10
    # -------------------------------------------------------------------------
    resolution_cfg = config.get("metadata").get("resolution", {})
    if resolution_cfg:
        x_res = resolution_cfg.get("x", None)
        y_res = resolution_cfg.get("y", None)

        # Buscar <gmd:MD_GridSpatialRepresentation>
        md_grid = root.find(".//gmd:MD_GridSpatialRepresentation", namespaces)
        if md_grid is not None and x_res is not None and y_res is not None:
            # Ejes row/column -> resoluciones
            axis_dims = md_grid.findall("gmd:axisDimensionProperties/gmd:MD_Dimension", namespaces)
            # Se asume que el primero corresponde a row, el segundo a column
            if len(axis_dims) >= 2:
                # row
                row_res_elem = axis_dims[0].find("gmd:resolution/gco:Measure", namespaces)
                if row_res_elem is not None:
                    row_res_elem.text = str(x_res)

                # column
                col_res_elem = axis_dims[1].find("gmd:resolution/gco:Measure", namespaces)
                if col_res_elem is not None:
                    col_res_elem.text = str(y_res)
        else:
            print("No se encontró gmd:MD_GridSpatialRepresentation o faltan x/y. No se actualiza resolución.")

    # -------------------------------------------------------------------------
    # (3) PALABRAS CLAVE (keywords)
    # YAML:
    # keywords:
    #   theme: "Water; NDWI"
    #   place: "Cuerda del Pozo; Soria"
    # -------------------------------------------------------------------------
    keywords_cfg = config.get("metadata").get("keywords", {})
    if keywords_cfg:
        theme_keywords = keywords_cfg.get("theme", "")
        place_keywords = keywords_cfg.get("place", "")

        
        # 3.1 Bloque "theme"
        if theme_keywords:
            # Crear <gmd:descriptiveKeywords> -> <gmd:MD_Keywords>
            theme_desc = ET.Element(f"{{{namespaces['gmd']}}}descriptiveKeywords")
            theme_md_kw = ET.SubElement(theme_desc, f"{{{namespaces['gmd']}}}MD_Keywords")
            for kw in theme_keywords:
                kw_elem = ET.SubElement(theme_md_kw, f"{{{namespaces['gmd']}}}keyword")
                char_kw = ET.SubElement(kw_elem, f"{{{namespaces['gco']}}}CharacterString")
                char_kw.text = kw
            # gmd:type
            kw_type = ET.SubElement(theme_md_kw, f"{{{namespaces['gmd']}}}type")
            type_code = ET.SubElement(kw_type, f"{{{namespaces['gmd']}}}MD_KeywordTypeCode", {
                "codeListValue": "theme",
                "codeList": "http://standards.iso.org/iso/19139/resources/gmxCodelists.xml#MD_KeywordTypeCode"
            })
            # Añadir al root o a la sección identificationInfo
            # Normalmente, se añade dentro de identificationInfo/gmd:MD_DataIdentification
            md_data_id = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification", namespaces)
            if md_data_id is not None:
                md_data_id.append(theme_desc)
            else:
                root.append(theme_desc)

        # 3.2 Bloque "place"
        if place_keywords:
            place_desc = ET.Element(f"{{{namespaces['gmd']}}}descriptiveKeywords")
            place_md_kw = ET.SubElement(place_desc, f"{{{namespaces['gmd']}}}MD_Keywords")
            for kw in place_keywords:
                kw_elem = ET.SubElement(place_md_kw, f"{{{namespaces['gmd']}}}keyword")
                char_kw = ET.SubElement(kw_elem, f"{{{namespaces['gco']}}}CharacterString")
                char_kw.text = kw
            # gmd:type
            kw_type = ET.SubElement(place_md_kw, f"{{{namespaces['gmd']}}}type")
            type_code = ET.SubElement(kw_type, f"{{{namespaces['gmd']}}}MD_KeywordTypeCode", {
                "codeListValue": "place",
                "codeList": "http://standards.iso.org/iso/19139/resources/gmxCodelists.xml#MD_KeywordTypeCode"
            })
            md_data_id = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification", namespaces)
            if md_data_id is not None:
                md_data_id.append(place_desc)
            else:
                root.append(place_desc)

    # -------------------------------------------------------------------------
    # (4) LICENCIA
    # YAML:
    # license:
    #   type: "Creative Commons 4.0"
    # -------------------------------------------------------------------------
    license_cfg = config.get("metadata").get("license", {})
    if license_cfg:
        lic_type = license_cfg.get("type", "")
        # Buscar o crear gmd:resourceConstraints/gmd:MD_LegalConstraints
        md_legal = root.find(".//gmd:resourceConstraints/gmd:MD_LegalConstraints", namespaces)
        if md_legal is not None and lic_type.strip() != "":
            # Actualizar <gmd:otherConstraints>
            other_constr = md_legal.find("gmd:otherConstraints/gco:CharacterString", namespaces)
            if other_constr is not None:
                other_constr.text = lic_type
            else:
                # Crear si no existe
                other_constr_block = md_legal.find("gmd:otherConstraints", namespaces)
                if not other_constr_block:
                    other_constr_block = ET.SubElement(md_legal, f"{{{namespaces['gmd']}}}otherConstraints", {"{http://www.isotc211.org/2005/gco}nilReason": "missing"})
                char_str = ET.SubElement(other_constr_block, f"{{{namespaces['gco']}}}CharacterString")
                char_str.text = lic_type
        else:
            print("No se encontró gmd:MD_LegalConstraints o 'license.type' está vacío, no se actualiza la licencia.")

    # Guardar el XML modificado en el fichero de salida
    try:
        tree.write(str(output_file), encoding="utf-8", xml_declaration=True)
    except Exception as e:
        print(e)

    # Utilizar minidom para obtener un XML con sangrado (pretty print)
    with open(str(output_file), 'r', encoding='utf-8') as f:
        xml_str = f.read()

    dom = md.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent="  ")

    # Eliminar líneas en blanco innecesarias
    pretty_xml = "\n".join([line for line in pretty_xml.splitlines() if line.strip()])

    with open(str(output_file), 'w', encoding='utf-8') as f:
        f.write(pretty_xml)

    print(f"XML actualizado y formateado guardado en: {output_file}")

if __name__ == "__main__":
    # Los archivos se reciben a través del objeto snakemake
    xml_file = snakemake.input.base_xml       # XML de entrada
    yaml_file = snakemake.input.config          # YAML con configuración (incluyendo title y processing)
    output_file = snakemake.output             # XML de salida modificado
    update_xml_with_processing(xml_file, yaml_file, output_file)
