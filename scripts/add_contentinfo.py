 
#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import sys
import yaml
import os
import rasterio
import xml.dom.minidom as md


# Definir namespaces (agregar gml además de los ya conocidos)
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


def insert_after(target_elem, new_elem, parent):
    """
    Inserta new_elem en parent inmediatamente después de target_elem,
    usando xml.etree.ElementTree estándar (sin getparent()).
    """
    children = list(parent)
    for i, child in enumerate(children):
        if child is target_elem:
            parent.insert(i + 1, new_elem)
            return
    # Si no se encontró el target_elem entre hijos, hacer append al final
    parent.append(new_elem)


def update_contentinfo(root, parameters_cfg):
    """
    Inserta un bloque <gmd:contentInfo> -> <gmi:MI_CoverageDescription> con
    rangeElementDescription para cada parámetro (name, definition, unit).
    Sin error/precision (para no invalidar la ISO en gmi:MI_Band).
    """
    if not parameters_cfg:
        return None  # No hay contentInfo si no hay parámetros

    # 1) Crear <gmd:contentInfo> -> <gmi:MI_CoverageDescription>
    gmd_contentinfo = ET.Element(f"{{{namespaces['gmd']}}}contentInfo")
    coverage_desc = ET.SubElement(gmd_contentinfo, f"{{{namespaces['gmi']}}}MI_CoverageDescription")

    # Bloques fijos de attributeDescription y contentType
    attr_desc = ET.SubElement(coverage_desc, f"{{{namespaces['gmd']}}}attributeDescription")
    record_type = ET.SubElement(attr_desc, f"{{{namespaces['gco']}}}RecordType")
    record_type.text = "Chl"  # O lo que quieras como genérico

    content_type = ET.SubElement(coverage_desc, f"{{{namespaces['gmd']}}}contentType")
    code_el = ET.SubElement(content_type, f"{{{namespaces['gmd']}}}MD_CoverageContentTypeCode", {
        "codeListValue": "physicalMeasurement",
        "codeList": "http://standards.iso.org/iso/19139/resources/gmxCodelists.xml#MD_CoverageContentTypeCode"
    })

    # 2) Para cada parámetro, creamos un <gmi:rangeElementDescription>
    for param_key, param_vals in parameters_cfg.items():
        name_val = param_vals.get("name", param_key)
        definition_val = param_vals.get("definition", "")
        unit_val = param_vals.get("unit", "unitless")
        type_val = param_vals.get("type")

        range_elem_desc = ET.SubElement(coverage_desc, f"{{{namespaces['gmi']}}}rangeElementDescription")
        mi_red = ET.SubElement(range_elem_desc, f"{{{namespaces['gmi']}}}MI_RangeElementDescription")

        # <gmi:name>
        gmi_name = ET.SubElement(mi_red, f"{{{namespaces['gmi']}}}name")
        char_name = ET.SubElement(gmi_name, f"{{{namespaces['gco']}}}CharacterString")
        char_name.text = name_val

        # <gmi:definition>
        gmi_def = ET.SubElement(mi_red, f"{{{namespaces['gmi']}}}definition")
        char_def = ET.SubElement(gmi_def, f"{{{namespaces['gco']}}}CharacterString")
        char_def.text = definition_val

        # <gmi:rangeElement> -> <gmi:MI_Band>
        gmi_range_elem = ET.SubElement(mi_red, f"{{{namespaces['gmi']}}}rangeElement")
        record = ET.SubElement(gmi_range_elem, f"{{{namespaces['gco']}}}Record")
        mi_band = ET.SubElement(record, f"{{{namespaces['gmi']}}}MI_Band")

        # sequenceIdentifier
        seq_id = ET.SubElement(mi_band, f"{{{namespaces['gmd']}}}sequenceIdentifier")
        member_name = ET.SubElement(seq_id, f"{{{namespaces['gco']}}}MemberName")

        a_name = ET.SubElement(member_name, f"{{{namespaces['gco']}}}aName")
        a_name_char = ET.SubElement(a_name, f"{{{namespaces['gco']}}}CharacterString")
        a_name_char.text = name_val

        attr_type = ET.SubElement(member_name, f"{{{namespaces['gco']}}}attributeType")
        type_name = ET.SubElement(attr_type, f"{{{namespaces['gco']}}}TypeName")
        tname_aaname = ET.SubElement(type_name, f"{{{namespaces['gco']}}}aName")
        tname_char = ET.SubElement(tname_aaname, f"{{{namespaces['gco']}}}CharacterString")
        tname_char.text = type_val

        # descriptor
        descriptor = ET.SubElement(mi_band, f"{{{namespaces['gmd']}}}descriptor")
        descriptor_char = ET.SubElement(descriptor, f"{{{namespaces['gco']}}}CharacterString")
        descriptor_char.text = definition_val

        # units (sin el slash en gml:id)
        valid_id = unit_val.replace("/", "_").replace(" ", "_")
        gmd_units = ET.SubElement(mi_band, f"{{{namespaces['gmd']}}}units")
        unit_def = ET.SubElement(gmd_units, f"{{{namespaces['gml']}}}UnitDefinition", {
            f"{{{namespaces['gml']}}}id": f"id_{valid_id}"
        })
        gml_identifier = ET.SubElement(unit_def, f"{{{namespaces['gml']}}}identifier", {
            "codeSpace": "http://www.opengis.net/def/uom/OGC/1.0"
        })
        gml_identifier.text = unit_val
        gml_name = ET.SubElement(unit_def, f"{{{namespaces['gml']}}}name")
        gml_name.text = unit_val
        gml_qtype = ET.SubElement(unit_def, f"{{{namespaces['gml']}}}quantityType")
        gml_qtype.text = type_val

    return gmd_contentinfo


    # 3) Insertar <gmd:contentInfo> justo después de <gmd:identificationInfo>
    #    Buscamos la etiqueta <gmd:identificationInfo> y la insertamos después
    identification_info = root.find(".//gmd:identificationInfo", namespaces)
    if identification_info is not None and parameters_cfg:
        insert_after(identification_info, content_info, root)

    # Si no se encontró identificationInfo o no hay parámetros, podrías appendeárselo a root
    elif parameters_cfg:
        root.append(content_info)

def create_scope_element(param_name):
    """
    Crea el bloque:
      <gmd:scope>
        <gmd:DQ_Scope>
          <gmd:level> ... </gmd:level>
          <gmd:levelDescription>
            <gmd:MD_ScopeDescription>
              <gmd:other>
                <gco:CharacterString>param_name</gco:CharacterString>
              </gmd:other>
            </gmd:MD_ScopeDescription>
          </gmd:levelDescription>
        </gmd:DQ_Scope>
      </gmd:scope>
    """
    gmd = namespaces["gmd"]
    gco = namespaces["gco"]

    scope_el = ET.Element(f"{{{gmd}}}scope")
    dq_scope = ET.SubElement(scope_el, f"{{{gmd}}}DQ_Scope")

    level_el = ET.SubElement(dq_scope, f"{{{gmd}}}level")
    scope_code = ET.SubElement(level_el, f"{{{gmd}}}MD_ScopeCode", {
        "codeList": "http://standards.iso.org/iso/19139/resources/gmxCodelists.xml#MD_ScopeCode",
        "codeListValue": "attribute"
    })

    # levelDescription con MD_ScopeDescription
    lvl_desc = ET.SubElement(dq_scope, f"{{{gmd}}}levelDescription")
    md_scope_desc = ET.SubElement(lvl_desc, f"{{{gmd}}}MD_ScopeDescription")
    other = ET.SubElement(md_scope_desc, f"{{{gmd}}}other")
    char_str = ET.SubElement(other, f"{{{gco}}}CharacterString")
    char_str.text = param_name

    return scope_el

def create_quantitative_result_block(parent_elem, id_str, val_str, name_val, type_val):
    """
    Crea un <gmd:report> con un <gmd:DQ_QuantitativeAttributeAccuracy> que
    contiene un <gmd:result> -> <gmd:DQ_QuantitativeResult>.
    """
    gmd = namespaces["gmd"]
    gco = namespaces["gco"]
    gml = namespaces["gml"]

    report_el = ET.SubElement(parent_elem, f"{{{gmd}}}report")
    dq_attr = ET.SubElement(report_el, f"{{{gmd}}}DQ_QuantitativeAttributeAccuracy")
    dq_result = ET.SubElement(dq_attr, f"{{{gmd}}}result")
    dq_qresult = ET.SubElement(dq_result, f"{{{gmd}}}DQ_QuantitativeResult")

    # valueUnit
    val_unit = ET.SubElement(dq_qresult, f"{{{gmd}}}valueUnit")
    unit_def = ET.SubElement(val_unit, f"{{{gml}}}UnitDefinition", {
        f"{{{gml}}}id": id_str + "_" + name_val
    })
    gml_id = ET.SubElement(unit_def, f"{{{gml}}}identifier", {
        "codeSpace": "http://www.opengis.net/def/uom/OGC/1.0"
    })
    gml_id.text = name_val
    gml_name = ET.SubElement(unit_def, f"{{{gml}}}name")
    gml_name.text = id_str
    gml_qtype = ET.SubElement(unit_def, f"{{{gml}}}quantityType")
    gml_qtype.text = type_val

    # value
    val_el = ET.SubElement(dq_qresult, f"{{{gmd}}}value")
    rec_el = ET.SubElement(val_el, f"{{{gco}}}Record")
    rec_el.text = val_str

def add_data_quality(root, config):
    """
    - Parsear el XML base y el YAML con 'parameters'.
    - Para cada parámetro con 'error' o 'precision', se crea <gmd:dataQualityInfo> con:
        <gmd:DQ_DataQuality>
          <gmd:scope> param_name </gmd:scope>
          <gmd:report> (uno para error) </gmd:report>
          <gmd:report> (otro para precision) </gmd:report>
    - Se ubica cada dataQualityInfo tras <gmd:identificationInfo>, si existe, o se appendea al root.
    """

    parameters = config.get("parameters", {})
    if not parameters:
        print("No se encontró 'parameters' en YAML. No se añaden dataQualityInfo.")
        tree.write(output_file, encoding="utf-8", xml_declaration=True)
        return

    # Ubicación donde insertamos dataQuality
    distributionInfo = root.find(".//{http://www.isotc211.org/2005/gmd}distributionInfo")
    parent = root

    for param_key, param_vals in parameters.items():
        error_val = param_vals.get("error")
        precision_val = param_vals.get("precision")
        name_val = param_vals.get("name")
        type_val = param_vals.get("type")

        if error_val is None and precision_val is None:
            continue

        param_name = param_vals.get("name", param_key)

        # gmd:dataQualityInfo -> gmd:DQ_DataQuality
        dq_info = ET.Element(f"{{{namespaces['gmd']}}}dataQualityInfo")
        dq_data_quality = ET.SubElement(dq_info, f"{{{namespaces['gmd']}}}DQ_DataQuality")

        # scope
        scope_el = create_scope_element(param_name)
        dq_data_quality.append(scope_el)

        # Crear un <gmd:report> para error
        if error_val is not None:
            create_quantitative_result_block(dq_data_quality, "ErrorValue", str(error_val), name_val, type_val)

        # Crear un <gmd:report> para precision
        if precision_val is not None:
            create_quantitative_result_block(dq_data_quality, "PrecisionValue", str(precision_val), name_val, type_val)

        # Insertar
        if distributionInfo is not None:
            insert_after(distributionInfo, dq_info, parent)
        else:
            parent.append(dq_info)


def main(xml_file, config, output_file):
    # Parsear XML y YAML
    registrar_espacios_de_nombres(xml_file)
    ET.register_namespace("gmi", "http://www.isotc211.org/2005/gmi")
    tree = ET.parse(xml_file)
    root = tree.getroot()

    parameters_cfg = config.get("parameters", {})

    # 1) Insertar contentInfo
    gmd_contentinfo = update_contentinfo(root, parameters_cfg)
    # Insertar <gmd:contentInfo> tras <gmd:identificationInfo>, si corresponde
    if gmd_contentinfo is not None:
        identification_info = root.find(".//{http://www.isotc211.org/2005/gmd}identificationInfo")
        if identification_info is not None:
            insert_after(identification_info, gmd_contentinfo, root)
        else:
            root.append(gmd_contentinfo)

    # 2) Insertar dataQualityInfo para error/precision
    add_data_quality(root, config)

    # Guardar
    tree.write(output_file, encoding='utf-8', xml_declaration=True)

    # Utilizar minidom para obtener un XML con sangrado (pretty print)
    with open(str(output_file), 'r', encoding='utf-8') as f:
        xml_str = f.read()

    dom = md.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent="  ")

    # Eliminar líneas en blanco innecesarias
    pretty_xml = "\n".join([line for line in pretty_xml.splitlines() if line.strip()])

    with open(str(output_file), 'w', encoding='utf-8') as f:
        f.write(pretty_xml)

    print(f"XML actualizado con contentInfo. Guardado en {output_file}")



if __name__ == '__main__':
    try:
        # Acceder a inputs y outputs vía snakemake
        xml_file = snakemake.input.xml              # Ej.: "updated_metadata.xml"
        yaml_file = snakemake.input.config           # Ej.: "config.yaml"
        output_file = snakemake.output[0]             # Ej.: "updated_metadata_with_coverage.xml"
    except NameError:
        # Ejecución manual
        if len(sys.argv) < 4:
            print("Uso: add_coverage.py <xml_file> <yaml_file> <output_file>")
            sys.exit(1)
        xml_file = sys.argv[1]
        yaml_file = sys.argv[2]
        output_file = sys.argv[3]

    # Leer el YAML de configuración
    with open(yaml_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    main(xml_file, config, output_file)
