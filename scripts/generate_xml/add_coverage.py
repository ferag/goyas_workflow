 
#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import sys
import yaml
import os
import rasterio

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

def add_coverage_to_xml(xml_file, config, tif_file, output_file):
    # Registrar namespaces (lee el XML para capturar los namespaces definidos)
    registrar_espacios_de_nombres(xml_file)
    tree = ET.parse(xml_file)
    root = tree.getroot()
    registrar_espacios_de_nombres(xml_file)


    # Buscar el elemento MD_DataIdentification dentro de identificationInfo
    data_ident = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification", namespaces)
    if data_ident is None:
        print("No se encontró el elemento gmd:MD_DataIdentification. Se agregará al final del documento.")

    # --- Bloque de cobertura temporal ---
    temporal = config.get("metadata").get('coverage', {}).get('temporal', {})
    begin_date = temporal.get('begin', '1900-01-01')
    end_date = temporal.get('end', '1900-01-01')

    # Buscar si ya existe un elemento gml:TimePeriod
    time_period = root.find('.//gml:TimePeriod', namespaces)

    if time_period is None:
        # Si no existe, crearlo con un id arbitrario
        time_period = ET.Element(f"{{{namespaces['gml']}}}TimePeriod", attrib={f"{{{namespaces['gml']}}}id": "tp1"})
        begin_elem = ET.SubElement(time_period, f"{{{namespaces['gml']}}}beginPosition")
        begin_elem.text = begin_date
        end_elem = ET.SubElement(time_period, f"{{{namespaces['gml']}}}endPosition")
        end_elem.text = end_date

        # Aquí debes insertar el time_period en su sitio correcto dentro del árbol,
        # por ejemplo, dentro de temporalElement si fuera parte del gmd:EX_TemporalExtent
        # Esto depende de la estructura de tu XML.
        # Ejemplo:
        # parent = root.find('.//gmd:extent', namespaces)
        # if parent is not None:
        #     parent.append(time_period)
    else:
        # Si ya existe, modificar los valores
        begin_elem = time_period.find('gml:beginPosition', namespaces)
        if begin_elem is not None:
            begin_elem.text = begin_date
        else:
            begin_elem = ET.SubElement(time_period, f"{{{namespaces['gml']}}}beginPosition")
            begin_elem.text = begin_date

        end_elem = time_period.find('gml:endPosition', namespaces)
        if end_elem is not None:
            end_elem.text = end_date
        else:
            end_elem = ET.SubElement(time_period, f"{{{namespaces['gml']}}}endPosition")
            end_elem.text = end_date

   # Crear el bloque de cobertura espacial
    spatial = config.get("metadata").get('coverage', {}).get('spatial', {})
    if spatial.get('auto', False):
        # Abre el archivo GeoTIFF
        print("Calculating bounding box")
        with rasterio.open(tif_file) as src:
            bounds = src.bounds
            west = bounds.left
            south = bounds.bottom
            east = bounds.right
            north = bounds.top

        print("Bounding Box:")
        print(f"west = {west}")
        print(f"east = {east}")
        print(f"south = {south}")
        print(f"north = {north}")
    else:
        west = spatial.get('west', -111.4564)
        east = spatial.get('east', -102.4564)
        south = spatial.get('south', 39.1428)
        north = spatial.get('north', 46.1428)

    # Crear la estructura de extensión espacial
    spatial_extent = ET.Element(f"{{{namespaces['gmd']}}}extent")
    ex_extent_spatial = ET.SubElement(spatial_extent, f"{{{namespaces['gmd']}}}EX_Extent")
    geographic_element = ET.SubElement(ex_extent_spatial, f"{{{namespaces['gmd']}}}geographicElement")
    bbox = ET.SubElement(geographic_element, f"{{{namespaces['gmd']}}}EX_GeographicBoundingBox")

    west_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}westBoundLongitude")
    west_decimal = ET.SubElement(west_elem, f"{{{namespaces['gco']}}}Decimal")
    west_decimal.text = str(west)

    east_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}eastBoundLongitude")
    east_decimal = ET.SubElement(east_elem, f"{{{namespaces['gco']}}}Decimal")
    east_decimal.text = str(east)

    south_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}southBoundLatitude")
    south_decimal = ET.SubElement(south_elem, f"{{{namespaces['gco']}}}Decimal")
    south_decimal.text = str(south)

    north_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}northBoundLatitude")
    north_decimal = ET.SubElement(north_elem, f"{{{namespaces['gco']}}}Decimal")
    north_decimal.text = str(north)

    # --- Insertar los bloques de cobertura ---
    if data_ident is not None:
        # Si se encontró MD_DataIdentification, se agregan allí
        data_ident.append(spatial_extent)
        print("Cobertura añadida dentro de MD_DataIdentification.")
    else:
        # Si no se encuentra, se agregan al final del documento
        root.append(spatial_extent)
        print("Cobertura añadida al final del XML, ya que no se encontró MD_DataIdentification.")

    # Guardar el XML actualizado
    tree.write(str(output_file), encoding='utf-8', xml_declaration=True)
    print(f"XML con cobertura actualizado guardado en: {output_file}")


if __name__ == '__main__':
    try:
        # Acceder a inputs y outputs vía snakemake
        xml_file = snakemake.input.xml              # Ej.: "updated_metadata.xml"
        yaml_file = snakemake.input.config           # Ej.: "config.yaml"
        tif_file = snakemake.input.data
        output_file = snakemake.output[0]             # Ej.: "updated_metadata_with_coverage.xml"
    except NameError:
        # Ejecución manual
        if len(sys.argv) < 4:
            print("Uso: add_coverage.py <xml_file> <yaml_file> <output_file>")
            sys.exit(1)
        xml_file = sys.argv[1]
        yaml_file = sys.argv[2]
        tif_file = sys.argv[3]
        output_file = sys.argv[4]

    # Leer el YAML de configuración
    with open(yaml_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    add_coverage_to_xml(xml_file, config, tif_file, output_file)
