#!/usr/bin/env python3
import calendar
import re
import sys
from datetime import datetime
import xml.etree.ElementTree as ET

import rasterio
from rasterio.warp import transform_bounds
import yaml

namespaces = {
    "gmd": "http://www.isotc211.org/2005/gmd",
    "gmx": "http://www.isotc211.org/2005/gmx",
    "gco": "http://www.isotc211.org/2005/gco",
    "gmi": "http://www.isotc211.org/2005/gmi",
    "gml": "http://www.opengis.net/gml/3.2",
}


def insert_before_first(md_data_id, new_elem, ordered_localnames):
    children = list(md_data_id)
    for idx, child in enumerate(children):
        localname = child.tag.split("}", 1)[-1]
        if localname in ordered_localnames:
            md_data_id.insert(idx, new_elem)
            return
    md_data_id.append(new_elem)


def registrar_espacios_de_nombres(xml_file):
    for event, (prefix, uri) in ET.iterparse(xml_file, events=["start-ns"]):
        try:
            ET.register_namespace(prefix, uri)
        except ValueError as exc:
            print(f"Skipping reserved namespace prefix: {prefix} -> {uri} ({exc})")


def _select_raster_source(dataset_path):
    with rasterio.open(dataset_path) as src:
        if src.count > 0:
            return dataset_path
        subdatasets = list(src.subdatasets or [])

    for subdataset in subdatasets:
        with rasterio.open(subdataset) as subsrc:
            if subsrc.count > 0:
                return subdataset

    raise ValueError(f"No se han encontrado bandas raster en: {dataset_path}")


def _parse_date_like(value, for_end=False):
    if not value:
        return None
    text = str(value).strip()
    if not text:
        return None

    if re.fullmatch(r"\d{4}-\d{2}$", text):
        year, month = map(int, text.split("-"))
        if for_end:
            day = calendar.monthrange(year, month)[1]
        else:
            day = 1
        return f"{year:04d}-{month:02d}-{day:02d}"

    if re.fullmatch(r"\d{4}-\d{2}-\d{2}$", text):
        return text

    iso_candidate = text.replace("Z", "+00:00")
    try:
        dt = datetime.fromisoformat(iso_candidate)
        return dt.date().isoformat()
    except ValueError:
        return None


def _infer_temporal_from_name(dataset_path):
    name = dataset_path.split("/")[-1]
    match_14 = re.search(r"(?<!\d)(\d{14})(?!\d)", name)
    if match_14:
        dt = datetime.strptime(match_14.group(1), "%Y%m%d%H%M%S")
        date_txt = dt.date().isoformat()
        return date_txt, date_txt

    match_8 = re.search(r"(?<!\d)(\d{8})(?!\d)", name)
    if match_8:
        dt = datetime.strptime(match_8.group(1), "%Y%m%d")
        date_txt = dt.date().isoformat()
        return date_txt, date_txt

    return None, None


def _infer_temporal_from_tags(dataset_path):
    tag_candidates = [
        "time_coverage_start",
        "time_coverage_end",
        "NC_GLOBAL#time_coverage_start",
        "NC_GLOBAL#time_coverage_end",
        "NC_GLOBAL#isodate",
        "isodate",
    ]

    begin = None
    end = None

    with rasterio.open(dataset_path) as src:
        tags = src.tags()
        begin = _parse_date_like(tags.get("time_coverage_start") or tags.get("NC_GLOBAL#time_coverage_start"))
        end = _parse_date_like(
            tags.get("time_coverage_end") or tags.get("NC_GLOBAL#time_coverage_end"), for_end=True
        )
        if begin and end:
            return begin, end

        for key in tag_candidates:
            parsed = _parse_date_like(tags.get(key), for_end=True)
            if parsed:
                return parsed, parsed

    return None, None


def _resolve_temporal_coverage(config, dataset_path):
    temporal_cfg = config.get("metadata", {}).get("coverage", {}).get("temporal", {}) or {}
    temporal_auto = temporal_cfg.get("auto", False) is True

    if temporal_auto:
        begin, end = _infer_temporal_from_name(dataset_path)
        if begin and end:
            print(f"Cobertura temporal inferida del nombre de archivo: {begin} -> {end}")
            return begin, end

        begin, end = _infer_temporal_from_tags(dataset_path)
        if begin and end:
            print(f"Cobertura temporal inferida de metadatos del raster: {begin} -> {end}")
            return begin, end

        print("No se pudo inferir la cobertura temporal automáticamente. Se usará la configuración manual.")

    begin = _parse_date_like(temporal_cfg.get("begin"))
    end = _parse_date_like(temporal_cfg.get("end"), for_end=True)
    begin = begin or "1900-01-01"
    end = end or begin
    return begin, end


def _resolve_spatial_coverage(config, dataset_path):
    spatial_cfg = config.get("metadata", {}).get("coverage", {}).get("spatial", {}) or {}
    if not spatial_cfg.get("auto", False):
        return (
            float(spatial_cfg.get("west", -111.4564)),
            float(spatial_cfg.get("south", 39.1428)),
            float(spatial_cfg.get("east", -102.4564)),
            float(spatial_cfg.get("north", 46.1428)),
        )

    source_path = _select_raster_source(dataset_path)
    with rasterio.open(source_path) as src:
        bounds = src.bounds
        src_crs = src.crs
        if src_crs is None:
            epsg = config.get("crs", {}).get("epsg")
            if epsg:
                src_crs = f"EPSG:{epsg}"

        if src_crs and str(src_crs).upper() != "EPSG:4326":
            west, south, east, north = transform_bounds(src_crs, "EPSG:4326", *bounds, densify_pts=21)
        else:
            west, south, east, north = bounds.left, bounds.bottom, bounds.right, bounds.top

    print("Bounding Box calculado:")
    print(f"  west={west} east={east} south={south} north={north}")
    return west, south, east, north


def _upsert_temporal_extent(root, data_ident, begin_date, end_date):
    time_period = root.find(".//gml:TimePeriod", namespaces)
    if time_period is None:
        temporal_extent = ET.Element(f"{{{namespaces['gmd']}}}extent")
        ex_extent = ET.SubElement(temporal_extent, f"{{{namespaces['gmd']}}}EX_Extent")
        temporal_element = ET.SubElement(ex_extent, f"{{{namespaces['gmd']}}}temporalElement")
        ex_temporal_extent = ET.SubElement(temporal_element, f"{{{namespaces['gmd']}}}EX_TemporalExtent")
        extent_elem = ET.SubElement(ex_temporal_extent, f"{{{namespaces['gmd']}}}extent")
        time_period = ET.SubElement(
            extent_elem, f"{{{namespaces['gml']}}}TimePeriod", attrib={f"{{{namespaces['gml']}}}id": "tp1"}
        )

        if data_ident is not None:
            insert_before_first(data_ident, temporal_extent, ["supplementalInformation"])
        else:
            root.append(temporal_extent)

    begin_elem = time_period.find("gml:beginPosition", namespaces)
    if begin_elem is None:
        begin_elem = ET.SubElement(time_period, f"{{{namespaces['gml']}}}beginPosition")
    begin_elem.text = begin_date

    end_elem = time_period.find("gml:endPosition", namespaces)
    if end_elem is None:
        end_elem = ET.SubElement(time_period, f"{{{namespaces['gml']}}}endPosition")
    end_elem.text = end_date


def add_coverage_to_xml(xml_file, config, dataset_file, output_file):
    registrar_espacios_de_nombres(xml_file)
    tree = ET.parse(xml_file)
    root = tree.getroot()

    data_ident = root.find(".//gmd:identificationInfo/gmd:MD_DataIdentification", namespaces)

    begin_date, end_date = _resolve_temporal_coverage(config, dataset_file)
    _upsert_temporal_extent(root, data_ident, begin_date, end_date)

    west, south, east, north = _resolve_spatial_coverage(config, dataset_file)

    spatial_extent = ET.Element(f"{{{namespaces['gmd']}}}extent")
    ex_extent_spatial = ET.SubElement(spatial_extent, f"{{{namespaces['gmd']}}}EX_Extent")
    geographic_element = ET.SubElement(ex_extent_spatial, f"{{{namespaces['gmd']}}}geographicElement")
    bbox = ET.SubElement(geographic_element, f"{{{namespaces['gmd']}}}EX_GeographicBoundingBox")

    west_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}westBoundLongitude")
    ET.SubElement(west_elem, f"{{{namespaces['gco']}}}Decimal").text = str(west)

    east_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}eastBoundLongitude")
    ET.SubElement(east_elem, f"{{{namespaces['gco']}}}Decimal").text = str(east)

    south_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}southBoundLatitude")
    ET.SubElement(south_elem, f"{{{namespaces['gco']}}}Decimal").text = str(south)

    north_elem = ET.SubElement(bbox, f"{{{namespaces['gmd']}}}northBoundLatitude")
    ET.SubElement(north_elem, f"{{{namespaces['gco']}}}Decimal").text = str(north)

    if data_ident is not None:
        insert_before_first(data_ident, spatial_extent, ["supplementalInformation"])
        print("Cobertura espacial añadida dentro de MD_DataIdentification.")
    else:
        root.append(spatial_extent)
        print("Cobertura espacial añadida al final del XML.")

    tree.write(str(output_file), encoding="utf-8", xml_declaration=True)
    print(f"XML con cobertura actualizado guardado en: {output_file}")


if __name__ == "__main__":
    try:
        xml_file = snakemake.input.xml
        yaml_file = snakemake.input.config
        dataset_file = snakemake.input.data
        output_file = snakemake.output[0]
    except NameError:
        if len(sys.argv) < 5:
            print("Uso: add_coverage.py <xml_file> <yaml_file> <dataset_file> <output_file>")
            sys.exit(1)
        xml_file = sys.argv[1]
        yaml_file = sys.argv[2]
        dataset_file = sys.argv[3]
        output_file = sys.argv[4]

    with open(yaml_file, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    add_coverage_to_xml(xml_file, config, dataset_file, output_file)
