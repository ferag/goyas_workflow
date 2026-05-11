#!/usr/bin/env python3
import json
import re
import sys
import unicodedata
from pathlib import Path
from urllib.parse import urlencode

import rasterio
from rasterio.crs import CRS
from rasterio.transform import Affine
from rasterio.warp import transform_bounds
import requests
import yaml


SUCCESS_CODES = {200, 201}


def eliminar_acentos(texto):
    forma_nfkd = unicodedata.normalize("NFKD", str(texto or ""))
    return "".join(c for c in forma_nfkd if not unicodedata.combining(c))


def slugify_name(text, fallback="coverage", max_len=48):
    normalized = eliminar_acentos(text).lower().strip()
    normalized = re.sub(r"[^a-z0-9_]+", "_", normalized)
    normalized = re.sub(r"_+", "_", normalized).strip("_")
    if not normalized:
        normalized = fallback
    return normalized[:max_len]


def _select_raster_source(dataset_path):
    with rasterio.open(dataset_path) as src:
        if src.count > 0:
            return dataset_path
        subdatasets = list(src.subdatasets or [])

    for subdataset in subdatasets:
        with rasterio.open(subdataset) as src:
            if src.count > 0:
                return subdataset

    raise ValueError(f"No se han encontrado bandas raster en: {dataset_path}")


def _select_netcdf_subdataset(dataset_path, config):
    dataset_cfg = config.get("dataset", {}) or {}
    geoserver_cfg = (config.get("services", {}) or {}).get("geoserver", {}) or {}
    metadata_cfg = config.get("metadata", {}) or {}
    preferred = (
        geoserver_cfg.get("netcdf_variable")
        or dataset_cfg.get("variable")
        or dataset_cfg.get("netcdf_variable")
    )
    if not preferred:
        parameters = metadata_cfg.get("parameters") or []
        if isinstance(parameters, list) and parameters:
            preferred = parameters[0].get("key") or parameters[0].get("name")

    with rasterio.open(dataset_path) as src:
        if src.count > 0:
            return dataset_path
        subdatasets = list(src.subdatasets or [])

    if preferred:
        preferred = str(preferred).strip()
        for subdataset in subdatasets:
            if subdataset.rsplit("//", 1)[-1] == preferred:
                return subdataset

    for subdataset in subdatasets:
        with rasterio.open(subdataset) as src:
            if src.count > 0:
                return subdataset

    raise ValueError(f"No se han encontrado variables raster publicables en: {dataset_path}")


def _netcdf_global_georef(dataset_path, config):
    with rasterio.open(dataset_path) as src:
        tags = src.tags()

    transform = None
    geotransform = tags.get("spatial_ref_GeoTransform")
    if geotransform:
        values = [float(value) for value in str(geotransform).split()]
        if len(values) == 6:
            transform = Affine.from_gdal(*values)

    crs = None
    epsg = (config.get("crs", {}) or {}).get("epsg")
    if epsg:
        crs = CRS.from_epsg(int(epsg))
    elif tags.get("proj4_string"):
        crs = CRS.from_string(tags["proj4_string"])

    return transform, crs


def _materialize_netcdf_as_geotiff(dataset_path, config, output_file):
    subdataset = _select_netcdf_subdataset(dataset_path, config)
    fallback_transform, fallback_crs = _netcdf_global_georef(dataset_path, config)

    with rasterio.open(subdataset) as src:
        profile = src.profile.copy()
        profile.update(
            driver="GTiff",
            compress="deflate",
            BIGTIFF="IF_SAFER",
        )
        if src.crs is None and fallback_crs is not None:
            profile["crs"] = fallback_crs
        if src.transform.is_identity and fallback_transform is not None:
            profile["transform"] = fallback_transform

        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.unlink(missing_ok=True)
        with rasterio.open(output_file, "w", **profile) as dst:
            for band_index in range(1, src.count + 1):
                dst.write(src.read(band_index), band_index)

    print(f"NetCDF materializado como GeoTIFF para GeoServer: {subdataset} -> {output_file}")
    return str(output_file)


def _format_bbox_value(value):
    return f"{float(value):.6f}".rstrip("0").rstrip(".")


def _resolve_wms_map_params(config, dataset_path, width=768, preview_width=2040):
    source_path = _select_raster_source(dataset_path)
    with rasterio.open(source_path) as src:
        bounds = src.bounds
        src_crs = src.crs
        if src_crs is None:
            epsg = (config.get("crs", {}) or {}).get("epsg")
            if epsg:
                src_crs = f"EPSG:{epsg}"

        if src_crs and str(src_crs).upper() != "EPSG:4326":
            west, south, east, north = transform_bounds(
                src_crs, "EPSG:4326", *bounds, densify_pts=21
            )
        else:
            west, south, east, north = bounds.left, bounds.bottom, bounds.right, bounds.top

        bbox_width = abs(east - west)
        bbox_height = abs(north - south)

        if bbox_width > 0 and bbox_height > 0:
            height = max(round(width * bbox_height / bbox_width), 1)
            preview_height = max(round(preview_width * bbox_height / bbox_width), 1)
        else:
            height = max(round(width * max(src.height, 1) / max(src.width, 1)), 1)
            preview_height = max(round(preview_width * max(src.height, 1) / max(src.width, 1)), 1)

    return {
        "bbox": ",".join(_format_bbox_value(value) for value in (west, south, east, north)),
        "preview_width": str(preview_width),
        "preview_height": str(preview_height),
        "width": str(width),
        "height": str(height),
    }


def _build_wms_layer_url(geoserver_url, workspace, layer_name, wms_params):
    query = urlencode(
        [
            ("SERVICE", "WMS"),
            ("SERVICE", "WMS"),
            ("VERSION", "1.3.0"),
            ("REQUEST", "GetCapabilities"),
            ("FORMAT", "image/png"),
            ("TRANSPARENT", "true"),
            ("LAYERS", layer_name),
            ("STYLES", ""),
            ("CRS", "EPSG:3857"),
            ("WIDTH", wms_params["preview_width"]),
            ("HEIGHT", wms_params["preview_height"]),
            ("BBOX", wms_params["bbox"]),
            ("width", wms_params["width"]),
            ("height", wms_params["height"]),
            ("srs", "EPSG:4326"),
        ]
    )
    return f"{geoserver_url.rstrip('/')}/{workspace}/ows?{query}"


def _write_wms_output(output_file, geoserver_url, workspace, layer_name, config, dataset_path):
    wms_params = _resolve_wms_map_params(config, dataset_path)
    wms_url = _build_wms_layer_url(geoserver_url, workspace, layer_name, wms_params)
    with open(str(output_file), "w", encoding="utf-8") as handle:
        handle.write(wms_url)
    print(f"URL WMS de la capa {workspace}:{layer_name} guardada en: {output_file}")


def _style_layer(geoserver_url, workspace, layer_name, style, auth):
    if not style:
        return
    url = f"{geoserver_url}/rest/layers/{workspace}:{layer_name}"
    payload = {"layer": {"defaultStyle": {"name": style}}}
    headers = {"accept": "application/json", "content-type": "application/json"}
    response = requests.put(url, headers=headers, data=json.dumps(payload), auth=auth, timeout=60)
    if response.status_code not in SUCCESS_CODES:
        print(
            f"Aviso: no se pudo asignar estilo '{style}' a {workspace}:{layer_name} "
            f"({response.status_code}): {response.text}"
        )


def _extract_coverages(payload):
    if not isinstance(payload, dict):
        return []

    coverages_block = payload.get("coverages", {})
    if isinstance(coverages_block, str) or coverages_block is None:
        return []
    if isinstance(coverages_block, list):
        raw_items = coverages_block
    elif isinstance(coverages_block, dict):
        raw_items = coverages_block.get("coverage", [])
    else:
        raw_items = []
    if isinstance(raw_items, dict):
        raw_items = [raw_items]

    names = []
    for item in raw_items:
        if isinstance(item, dict):
            name = item.get("name")
        else:
            name = str(item)
        if name:
            names.append(name.split(":")[-1])
    return names


def _create_geotiff_store(geoserver_url, workspace, dataset_path, store_name, auth):
    upload_url = (
        f"{geoserver_url}/rest/workspaces/{workspace}/coveragestores/{store_name}/file.geotiff"
        f"?configure=first&coverageName={store_name}"
    )
    headers = {"Content-Type": "image/tiff", "accept": "application/json"}
    with open(dataset_path, "rb") as handle:
        response = requests.put(upload_url, headers=headers, data=handle, auth=auth, timeout=300)

    if response.status_code not in SUCCESS_CODES:
        raise RuntimeError(
            f"Error al subir GeoTIFF al store '{store_name}' ({response.status_code}): {response.text}"
        )
    print(f"Store GeoTIFF publicado: {store_name}")
    return store_name


def _create_netcdf_store(geoserver_url, workspace, netcdf_url, base_store_name, auth):
    file_uri = str(netcdf_url)

    for suffix in range(0, 20):
        store_name = base_store_name if suffix == 0 else f"{base_store_name}_{suffix}"
        create_url = f"{geoserver_url}/rest/workspaces/{workspace}/coveragestores"
        payload = {
            "coverageStore": {
                "name": store_name,
                "type": "NetCDF",
                "enabled": True,
                "workspace": {"name": workspace},
                "url": file_uri,
            }
        }
        headers = {"accept": "application/json", "content-type": "application/json"}
        response = requests.post(create_url, headers=headers, data=json.dumps(payload), auth=auth, timeout=120)

        if response.status_code in SUCCESS_CODES:
            print(f"Store NetCDF creado: {store_name} -> {file_uri}")
            return store_name

        response_text = (response.text or "").lower()
        if response.status_code == 409 or "already exists" in response_text:
            continue

        raise RuntimeError(
            "Error al crear store NetCDF "
            f"'{store_name}' ({response.status_code}): {response.text}"
        )

    raise RuntimeError("No se pudo crear un nombre único para el store NetCDF.")


def _list_store_coverages(geoserver_url, workspace, store_name, auth):
    url = f"{geoserver_url}/rest/workspaces/{workspace}/coveragestores/{store_name}/coverages.json"
    headers = {"accept": "application/json"}
    response = requests.get(url, headers=headers, auth=auth, timeout=120)
    if response.status_code not in SUCCESS_CODES:
        raise RuntimeError(
            f"No se pudieron listar coberturas de '{store_name}' ({response.status_code}): {response.text}"
        )

    try:
        payload = response.json()
    except ValueError:
        payload = {}
    names = _extract_coverages(payload)
    if not names:
        print(f"Aviso: coverages vacíos para store '{store_name}'. Payload={payload}")
    return names


def _upload_netcdf_store(geoserver_url, workspace, dataset_path, store_name, auth):
    upload_url = (
        f"{geoserver_url}/rest/workspaces/{workspace}/coveragestores/{store_name}/file.netcdf"
        "?configure=all"
    )
    headers = {"accept": "application/json", "content-type": "application/x-netcdf"}
    with open(dataset_path, "rb") as handle:
        response = requests.put(upload_url, headers=headers, data=handle, auth=auth, timeout=300)

    return response


def _create_netcdf_store_from_upload(geoserver_url, workspace, dataset_path, base_store_name, auth):
    last_error = ""
    for suffix in range(0, 20):
        store_name = base_store_name if suffix == 0 else f"{base_store_name}_{suffix}"
        response = _upload_netcdf_store(geoserver_url, workspace, dataset_path, store_name, auth)

        if response.status_code in SUCCESS_CODES:
            print(f"Store NetCDF publicado por subida directa: {store_name}")
            return store_name

        last_error = f"{response.status_code}: {response.text}"
        response_text = (response.text or "").lower()
        if response.status_code == 409 or "already exists" in response_text:
            continue

        raise RuntimeError(
            f"Error al subir NetCDF al store '{store_name}' ({response.status_code}): {response.text}"
        )

    raise RuntimeError(
        "No se pudo crear un nombre único para subir el NetCDF por API REST. "
        f"Último error: {last_error}"
    )


def create_coverage_store(yaml_file, output_file):
    with open(yaml_file, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    dataset_cfg = config.get("dataset", {}) or {}
    geoserver_cfg = (config.get("services", {}) or {}).get("geoserver", {}) or {}
    metadata_cfg = config.get("metadata", {}) or {}

    resource_format = str(dataset_cfg.get("resourceFormat", "")).strip()
    dataset_path = dataset_cfg.get("file")
    publish_enabled = geoserver_cfg.get("publish_geoserver")
    if publish_enabled is None:
        publish_enabled = geoserver_cfg.get("enabled", False)

    if not publish_enabled:
        print("Publicación en GeoServer omitida (publish_geoserver=false).")
        Path(output_file).write_text("", encoding="utf-8")
        return

    geoserver_url = str(geoserver_cfg.get("url", "")).rstrip("/")
    workspace = geoserver_cfg.get("workspace")
    username = geoserver_cfg.get("username")
    password = geoserver_cfg.get("password")
    style = geoserver_cfg.get("style")
    netcdf_store_url = geoserver_cfg.get("netcdf_store_url")

    if not geoserver_url or not workspace or not username or not password:
        raise ValueError("Faltan credenciales/URL/workspace de GeoServer en services.geoserver.")
    if not dataset_path:
        raise ValueError("Falta dataset.file en la configuración.")
    if not Path(dataset_path).exists():
        raise FileNotFoundError(f"No existe el dataset para publicar en GeoServer: {dataset_path}")

    title = metadata_cfg.get("title") or Path(dataset_path).stem
    base_store = slugify_name(title, fallback=slugify_name(Path(dataset_path).stem))
    auth = (username, password)

    if resource_format == "GeoTIFF":
        layer_name = _create_geotiff_store(geoserver_url, workspace, dataset_path, base_store, auth)
    elif resource_format == "NetCDF":
        # Preferimos subir el fichero por REST para no depender de rutas file:// visibles desde GeoServer.
        netcdf_upload_error = None
        try:
            store_name = _create_netcdf_store_from_upload(
                geoserver_url, workspace, dataset_path, base_store, auth
            )
            coverage_names = _list_store_coverages(geoserver_url, workspace, store_name, auth)
        except RuntimeError as exc:
            netcdf_upload_error = exc
            coverage_names = []

            if netcdf_store_url:
                print(
                    "Aviso: falló la subida directa del NetCDF. "
                    "Se intenta publicar usando services.geoserver.netcdf_store_url."
                )
                store_name = _create_netcdf_store(
                    geoserver_url, workspace, str(netcdf_store_url), base_store, auth
                )
                coverage_names = _list_store_coverages(geoserver_url, workspace, store_name, auth)

        if not coverage_names:
            print(
                "Aviso: GeoServer no pudo publicar el NetCDF como cobertura nativa. "
                f"Se publicará una capa WMS GeoTIFF derivada. Error original: {netcdf_upload_error}"
            )
            derived_tif = Path(output_file).with_name(f"{base_store}_geoserver.tif")
            geotiff_path = _materialize_netcdf_as_geotiff(dataset_path, config, derived_tif)
            layer_name = _create_geotiff_store(
                geoserver_url, workspace, geotiff_path, base_store, auth
            )
        else:
            layer_name = coverage_names[0]
            print(f"Coberturas publicadas en store '{store_name}': {coverage_names}")
    else:
        raise ValueError(f"Formato no soportado para publicación GeoServer: {resource_format}")

    _style_layer(geoserver_url, workspace, layer_name, style, auth)
    _write_wms_output(output_file, geoserver_url, workspace, layer_name, config, dataset_path)


if __name__ == "__main__":
    try:
        yaml_file = snakemake.input.config
        output_file = snakemake.output[0]
    except NameError:
        if len(sys.argv) < 3:
            print("Uso: publish_geoserver.py <config_yaml> <output_file>")
            sys.exit(1)
        yaml_file = sys.argv[1]
        output_file = sys.argv[2]

    create_coverage_store(yaml_file, output_file)
