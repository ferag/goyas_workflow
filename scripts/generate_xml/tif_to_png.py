#!/usr/bin/env python3
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib.pyplot as plt
import numpy as np
import rasterio
import yaml


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


def _thumbnail_config(config):
    metadata_cfg = config.get("metadata", {}) if isinstance(config, dict) else {}
    thumb_cfg = metadata_cfg.get("thumbnail", {}) if isinstance(metadata_cfg, dict) else {}
    if not isinstance(thumb_cfg, dict):
        thumb_cfg = {}
    colormap = str(thumb_cfg.get("colormap", "turbo")).strip() or "turbo"
    if colormap not in plt.colormaps():
        raise ValueError(f"Colormap no soportado por matplotlib: {colormap}")

    return {
        "mode": str(thumb_cfg.get("mode", "scalar")).strip().lower(),
        "colormap": colormap,
        "percentile_min": float(thumb_cfg.get("percentile_min", 2)),
        "percentile_max": float(thumb_cfg.get("percentile_max", 98)),
        "band": int(thumb_cfg.get("band", 1)),
        "rgb_bands": thumb_cfg.get("rgb_bands") or [1, 2, 3],
        "crop_to_data": thumb_cfg.get("crop_to_data", True) is not False,
        "dpi": int(thumb_cfg.get("dpi", 120)),
    }


def _normalize_band(band, percentile_min=2, percentile_max=98):
    band = band.astype("float32")
    mask = np.isfinite(band)
    if not mask.any():
        return np.zeros_like(band, dtype="float32")
    values = band[mask]
    lo, hi = np.percentile(values, [percentile_min, percentile_max])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(values.min()), float(values.max())
    if hi <= lo:
        return np.zeros_like(band, dtype="float32")
    scaled = (band - lo) / (hi - lo)
    scaled = np.clip(scaled, 0.0, 1.0)
    scaled[~mask] = np.nan
    return scaled


def _crop_to_valid(data):
    if data.ndim == 2:
        valid = np.isfinite(data)
    else:
        valid = np.any(np.isfinite(data), axis=0)

    rows, cols = np.where(valid)
    if rows.size == 0 or cols.size == 0:
        return data

    row_min, row_max = rows.min(), rows.max() + 1
    col_min, col_max = cols.min(), cols.max() + 1
    if data.ndim == 2:
        return data[row_min:row_max, col_min:col_max]
    return data[:, row_min:row_max, col_min:col_max]


def generate_png(dataset_path, output_png, config=None):
    """
    Genera un PNG de vista previa a partir de un raster (GeoTIFF o NetCDF).
    Si el NetCDF no tiene bandas en el dataset raíz, usa el primer subdataset válido.
    """
    thumb_cfg = _thumbnail_config(config or {})
    source_path = _select_raster_source(dataset_path)
    with rasterio.open(source_path) as src:
        num_bands = src.count

        if num_bands == 1 or thumb_cfg["mode"] == "scalar":
            band_idx = min(max(thumb_cfg["band"], 1), num_bands)
            raw_band = src.read(band_idx)
            if thumb_cfg["crop_to_data"]:
                raw_band = _crop_to_valid(raw_band)
            band = _normalize_band(
                raw_band,
                thumb_cfg["percentile_min"],
                thumb_cfg["percentile_max"],
            )
            cmap = plt.get_cmap(thumb_cfg["colormap"]).copy()
            cmap.set_bad(alpha=0)
            image_data = band
            image_kwargs = {"cmap": cmap, "vmin": 0, "vmax": 1}
        elif num_bands >= 3:
            rgb_bands = thumb_cfg["rgb_bands"]
            valid_bands = [idx for idx in rgb_bands if 1 <= idx <= num_bands]
            if len(valid_bands) < 3:
                valid_bands = [1, min(2, num_bands), min(3, num_bands)]
            rgb = src.read(valid_bands[:3]).astype("float32")
            if thumb_cfg["crop_to_data"]:
                rgb = _crop_to_valid(rgb)
            rgb = np.stack(
                [
                    _normalize_band(
                        rgb[i],
                        thumb_cfg["percentile_min"],
                        thumb_cfg["percentile_max"],
                    )
                    for i in range(3)
                ],
                axis=0,
            )
            rgb = np.nan_to_num(rgb, nan=0.0)
            image_data = rgb.transpose((1, 2, 0))
            image_kwargs = {}
        else:
            raw_band = src.read(1)
            if thumb_cfg["crop_to_data"]:
                raw_band = _crop_to_valid(raw_band)
            image_data = _normalize_band(
                raw_band,
                thumb_cfg["percentile_min"],
                thumb_cfg["percentile_max"],
            )
            image_kwargs = {"cmap": thumb_cfg["colormap"]}

        height, width = image_data.shape[:2]
        aspect_ratio = max(height, 1) / max(width, 1)
        fig_width = 8
        fig_height = max(fig_width * aspect_ratio, 1.0)

        plt.figure(figsize=(fig_width, fig_height))
        plt.imshow(image_data, interpolation="none", aspect="equal", **image_kwargs)

        plt.axis("off")
        plt.savefig(str(output_png), bbox_inches="tight", pad_inches=0, dpi=thumb_cfg["dpi"])
        plt.close()
        print(
            "PNG snapshot generado "
            f"desde {source_path} con colormap={thumb_cfg['colormap']}: {output_png}"
        )


if __name__ == "__main__":
    try:
        input_dataset = snakemake.input.data
        config_file = snakemake.input.config
        output_png = snakemake.output[0]
    except NameError:
        if len(sys.argv) < 3:
            print("Uso: tif_to_png.py <input_raster> <output_png> [config.yaml]")
            sys.exit(1)
        input_dataset = sys.argv[1]
        output_png = sys.argv[2]
        config_file = sys.argv[3] if len(sys.argv) > 3 else None

    config = {}
    if config_file:
        with open(config_file, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
    generate_png(input_dataset, output_png, config)
