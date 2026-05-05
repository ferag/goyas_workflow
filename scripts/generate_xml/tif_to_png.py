#!/usr/bin/env python3
import sys

import matplotlib.pyplot as plt
import numpy as np
import rasterio


def _normalize_band(band):
    band = band.astype("float32")
    mask = np.isfinite(band)
    if not mask.any():
        return np.zeros_like(band, dtype="float32")
    values = band[mask]
    lo, hi = np.percentile(values, [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(values.min()), float(values.max())
    if hi <= lo:
        return np.zeros_like(band, dtype="float32")
    scaled = (band - lo) / (hi - lo)
    scaled = np.clip(scaled, 0.0, 1.0)
    scaled[~mask] = 0.0
    return scaled


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


def generate_png(dataset_path, output_png, rgb_bands=None):
    """
    Genera un PNG de vista previa a partir de un raster (GeoTIFF o NetCDF).
    Si el NetCDF no tiene bandas en el dataset raíz, usa el primer subdataset válido.
    """
    source_path = _select_raster_source(dataset_path)
    with rasterio.open(source_path) as src:
        num_bands = src.count
        width = max(src.width, 1)
        height = max(src.height, 1)
        aspect_ratio = height / width
        fig_width = 8
        fig_height = max(fig_width * aspect_ratio, 1.0)

        plt.figure(figsize=(fig_width, fig_height))

        if num_bands == 1:
            band = _normalize_band(src.read(1))
            plt.imshow(band, cmap="gray", interpolation="none", aspect="equal")
        elif num_bands >= 3:
            if rgb_bands is None:
                rgb_bands = [1, 2, 3]
            valid_bands = [idx for idx in rgb_bands if 1 <= idx <= num_bands]
            if len(valid_bands) < 3:
                valid_bands = [1, min(2, num_bands), min(3, num_bands)]
            rgb = src.read(valid_bands[:3]).astype("float32")
            rgb = np.stack([_normalize_band(rgb[i]) for i in range(3)], axis=0)
            rgb = rgb.transpose((1, 2, 0))
            plt.imshow(rgb, interpolation="none", aspect="equal")
        else:
            band = _normalize_band(src.read(1))
            plt.imshow(band, cmap="gray", interpolation="none", aspect="equal")

        plt.axis("off")
        plt.savefig(str(output_png), bbox_inches="tight", pad_inches=0)
        plt.close()
        print(f"PNG snapshot generado desde {source_path}: {output_png}")


if __name__ == "__main__":
    try:
        input_dataset = snakemake.input.data
        output_png = snakemake.output[0]
    except NameError:
        if len(sys.argv) < 3:
            print("Uso: tif_to_png.py <input_raster> <output_png>")
            sys.exit(1)
        input_dataset = sys.argv[1]
        output_png = sys.argv[2]
    generate_png(input_dataset, output_png)
