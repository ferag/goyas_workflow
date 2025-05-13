#!/usr/bin/env python3
import rasterio
import matplotlib.pyplot as plt
import sys

import rasterio
import matplotlib.pyplot as plt

def generate_png(input_tif, output_png, rgb_bands=None):
    """
    Genera un snapshot PNG a partir de un GeoTIFF:
      - Si tiene 1 banda, se muestra en escala de grises.
      - Si tiene 3 bandas, se muestra en composición RGB.
      - Si tiene más de 3 bandas, se utiliza la composición definida en rgb_bands (lista de 3 índices).
    """
    with rasterio.open(input_tif) as src:
        num_bands = src.count
        width = src.width
        height = src.height
        aspect_ratio = height / width
        fig_width = 8  # ancho en pulgadas (valor arbitrario)
        fig_height = fig_width * aspect_ratio

        plt.figure(figsize=(fig_width, fig_height))

        if num_bands == 1:
            # Imagen en escala de grises
            band1 = src.read(1)
            plt.imshow(band1, cmap='gray', interpolation='none', aspect='equal')
        elif num_bands == 3:
            # Composición RGB utilizando las tres bandas
            rgb = src.read([1, 2, 3])
            # Transponer de (bands, height, width) a (height, width, bands)
            rgb = rgb.transpose((1, 2, 0))
            plt.imshow(rgb, interpolation='none', aspect='equal')
        elif num_bands > 3:
            # Si se especifican índices para la composición RGB, se usan; sino, se usan las tres primeras bandas
            if rgb_bands is None:
                rgb_bands = [1, 2, 3]
            rgb = src.read(rgb_bands)
            rgb = rgb.transpose((1, 2, 0))
            plt.imshow(rgb, interpolation='none', aspect='equal')
        else:
            raise ValueError("El GeoTIFF debe tener al menos 1 banda.")

        plt.axis('off')  # Oculta ejes y leyenda
        plt.savefig(str(output_png), bbox_inches='tight', pad_inches=0)
        plt.close()
        print(f"PNG snapshot generado: {output_png}")

if __name__ == "__main__":
    try:
        # Si se ejecuta desde Snakemake, se accede a los inputs/outputs mediante snakemake.
        input_tif = snakemake.input.data
        output_png = snakemake.output
    except NameError:
        # Modo de ejecución manual: se esperan argumentos de línea de comandos
        if len(sys.argv) < 3:
            print("Uso: tif_to_png.py <input_tif> <output_png>")
            sys.exit(1)
        input_tif = sys.argv[1]
        output_png = sys.argv[2]
    generate_png(input_tif, output_png)
