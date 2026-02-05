#!/usr/bin/env python3
import requests
import sys
import yaml
import json
import rasterio
from rasterio.warp import transform_bounds
import unicodedata

def eliminar_acentos(texto):
    # Descompone los caracteres (la 'é' se convierte en 'e' + '´')
    forma_nfkd = unicodedata.normalize('NFKD', texto)
    # Filtra solo los caracteres que no sean marcas de acentuación
    solo_ascii = "".join([c for c in forma_nfkd if not unicodedata.combining(c)])
    return solo_ascii

def create_coverage_store(yaml_file, output_file):
    # Parámetros de conexión y configuración
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)

    geoserver_url = config.get("services")['geoserver']['url']
    workspace = config.get("services")['geoserver']['workspace']
    coveragestore = eliminar_acentos(config.get("metadata").get("title"))[0:10]
    style = config.get("services")['geoserver']["style"]
    username = config.get("services")['geoserver']['username']
    password = config.get("services")['geoserver']['password']
    tif_path = config.get("dataset")['file']

    status_code = 0
    while(status_code not in [200, 201]):
      # Payload XML, utilizando la variable coveragestore en el <name> y en la URL de coverages
      xml_payload = f"""<?xml version="1.0" encoding="UTF-8"?>
    <coverageStore>
        <name>{coveragestore}</name>
        <description>First test API</description>
        <type>GeoTIFF</type>
        <enabled>true</enabled>
        <workspace>
            <name>{workspace}</name>
            <link>https://goyas.csic.es/geoserver/rest/workspaces/{workspace}</link>
        </workspace>
        <__default__>true</__default__>
        <url>file:data/{tif_path}</url>
        <coverages>
            <link>https://goyas.csic.es/geoserver/rest/workspaces/{workspace}/coveragestores/{coveragestore}/coverages.json</link>
        </coverages>
    </coverageStore>"""

      # URL para crear el coverage store
      url = f"{geoserver_url}/rest/workspaces/{workspace}/coveragestores"
      headers = {
          "accept": "application/json",
          "content-type": "application/xml"
      }
      auth = (username, password)

      response = requests.post(url, headers=headers, data=xml_payload, auth=auth)
      print("Código de respuesta:", response.status_code)
      if response.status_code in [200, 201]:
          print("Coverage store creado correctamente.")
          status_code = response.status_code
          print(response.text)
      else:
          print("Error al crear coverage store:")
          print(response.text)
          coveragestore = coveragestore + "_" + coveragestore[0]
          status_code = response.status_code
          if len(coveragestore) > 111:
            sys.exit(1)

    # 2 subo el fichero
    upload_url = f"{geoserver_url}/rest/workspaces/{workspace}/coveragestores/{coveragestore}/file.geotiff?filename={tif_path}&updateBBox=true"

    # Encabezado para indicar que se envía un TIFF
    headers = {
        "Content-Type": "image/tiff",
        "accept": "application/json"
    }

    with open(tif_path, "rb") as f:
        tif_data = f.read()

    response = requests.put(upload_url, headers=headers, data=tif_data, auth=auth)

    print("Código de respuesta:", response.status_code)
    if response.status_code in [200, 201]:
        print("GeoTIFF subido y coverage store creado correctamente.")
        print(response.text)
    else:
        print("Error al subir el GeoTIFF:")
        print(response.text)
    print(response.text)

    coverage_name = coveragestore
    coverage_title = tif_path

    with rasterio.open(tif_path) as src:
        width = src.width
        height = src.height
        bounds = src.bounds  # left, bottom, right, top
        crs = src.crs.to_string() if src.crs else "EPSG:4326"
        print(f"CRS: {crs}")
        transform = src.transform

        # Extraer parámetros de la transformación affine
        scaleX = transform.a
        shearX = transform.b
        translateX = transform.c
        shearY = transform.d
        scaleY = transform.e
        translateY = transform.f

        # Calcular el bounding box en EPSG:4326 (lat-lon)
        if src.crs:
            ll_bounds = transform_bounds(src.crs, "EPSG:4326", *bounds, densify_pts=21)
            crs = "EPSG:4326"
        else:
            ll_bounds = bounds

    print(f"Bounds calculados: {bounds} - {ll_bounds}")

    print("Número de bandas:", src.count)

    if src.count == 3:
        dimensions = {
          "coverageDimension":[
            {
              "name":"RED_BAND",
              "description":"GridSampleDimension[0.0,0.0]",
              "range":{
                "min":0,
                "max":0
              },
              "nullValues":{
                "double":0
              }
            },
            {
              "name":"GREEN_BAND",
              "description":"GridSampleDimension[0.0,0.0]",
              "range":{
                "min":0,
                "max":0
              },
              "nullValues":{
                "double":0
              }
            },
            {
              "name":"BLUE_BAND",
              "description":"GridSampleDimension[0.0,0.0]",
              "range":{
                "min":0,
                "max":0
              },
              "nullValues":{
                "double":0
              }
            }
          ]
        }
    else:
       dimensions = {
          "coverageDimension":[
            {
              "name":"GREY_BAND",
              "description":"GridSampleDimension[0.0,0.0]",
              "range":{
                "min":0,
                "max":0
              },
              "nullValues":{
                "double":0
              }
            }
          ]
        }

    payload = {
    "coverage": {
        "name": coverage_name,
        "title": coverage_title,
        "abstract": "Abstract generated from GeoTIFF metadata",
        "defaultInterpolationMethod": "nearest neighbor",
        "description": "Coverage generated from GeoTIFF file",
        "enabled": True,
        "nativeFormat": "GeoTIFF",
        "nativeName": coverage_name,
        "nativeCRS": {
            "$": crs,
            "@class": "projected"
        },
        "grid": {
        "@dimension": "2",
        "crs": crs,
        "range": {
            "low": "0 0",
            "high": f"{width} {height}"
        },
        "transform": {
            "scaleX": scaleX,
            "shearX": shearX,
            "translateX": translateX,
            "shearY": shearY,
            "scaleY": scaleY,
            "translateY": translateY
        }
        },
        "latLonBoundingBox": {
        "crs": "EPSG:4326",
        "minx": ll_bounds[0],
        "miny": ll_bounds[1],
        "maxx": ll_bounds[2],
        "maxy": ll_bounds[3]
        },
        "nativeBoundingBox": {
        "crs": crs,
        "minx": bounds.left,
        "miny": bounds.bottom,
        "maxx": bounds.right,
        "maxy": bounds.top
        },
    "dimensions": dimensions,
    "supportedFormats":{
      "string":[
        "GIF",
        "PNG",
        "JPEG",
        "TIFF",
        "ImageMosaic",
        "GEOTIFF",
        "ArcGrid"
      ]
    },
    "requestSRS":{
      "string":"EPSG:4326"
    },
    "responseSRS":{
      "string":"EPSG:4326"
    },
    "parameters":{
      "entry":[
        {
          "string":"BackgroundValues",
          "null":""
        },
        {
          "string":"Filter",
          "null":""
        },
        {
          "string":[
            "USE_JAI_IMAGEREAD",
            "true"
          ]
        }
      ]
    }
    }
    }

    url = f'https://goyas.csic.es/geoserver/rest/workspaces/goyas/coveragestores/{coveragestore}/coverages'

    auth = (username, password)
    headers = {
        "accept": "application/json",
        "content-type": "application/json"
    }
    response = requests.post(url, headers=headers, data=json.dumps(payload), auth=auth)

    print("Código de respuesta:", response.status_code)
    if response.status_code in [200, 201]:
        print("Coverage store creado correctamente.")
        print(response.text)
    else:
        print("Error al crear coverage store:")
        print(response.text)

    url = f'https://goyas.csic.es/geoserver/rest/layers/{workspace}:{coveragestore}'

    payload = {
        "layer": {
          "defaultStyle": {
            "name": f"{style}"
          }
        }
    }
    headers = {
        "accept": "application/json",
        "content-type": "application/json"
    }
    requests.put(url, headers=headers, data=json.dumps(payload), auth=auth)


    with open(str(output_file), 'w', encoding='utf-8') as f:
        ogc_url = f"https://goyas.csic.es/geoserver/{workspace}/ows?SERVICE=WMS&SERVICE=WMS&VERSION=1.3.0&REQUEST=GetCapabilities&FORMAT=image%2Fpng&TRANSPARENT=true&LAYERS={coveragestore}1&STYLES=&CRS=EPSG%3A3857&WIDTH=2040&HEIGHT=892&BBOX={bounds.left}%2C{bounds.bottom}%2C{bounds.right}%2C{bounds.top}&width=768&height=448&srs=EPSG%3A4326"
        f.write(ogc_url)

if __name__ == "__main__":
    try:
        yaml_file = snakemake.input.config
        output_file = snakemake.output
    except NameError:
        # Modo de ejecución manual
        if len(sys.argv) < 4:
            print("Uso: upload_metadata_initial.py <session_file> <xml_file_path> <output_file>")
            sys.exit(1)
            session_file = sys.argv[1]
            xml_file_path = sys.argv[2]
            output_file = sys.argv[3]
    create_coverage_store(yaml_file, output_file)
