import re
import sys
from pathlib import Path
from datetime import datetime, date, timezone
import yaml
import json
import traceback

try:
    from zoneinfo import ZoneInfo
    zone_info = True
except Exception:
    zone_info = False

lang_enum = {"spa", "eng", "cat", "fra"}
status_enum = {"completed", "ongoing", "planned", "deprecated"}
freq_enum = {"notPlanned", "daily", "weekly", "monthly", "annually", "asNeeded"}
role_enum = {"author", "coAuthor", "collaborator", "contributor", "custodian", "distributor", "editor", "funder", "mediator", "originator",
              "owner", "pointOfContact", "principalInvestigator", "processor", "publisher", "resourceProvider", "rightsHolder", "sponsor", "stakeholder", "user"}
resource_enum = {"GeoTIFF", "NetCDF"}

re_email = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
re_url = re.compile(r"^https?://", re.IGNORECASE)

root_dir = Path(__file__).resolve().parents[2]

def _as_date(s):
    return date.fromisoformat(s)

def _is_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)

def validate_config_yaml(cfg_path):
    """Validate the configuration YAML file."""

    warnings, errors = [], []
    path = Path(cfg_path)
    if not path.exists():
        return {"ok": False, "errors": [f"File not found: {cfg_path}"], "warnings": []}
    
    try:
        cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    except Exception as e:
        return {"ok": False, "errors": [f"YAML parsing error: {e}"], "warnings": []}
    
    def g(*keys, default=None):
        d = cfg
        for key in keys:
            if not isinstance(d, dict) or key not in d:
                return default
            d = d[key]
        return d
    
    # Dataset
    dataset_file = g("dataset", "file")
    resource_format = g("dataset", "resourceFormat")
    if not dataset_file or not isinstance(dataset_file, str):
        errors.append("dataset.file is required and must be a string.")
    elif not Path(dataset_file).exists():
        errors.append(f"dataset.file does not exist: {dataset_file}")
    if resource_format not in resource_enum:
        errors.append(f"dataset.resourceFormat must be one of {resource_enum}.")

    if dataset_file and isinstance(dataset_file, str):
        suf = Path(dataset_file).suffix.lower()
        if resource_format == "GeoTIFF" and suf not in {".tif", ".tiff"}:
            warnings.append(f"dataset.file extension does not match GeoTIFF format: {suf}")
        if resource_format == "NetCDF" and suf not in {".nc", ".cdf"}:
            warnings.append(f"dataset.file extension does not match NetCDF format: {suf}")

    # Metadata
    for k in ("title", "abstract", "identifier", "language"):
        v = g("metadata", k)
        if not v or not isinstance(v, str):
            errors.append(f"metadata.{k} is required and must be a string.")

    language = g("metadata", "language")
    if language not in lang_enum:
        errors.append(f"metadata.language must be one of {lang_enum}.")

    # Contact
    contact = g("metadata", "contact")
    if not contact:
        errors.append("metadata.contact is required.")
    else:
        for k in ("individualName", "organisationName", "role", "email"):
            if not contact.get(k):
                errors.append(f"metadata.contact.{k} is required.")
        email = contact.get("email")
        if email and not re_email.match(email):
            errors.append(f"metadata.contact.email is not a valid email address: {email}")
        role = contact.get("role")
        if role and role not in role_enum:
            errors.append(f"metadata.contact.role must be one of {role_enum}.")

    # temporal Coverage
    begin = g("metadata", "coverage", "temporal", "begin")
    end = g("metadata", "coverage", "temporal", "end")
    timezone = g("metadata", "coverage", "temporal", "timezone", "UTC")

    d1 = _as_date(begin) if begin else None
    d2 = _as_date(end) if end else None
    if (not d1 or not d2) or (d1 > d2):
        errors.append("metadata.coverage.temporal.begin and end must be valid dates with begin <= end.")

    if timezone:
        if not isinstance(timezone, str):
            errors.append("metadata.coverage.temporal.timezone must be a string.")
        elif zone_info:
            try:
                ZoneInfo(timezone)
            except Exception:
                warnings.append(f"metadata.coverage.temporal.timezone is not a valid timezone: {timezone}")

    # spatial Coverage
    spatial = g("metadata", "coverage", "spatial")
    if not spatial:
        warnings.append("metadata.coverage.spatial is not defined.")
    elif not spatial.get("auto", False):
        for k, (mn, mx) in {"west": (-180, 180), "east": (-180, 180), "south": (-90, 90), "north": (-90, 90)}.items():
            v = spatial.get(k)
            if v is None:
                errors.append(f"metadata.coverage.spatial.{k} is required if auto=False.")
            elif not _is_number(v) or not (mn <= float(v) <= mx):
                errors.append(f"metadata.coverage.spatial.{k} out of range. Must be a number between {mn} and {mx}.")

    # crs
    epsg = g("crs", "epsg")
    if not isinstance(epsg, int) or epsg <= 0:
        errors.append("metadata.epsg must be a positive integer (Ej., 4326, 32630).")

    # Licence
    licence = g("metadata", "licence")
    if not licence or not isinstance(licence, dict):
        errors.append("metadata.licence is required.")
    else:
        if not licence.get("id"):
            errors.append("metadata.licence.id is required (SPDX).")
        if not (isinstance(licence.get("url"), str) and re_url.match(licence.get("url"))):
            errors.append("metadata.licence.url is required and must be a valid URL.")

    # status
    status = g("metadata", "status")
    if status and status not in status_enum:
        errors.append(f"metadata.status must be one of {status_enum}.")

    # frequency
    frequency = g("metadata", "maintenanceFrequency")
    if frequency and frequency not in freq_enum:
        errors.append(f"metadata.maintenanceFrequency must be one of {freq_enum}.")

    # resolution
    res = g("metadata", "resolution")
    if isinstance(res, dict):
         for k in ("x", "y"):
             v = res.get(k)
             if v is not None and (not _is_number(v) or float(v) <= 0):
                 errors.append(f"metadata.resolution.{k} must be a positive number if defined.")

    # Services
    geonetwork = g("services", "geonetwork")
    if isinstance(geonetwork, dict):
        url = geonetwork.get("url")
        if url and not re_url.match(url):
            errors.append(f"services.geonetwork.url is not a valid URL: {url}")

    geoserver = g("services", "geoserver")
    if isinstance(geoserver, dict):
        url = geoserver.get("url")
        if url and not re_url.match(url):
            errors.append(f"services.geoserver.url is not a valid URL: {url}")
        pg = geoserver.get("publish_geoserver")
        if pg is not None and not isinstance(pg, bool):
            errors.append("services.geoserver.publish_geoserver must be a boolean if defined.")

    # Storage
    storage = g("storage")
    if isinstance(storage, dict):
        eu = storage.get("endpoint_url")
        if eu and not re_url.match(eu):
            errors.append(f"storage.endpoint_url is not a valid URL: {eu}")
        if not storage.get("bucket"):
            warnings.append("storage.bucket is not defined.")

    return {"ok": len(errors) == 0, "errors": errors, "warnings": warnings}

if __name__ == "__main__":
    # When executed via Snakemake `script:`, a `snakemake` object is injected in globals().
    if "snakemake" in globals():
        cfg_path = Path(snakemake.input[0]).resolve()
        log_path = Path(snakemake.output[0]).resolve()
    else:
        cfg_arg = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
        cfg_path = (root_dir / cfg_arg).resolve()
        log_path = (root_dir / "logs" / "config_validation.json").resolve()

    try:
        result = validate_config_yaml(cfg_path)
    except Exception as e:
        # Unexpected crash: log it and fail hard so the workflow stops.
        result = {"ok": False, "errors": [f"Unhandled exception: {e}"], "warnings": []}
        print("ERROR:", e, file=sys.stderr)
        traceback.print_exc()

    log = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "config_file": str(cfg_path),
        "ok": result["ok"],
        "errors": result["errors"],
        "warnings": result["warnings"],
    }

    # Atomic log write (avoid leaving a partial file on crash)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = log_path.with_suffix(log_path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp_path.replace(log_path)

    if result["ok"]:
        print("Configuration is valid.")
    else:
        print("Configuration has errors:", file=sys.stderr)
        for err in result["errors"]:
            print(f"  - {err}", file=sys.stderr)

    if result["warnings"]:
        print("Warnings:")
        for warn in result["warnings"]:
            print(f"  - {warn}")

    # Critical: make Snakemake stop on validation errors
    if not result["ok"]:
        sys.exit(1)
