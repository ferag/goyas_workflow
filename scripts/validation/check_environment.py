import importlib
import json
import sys
from pathlib import Path
from datetime import datetime, UTC

min_python_version = (3, 9)
root_dir = Path(__file__).resolve().parents[2]
requirements_file = root_dir / "requirements.txt" # Check if requirements.txt in the parent directory

# Directorios requeridos
required_dirs = [root_dir / "logs",
                 root_dir / "tmp",
                 root_dir / "output"]

def read_requirements(file_path):
    if not file_path.exists():
        raise FileNotFoundError(f"Requirements file not found: {file_path}")
    
    packages = []
    for line in file_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        pkg = line.split("==")[0].split(">=")[0].split("<=")[0].strip()
        if pkg:
            packages.append(pkg)
    return packages

def check_package(packages):
    results = {}
    PACKAGE_ALIASES = {
        "pyyaml": "yaml",
        "python-dotenv": "dotenv",
    }
    for package in packages:
        try:
            importlib.import_module(PACKAGE_ALIASES.get(package, package))
            results[package] = True
        except ModuleNotFoundError:
            results[package] = "NOT FOUND"
        except Exception as e:
            results[package] = f"ERROR: {str(e)}"
    return results

def check_python_version(min_version):
    current_version = sys.version_info[:3]
    if current_version < min_version:
        return False, f"Python {min_version[0]}.{min_version[1]} or higher is required, but you have {current_version[0]}.{current_version[1]}.{current_version[2]}"
    return True, f"Python version {current_version[0]}.{current_version[1]}.{current_version[2]} is sufficient."

def setup_directories(directories):
    created, existing, errors = [], [], []
    for dir_name in directories:
        dir_path = Path(dir_name)
        try:
            if not dir_path.exists():
                dir_path.mkdir(parents=True, exist_ok=True)
                created.append(str(dir_path.relative_to(root_dir)))
            else:
                existing.append(str(dir_path.relative_to(root_dir)))
        except Exception as e:
            errors.append(f"Failed to create directory '{dir_name}': {str(e)}")
    return {"created": created, "existing": existing, "errors": errors}

def check_environment():
    report = {
        "timestamp": datetime.now(UTC).isoformat(),
        "python_version": {},
        "packages": {},
        "directories": {},
        "ok": True,
        "errors": [],
        "warnings": []
    }

    # Check Python version
    py_check, py_message = check_python_version(min_python_version)
    report["python_version"]["check"] = py_check
    report["python_version"]["message"] = py_message
    if not py_check:
        report["ok"] = False
        report["errors"].append(py_message)

    # Check required packages
    try:
        packages = read_requirements(requirements_file)
    except Exception as e:
        report["ok"] = False
        report["errors"].append(f"Failed to read requirements: {str(e)}")
        return report
    
    if not packages:
        report["warnings"].append("No packages found in requirements file.")
        return report
    
    package_results = check_package(packages)
    report["packages"] = package_results
    for pkg, status in package_results.items():
        if status is not True:
            report["ok"] = False
            report["errors"].append(f"Package '{pkg}': {status}")

    # Setup required directories
    dir_status = setup_directories(required_dirs)
    report["directories"] = dir_status
    if dir_status["errors"]:
        report["ok"] = False
        for err in dir_status["errors"]:
            report["errors"].append(err)

    return report

if __name__ == "__main__":
    env_report = check_environment()

    logs_dir = root_dir / "logs"
    output_path = Path(logs_dir / "env_check.json")
    output_path.write_text(json.dumps(env_report, indent=2, ensure_ascii=False))
    print(json.dumps(env_report, indent=2, ensure_ascii=False))

    if not env_report["ok"]:
        sys.exit(1)