#!/usr/bin/env python3
import argparse
import re
import subprocess
import sys
from pathlib import Path


PLACEHOLDER_PATTERNS = [
    re.compile(r"^\*+$"),
    re.compile(r"^x+$", re.IGNORECASE),
    re.compile(r"^<.*>$"),
    re.compile(r"^(replace|change|your|example|dummy|test)", re.IGNORECASE),
]

KEY_VALUE_RE = re.compile(
    r"""(?ix)
    (?<![A-Za-z0-9_])
    (?P<key>
      password|passwd|secret|api[_-]?key|access[_-]?key|client[_-]?secret|
      private[_-]?key|access[_-]?token|refresh[_-]?token|bearer[_-]?token
    )
    (?![A-Za-z0-9_])
    \s*[:=]\s*
    (?P<quote>["'])
    (?P<value>[^"']+)
    (?P=quote)
    """
)

URL_CREDENTIAL_RE = re.compile(r"(?i)https?://[^/\s:@]+:[^/\s@]+@")
AWS_ACCESS_KEY_RE = re.compile(r"\bAKIA[0-9A-Z]{16}\b")
PRIVATE_KEY_BLOCK_RE = re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----")


def run_git(args):
    return subprocess.run(
        ["git", *args],
        check=False,
        text=True,
        capture_output=True,
    )


def list_staged_files():
    proc = run_git(["diff", "--cached", "--name-only", "--diff-filter=ACMRTUXB"])
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "No se pudieron listar archivos staged.")
    return [line.strip() for line in proc.stdout.splitlines() if line.strip()]


def list_tracked_files():
    proc = run_git(["ls-files"])
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "No se pudieron listar archivos trackeados.")
    return [line.strip() for line in proc.stdout.splitlines() if line.strip()]


def looks_like_placeholder(value):
    normalized = value.strip()
    if not normalized:
        return True
    for pattern in PLACEHOLDER_PATTERNS:
        if pattern.search(normalized):
            return True
    return False


def is_probably_text(path: Path):
    try:
        chunk = path.read_bytes()[:4096]
    except OSError:
        return False
    return b"\x00" not in chunk


def scan_file(path: Path):
    issues = []
    if not path.exists() or not path.is_file():
        return issues
    if not is_probably_text(path):
        return issues

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        try:
            lines = path.read_text(encoding="latin-1").splitlines()
        except Exception:
            return issues

    for idx, line in enumerate(lines, start=1):
        lowered = line.lower()
        if "secret-scan: allow" in lowered:
            continue

        if PRIVATE_KEY_BLOCK_RE.search(line):
            issues.append((idx, "Bloque de clave privada detectado."))
            continue

        if URL_CREDENTIAL_RE.search(line):
            issues.append((idx, "URL con credenciales embebidas detectada."))

        if AWS_ACCESS_KEY_RE.search(line):
            issues.append((idx, "AWS Access Key detectada."))

        for match in KEY_VALUE_RE.finditer(line):
            value = match.group("value")
            if looks_like_placeholder(value):
                continue
            issues.append(
                (
                    idx,
                    f"Posible secreto en clave '{match.group('key')}' (valor no enmascarado).",
                )
            )

    return issues


def main():
    parser = argparse.ArgumentParser(
        description="Bloquea commits/push con secretos potenciales."
    )
    parser.add_argument(
        "--mode",
        choices=["staged", "tracked"],
        default="staged",
        help="staged para pre-commit, tracked para pre-push.",
    )
    args = parser.parse_args()

    try:
        files = list_staged_files() if args.mode == "staged" else list_tracked_files()
    except RuntimeError as exc:
        print(f"[secret-guard] {exc}", file=sys.stderr)
        return 2

    if not files:
        return 0

    root = Path.cwd()
    all_issues = []
    for rel in files:
        path = root / rel
        issues = scan_file(path)
        if issues:
            all_issues.append((rel, issues))

    if not all_issues:
        return 0

    print("[secret-guard] Se detectaron posibles secretos. Commit/push bloqueado.\n")
    for rel, issues in all_issues:
        print(f"- {rel}")
        for line_no, msg in issues:
            print(f"  L{line_no}: {msg}")

    print(
        "\nAcciones recomendadas:\n"
        "1) Sustituir el valor por placeholder o variable de entorno.\n"
        "2) Si es un falso positivo, añade comentario 'secret-scan: allow' en esa línea."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
