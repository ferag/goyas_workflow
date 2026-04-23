#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
HOOKS_DIR="$REPO_ROOT/.git/hooks"
GUARD_SCRIPT="$REPO_ROOT/scripts/security/secret_guard.py"

if [[ ! -f "$GUARD_SCRIPT" ]]; then
  echo "No se encuentra $GUARD_SCRIPT"
  exit 1
fi

cat > "$HOOKS_DIR/pre-commit" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
python3 scripts/security/secret_guard.py --mode staged
EOF

cat > "$HOOKS_DIR/pre-push" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
python3 scripts/security/secret_guard.py --mode tracked
EOF

chmod +x "$HOOKS_DIR/pre-commit" "$HOOKS_DIR/pre-push"

echo "Hooks instalados:"
echo " - .git/hooks/pre-commit"
echo " - .git/hooks/pre-push"
