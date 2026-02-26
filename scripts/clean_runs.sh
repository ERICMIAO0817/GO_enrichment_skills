#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_DIRS=(
  "$ROOT_DIR/runs"
  "$ROOT_DIR/skills/go-enrichment-reliability/scripts/runs"
)

usage() {
  cat <<'EOF'
Usage:
  bash scripts/clean_runs.sh [--dry-run]

Options:
  --dry-run   Print what would be removed, but do not delete anything.
EOF
}

DRY_RUN=false
if [[ $# -gt 1 ]]; then
  usage
  exit 2
fi
if [[ $# -eq 1 ]]; then
  case "$1" in
    --dry-run) DRY_RUN=true ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage
      exit 2
      ;;
  esac
fi

for dir in "${TARGET_DIRS[@]}"; do
  mkdir -p "$dir"
  echo "Target: $dir"
  if [[ "$DRY_RUN" == "true" ]]; then
    find "$dir" -mindepth 1 -maxdepth 1 -print | sed 's#^#  would remove: #'
    continue
  fi

  find "$dir" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  echo "  cleared"
done

if [[ "$DRY_RUN" == "true" ]]; then
  echo "Dry run complete."
else
  echo "Run directories cleared."
fi
