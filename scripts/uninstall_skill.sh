#!/usr/bin/env bash
set -euo pipefail

SKILL_NAME="go-enrichment-reliability"
DEST_SKILLS_DIR=""

usage() {
  cat <<'EOF'
Uninstall GO enrichment skill from Codex skills directory.

Usage:
  bash uninstall.sh [options]
  bash scripts/uninstall_skill.sh [options]

Options:
  -d, --dest <dir>   Target skills directory. Default: $CODEX_HOME/skills or ~/.codex/skills
  -h, --help         Show this help message

Examples:
  bash uninstall.sh
  bash uninstall.sh --dest ~/.codex/skills
EOF
}

detect_default_skills_dir() {
  if [[ -n "${CODEX_HOME:-}" ]]; then
    echo "$CODEX_HOME/skills"
    return
  fi
  echo "$HOME/.codex/skills"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -d|--dest)
      if [[ $# -lt 2 ]]; then
        echo "Error: --dest requires a path." >&2
        exit 2
      fi
      DEST_SKILLS_DIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown option '$1'." >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "$DEST_SKILLS_DIR" ]]; then
  DEST_SKILLS_DIR="$(detect_default_skills_dir)"
fi

DEST_SKILLS_DIR="${DEST_SKILLS_DIR/#\~/$HOME}"
TARGET="$DEST_SKILLS_DIR/$SKILL_NAME"

if [[ ! -e "$TARGET" && ! -L "$TARGET" ]]; then
  echo "Skill not found: $TARGET"
  exit 0
fi

rm -rf "$TARGET"
echo "Removed skill: $TARGET"
