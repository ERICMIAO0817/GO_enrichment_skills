#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SKILL_NAME="go-enrichment-reliability"
SRC_DIR="$ROOT_DIR/skills/$SKILL_NAME"

MODE="copy"
FORCE=false
DEST_SKILLS_DIR=""

usage() {
  cat <<'EOF'
Install GO enrichment skill into Codex skills directory.

Usage:
  bash install.sh [options]
  bash scripts/install_skill.sh [options]

Options:
  -d, --dest <dir>   Target skills directory. Default: $CODEX_HOME/skills or ~/.codex/skills
  --copy             Copy skill directory (default)
  --link             Symlink skill directory (useful for local development)
  --force            Overwrite existing installed skill
  -h, --help         Show this help message

Examples:
  bash install.sh
  bash install.sh --force
  bash install.sh --link --force
  bash install.sh --dest ~/.codex/skills --force
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
    --copy)
      MODE="copy"
      shift
      ;;
    --link)
      MODE="link"
      shift
      ;;
    --force)
      FORCE=true
      shift
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

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Error: skill source not found: $SRC_DIR" >&2
  exit 1
fi

if [[ -z "$DEST_SKILLS_DIR" ]]; then
  DEST_SKILLS_DIR="$(detect_default_skills_dir)"
fi

DEST_SKILLS_DIR="${DEST_SKILLS_DIR/#\~/$HOME}"
DEST_DIR="$DEST_SKILLS_DIR/$SKILL_NAME"

mkdir -p "$DEST_SKILLS_DIR"

if [[ -e "$DEST_DIR" || -L "$DEST_DIR" ]]; then
  if [[ "$FORCE" != "true" ]]; then
    echo "Skill already exists: $DEST_DIR" >&2
    echo "Re-run with --force to overwrite." >&2
    exit 1
  fi
  rm -rf "$DEST_DIR"
fi

if [[ "$MODE" == "link" ]]; then
  ln -s "$SRC_DIR" "$DEST_DIR"
else
  cp -R "$SRC_DIR" "$DEST_DIR"
fi

if [[ ! -f "$DEST_DIR/SKILL.md" ]]; then
  echo "Error: install failed, SKILL.md missing at $DEST_DIR" >&2
  exit 1
fi

echo "Installed skill: $SKILL_NAME"
echo "Location: $DEST_DIR"
if [[ "$MODE" == "link" ]]; then
  echo "Mode: symlink"
else
  echo "Mode: copy"
fi
