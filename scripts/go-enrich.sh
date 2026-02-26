#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENTRY="$ROOT_DIR/skills/go-enrichment-reliability/scripts/run_go_enrichment.py"
REQUEST_PARSER="$ROOT_DIR/skills/go-enrichment-reliability/scripts/parse_request_text.py"

if [[ ! -f "$ENTRY" ]]; then
  echo "Error: entry script not found: $ENTRY" >&2
  exit 1
fi

GENES=""
GENES_SET=false
SIG_THRESHOLD="0.05"
SIG_THRESHOLD_SET=false
SPECIES="human"
SPECIES_SET=false
ONTOLOGY="ALL"
ONTOLOGY_SET=false
SCORES=""
MODE="auto"
ENGINE="auto"
BACKGROUND=""
LANGUAGE="auto"
OUTDIR=""
FOCUS=""
FOCUS_SET=false
REQUEST_TEXT=""
BENCHMARK_MODE=false
SKIP_BASELINES=false
SKIP_CELLTYPE_LINK=false

usage() {
  cat <<'USAGE'
Conversation-friendly GO enrichment CLI wrapper.

Usage:
  bash scripts/go-enrich.sh --genes "TP53,EGFR" [options]
  bash scripts/go-enrich.sh --request "用这组基因[...]，sig_threshold=0.05，species=mouse，ontology=ALL，重点看beta cell"

Required:
  one of:
    --genes <inline_list_or_file>
    --request <natural_language_text>

Options:
  --request <free_text>            parse genes/sig/species/ontology/focus from natural language
  --sig-threshold <float>         default: 0.05
  --species <human|mouse>         default: human
  --ontology <ALL|BP,MF,CC>       default: ALL
  --scores <file_or_inline_rank>  optional, enables GSEA in auto mode
  --mode <auto|ora|gsea>          default: auto
  --engine <auto|python|r>        default: auto
  --background <file>             optional
  --language <zh|en|auto>         default: auto
  --outdir <path>                 optional
  --focus <comma_separated_cells> maps to --celltype-candidates
  --benchmark-mode                strict reliability gate
  --skip-baselines                skip baseline collection
  --skip-celltype-link            skip cell-type inference
  -h, --help                      show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genes) GENES="$2"; GENES_SET=true; shift 2 ;;
    --request) REQUEST_TEXT="$2"; shift 2 ;;
    --sig-threshold) SIG_THRESHOLD="$2"; SIG_THRESHOLD_SET=true; shift 2 ;;
    --species) SPECIES="$2"; SPECIES_SET=true; shift 2 ;;
    --ontology) ONTOLOGY="$2"; ONTOLOGY_SET=true; shift 2 ;;
    --scores) SCORES="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --engine) ENGINE="$2"; shift 2 ;;
    --background) BACKGROUND="$2"; shift 2 ;;
    --language) LANGUAGE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --focus) FOCUS="$2"; FOCUS_SET=true; shift 2 ;;
    --benchmark-mode) BENCHMARK_MODE=true; shift ;;
    --skip-baselines) SKIP_BASELINES=true; shift ;;
    --skip-celltype-link) SKIP_CELLTYPE_LINK=true; shift ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -n "$REQUEST_TEXT" ]]; then
  if [[ ! -f "$REQUEST_PARSER" ]]; then
    echo "Error: request parser not found: $REQUEST_PARSER" >&2
    exit 1
  fi
  set +e
  PARSED_VARS="$(python3 "$REQUEST_PARSER" --request "$REQUEST_TEXT" --format shell)"
  PARSE_RC=$?
  set -e
  if [[ $PARSE_RC -ne 0 ]]; then
    echo "Warning: --request parsing incomplete; explicit CLI args will be preferred." >&2
  fi
  if [[ -n "$PARSED_VARS" ]]; then
    eval "$PARSED_VARS"
  fi
  if [[ "$GENES_SET" == false ]]; then
    GENES="${PARSED_GENES:-$GENES}"
  fi
  if [[ "$SIG_THRESHOLD_SET" == false && -n "${PARSED_SIG_THRESHOLD:-}" ]]; then
    SIG_THRESHOLD="$PARSED_SIG_THRESHOLD"
  fi
  if [[ "$SPECIES_SET" == false && -n "${PARSED_SPECIES:-}" ]]; then
    SPECIES="$PARSED_SPECIES"
  fi
  if [[ "$ONTOLOGY_SET" == false && -n "${PARSED_ONTOLOGY:-}" ]]; then
    ONTOLOGY="$PARSED_ONTOLOGY"
  fi
  if [[ "$FOCUS_SET" == false && -n "${PARSED_FOCUS:-}" ]]; then
    FOCUS="$PARSED_FOCUS"
  fi
fi

if [[ -z "$GENES" ]]; then
  echo "Error: genes are required. Provide --genes or --request." >&2
  usage
  exit 2
fi

CMD=(
  python3 "$ENTRY"
  --genes "$GENES"
  --sig-threshold "$SIG_THRESHOLD"
  --species "$SPECIES"
  --ontology "$ONTOLOGY"
  --mode "$MODE"
  --engine "$ENGINE"
  --language "$LANGUAGE"
)

if [[ -n "$SCORES" ]]; then
  CMD+=(--scores "$SCORES")
fi
if [[ -n "$BACKGROUND" ]]; then
  CMD+=(--background "$BACKGROUND")
fi
if [[ -n "$OUTDIR" ]]; then
  CMD+=(--outdir "$OUTDIR")
fi
if [[ -n "$FOCUS" ]]; then
  CMD+=(--celltype-candidates "$FOCUS")
fi
if [[ "$BENCHMARK_MODE" == true ]]; then
  CMD+=(--benchmark-mode)
fi
if [[ "$SKIP_BASELINES" == true ]]; then
  CMD+=(--skip-baselines)
fi
if [[ "$SKIP_CELLTYPE_LINK" == true ]]; then
  CMD+=(--skip-celltype-link)
fi

"${CMD[@]}"
