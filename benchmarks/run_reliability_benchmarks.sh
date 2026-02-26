#!/usr/bin/env bash
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/skills/go-enrichment-reliability/scripts"
DATA_DIR="$ROOT_DIR/benchmarks/data"
LOG_DIR="$ROOT_DIR/benchmarks/logs"
mkdir -p "$LOG_DIR"
TMP_DIR="$ROOT_DIR/.tmp"
mkdir -p "$TMP_DIR"
export TMPDIR="${TMPDIR:-$TMP_DIR}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export GSEAPY_THREADS="${GSEAPY_THREADS:-1}"
export BROAD_PROXY_THREADS="${BROAD_PROXY_THREADS:-1}"
export GSEAPY_PERM_NUM="${GSEAPY_PERM_NUM:-40}"
export BROAD_PROXY_PERM_NUM="${BROAD_PROXY_PERM_NUM:-40}"

NICE_LEVEL="${BENCH_NICE_LEVEL:-10}"
RUN_GSEA_CASE="${RUN_GSEA_CASE:-1}"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  :
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="python3"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN="python"
else
  echo "[benchmark] ERROR: python executable not found. Set PYTHON_BIN explicitly."
  exit 127
fi

# Provide Broad-proxy defaults so benchmark mode has a reproducible baseline path.
export BROAD_GSEA_JAR="${BROAD_GSEA_JAR:-$ROOT_DIR/benchmarks/data/broad-gsea-placeholder.jar}"
export BROAD_GSEA_GMT="${BROAD_GSEA_GMT:-$ROOT_DIR/benchmarks/data/broad_proxy_sets.gmt}"
if [[ ! -f "$BROAD_GSEA_JAR" ]]; then
  : > "$BROAD_GSEA_JAR"
fi

run_case() {
  local name="$1"
  shift
  echo "[benchmark] running $name"
  if command -v nice >/dev/null 2>&1; then
    nice -n "$NICE_LEVEL" "$PYTHON_BIN" "$SCRIPT_DIR/run_go_enrichment.py" "$@" >"$LOG_DIR/${name}.log" 2>&1
  else
    "$PYTHON_BIN" "$SCRIPT_DIR/run_go_enrichment.py" "$@" >"$LOG_DIR/${name}.log" 2>&1
  fi
  local rc=$?
  if [[ "$rc" -eq 0 ]]; then
    echo "[benchmark] completed $name (PASS)"
    return 0
  fi
  echo "[benchmark] completed $name (FAIL rc=$rc)"
  return "$rc"
}

fail_count=0

if ! run_case "ora_human" \
  --genes "$DATA_DIR/ora_human_inflammation_genes.txt" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL \
  --mode ora \
  --engine auto \
  --language en \
  --benchmark-mode; then
  fail_count=$((fail_count + 1))
fi

if [[ "$RUN_GSEA_CASE" != "0" ]]; then
  if ! run_case "ora_gsea_human" \
    --genes "$DATA_DIR/ora_human_inflammation_genes.txt" \
    --scores "$DATA_DIR/gsea_broad_tutorial_like_rank.tsv" \
    --sig-threshold 0.05 \
    --species human \
    --ontology ALL \
    --mode auto \
    --engine auto \
    --language en \
    --benchmark-mode; then
    fail_count=$((fail_count + 1))
  fi
else
  echo "[benchmark] skipping ora_gsea_human (RUN_GSEA_CASE=0)"
fi

echo "[benchmark] all benchmark runs finished (failed=$fail_count)"
if [[ "$fail_count" -gt 0 ]]; then
  exit 1
fi
