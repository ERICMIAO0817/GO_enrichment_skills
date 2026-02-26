#!/usr/bin/env python3
"""Compute consistency metrics between run outputs and baselines."""

from __future__ import annotations

import argparse
from collections import defaultdict
import re
import statistics
from pathlib import Path
from typing import Dict, List, Tuple

from enrichment_common import read_csv, write_json


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Compare run results against baseline tools")
    p.add_argument("--outdir", required=True)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--ora-threshold", type=float, default=0.6)
    p.add_argument("--gsea-threshold", type=float, default=0.6)
    p.add_argument("--benchmark-mode", default="false")
    return p


def parse_bool(v: str) -> bool:
    return str(v).lower() in {"1", "true", "yes", "y"}


_GO_ID_RE = re.compile(r"\bGO:\d{7}\b", flags=re.IGNORECASE)


def _extract_term_id(row: dict) -> str:
    tid = (row.get("term_id") or "").strip()
    if tid:
        m = _GO_ID_RE.search(tid)
        return m.group(0).upper() if m else tid.lower()
    tname = (row.get("term_name") or "").strip()
    m = _GO_ID_RE.search(tname)
    return m.group(0).upper() if m else ""


def _normalize_term_name(text: str) -> str:
    text = _GO_ID_RE.sub(" ", text or "")
    text = re.sub(r"[\W_]+", " ", text.lower())
    return re.sub(r"\s+", " ", text).strip()


def _normalize_term(row: dict) -> str:
    tid = _extract_term_id(row)
    tname = _normalize_term_name((row.get("term_name") or "").strip())
    if tid:
        return tid
    return tname


def _safe_float(v, default: float = 1.0) -> float:
    try:
        return float(v)
    except Exception:
        return default


def _top_terms_ora(rows: List[dict], top_n: int) -> List[dict]:
    sorted_rows = sorted(rows, key=lambda r: (_safe_float(r.get("padj"), 1.0), _safe_float(r.get("pvalue"), 1.0)))
    out = []
    seen = set()
    for r in sorted_rows:
        key = _normalize_term(r)
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(r)
        if len(out) >= top_n:
            break
    return out


def _top_terms_gsea(rows: List[dict], top_n: int) -> List[dict]:
    sorted_rows = sorted(rows, key=lambda r: (_safe_float(r.get("padj"), 1.0), abs(_safe_float(r.get("nes"), 0.0)) * -1))
    out = []
    seen = set()
    for r in sorted_rows:
        key = _normalize_term(r)
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(r)
        if len(out) >= top_n:
            break
    return out


def _direction(row: dict) -> str:
    d = str(row.get("direction", "")).strip().lower()
    if d in {"up", "down"}:
        return d
    nes = _safe_float(row.get("nes"), 0.0)
    return "up" if nes >= 0 else "down"


def compute_ora_overlap(run_rows: List[dict], baseline_rows: List[dict], top_n: int) -> Dict:
    run_top = _top_terms_ora(run_rows, top_n)
    base_top = _top_terms_ora(baseline_rows, top_n)
    run_set = {_normalize_term(r) for r in run_top}
    base_set = {_normalize_term(r) for r in base_top}

    if not run_set:
        return {"overlap": 0.0, "common": 0, "denominator": 0}

    common = run_set & base_set
    denom = len(run_set)
    overlap = len(common) / denom if denom else 0.0
    return {"overlap": overlap, "common": len(common), "denominator": denom}


def build_ora_consensus(baseline_map: Dict[str, List[dict]], top_n: int) -> List[str]:
    rank_lists = {name: _top_terms_ora(rows, top_n) for name, rows in baseline_map.items()}
    if not any(rank_lists.values()):
        return []

    min_rank = {}
    support = defaultdict(int)
    mean_rank = defaultdict(list)
    for rows in rank_lists.values():
        for idx, row in enumerate(rows, start=1):
            key = _normalize_term(row)
            if not key:
                continue
            support[key] += 1
            mean_rank[key].append(idx)
            prev = min_rank.get(key, 10**9)
            if idx < prev:
                min_rank[key] = idx

    ranked = sorted(
        min_rank.keys(),
        key=lambda k: (min_rank.get(k, 10**9), -support.get(k, 0), statistics.mean(mean_rank.get(k, [10**9]))),
    )
    return ranked[:top_n]


def compute_ora_overlap_to_consensus(run_rows: List[dict], consensus_terms: List[str], top_n: int) -> Dict:
    run_top = _top_terms_ora(run_rows, top_n)
    run_set = {_normalize_term(r) for r in run_top}
    run_set = {x for x in run_set if x}
    if not run_set:
        return {"overlap": 0.0, "common": 0, "denominator": 0}

    base_set = {x for x in consensus_terms if x}
    if not base_set:
        return {"overlap": 0.0, "common": 0, "denominator": len(run_set)}

    common = run_set & base_set
    denom = len(run_set)
    return {"overlap": len(common) / denom if denom else 0.0, "common": len(common), "denominator": denom}


def compute_gsea_direction_consistency(run_rows: List[dict], baseline_rows: List[dict], top_n: int) -> Dict:
    run_top = _top_terms_gsea(run_rows, top_n)
    base_top = _top_terms_gsea(baseline_rows, top_n)

    run_map = {_normalize_term(r): _direction(r) for r in run_top}
    base_map = {_normalize_term(r): _direction(r) for r in base_top}
    common_keys = [k for k in run_map.keys() if k in base_map]

    if not common_keys:
        return {"consistency": 0.0, "common": 0, "denominator": 0}

    same = sum(1 for k in common_keys if run_map[k] == base_map[k])
    denom = len(common_keys)
    return {
        "consistency": same / denom if denom else 0.0,
        "common": same,
        "denominator": denom,
    }


def build_gsea_consensus_direction(baseline_map: Dict[str, List[dict]], top_n: int) -> Dict[str, str]:
    direction_votes: Dict[str, List[str]] = defaultdict(list)
    for rows in baseline_map.values():
        top = _top_terms_gsea(rows, top_n)
        for row in top:
            key = _normalize_term(row)
            if not key:
                continue
            direction_votes[key].append(_direction(row))

    out: Dict[str, str] = {}
    for key, dirs in direction_votes.items():
        if not dirs:
            continue
        up_n = sum(1 for d in dirs if d == "up")
        down_n = sum(1 for d in dirs if d == "down")
        out[key] = "up" if up_n >= down_n else "down"
    return out


def compute_gsea_consistency_to_consensus(run_rows: List[dict], consensus_map: Dict[str, str], top_n: int) -> Dict:
    run_top = _top_terms_gsea(run_rows, top_n)
    run_map = {_normalize_term(r): _direction(r) for r in run_top}
    run_map = {k: v for k, v in run_map.items() if k}

    common_keys = [k for k in run_map.keys() if k in consensus_map]
    if not common_keys:
        return {"consistency": 0.0, "common": 0, "denominator": 0}

    same = sum(1 for k in common_keys if run_map[k] == consensus_map[k])
    denom = len(common_keys)
    return {"consistency": same / denom if denom else 0.0, "common": same, "denominator": denom}


def pick_primary(run_py: List[dict], run_r: List[dict]) -> Tuple[str, List[dict]]:
    if run_py:
        return "python", run_py
    if run_r:
        return "r", run_r
    return "none", []


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir).expanduser().resolve()
    benchmark_mode = parse_bool(args.benchmark_mode)

    run_ora_py = read_csv(outdir / "ora_results_python.csv")
    run_ora_r = read_csv(outdir / "ora_results_r.csv")
    run_gsea_py = read_csv(outdir / "gsea_results_python.csv")
    run_gsea_r = read_csv(outdir / "gsea_results_r.csv")

    ora_engine, run_ora = pick_primary(run_ora_py, run_ora_r)
    gsea_engine, run_gsea = pick_primary(run_gsea_py, run_gsea_r)

    ora_baselines = {
        "enrichr": read_csv(outdir / "baseline_enrichr.csv"),
        "gprofiler": read_csv(outdir / "baseline_gprofiler.csv"),
        "clusterprofiler": read_csv(outdir / "baseline_clusterprofiler.csv"),
    }

    gsea_baselines = {
        "clusterprofiler_gsea": read_csv(outdir / "baseline_clusterprofiler_gsea.csv"),
        "broad_gsea": read_csv(outdir / "baseline_broad_gsea.csv"),
    }

    ora_details = {}
    for name, rows in ora_baselines.items():
        metric = compute_ora_overlap(run_ora, rows, args.top_n)
        ora_details[name] = metric
    ora_consensus_terms = build_ora_consensus(ora_baselines, args.top_n)
    ora_consensus_metric = compute_ora_overlap_to_consensus(run_ora, ora_consensus_terms, args.top_n)
    ora_details["consensus"] = ora_consensus_metric

    ora_score = ora_consensus_metric["overlap"] if ora_consensus_metric["denominator"] > 0 else 0.0
    ora_pass = ora_score >= args.ora_threshold if run_ora else False

    gsea_details = {}
    for name, rows in gsea_baselines.items():
        metric = compute_gsea_direction_consistency(run_gsea, rows, args.top_n)
        gsea_details[name] = metric
    gsea_consensus_map = build_gsea_consensus_direction(gsea_baselines, args.top_n)
    gsea_consensus_metric = compute_gsea_consistency_to_consensus(run_gsea, gsea_consensus_map, args.top_n)
    gsea_details["consensus"] = gsea_consensus_metric

    gsea_applicable = len(run_gsea) > 0
    if gsea_applicable:
        gsea_score = gsea_consensus_metric["consistency"] if gsea_consensus_metric["denominator"] > 0 else 0.0
        gsea_pass = gsea_score >= args.gsea_threshold
    else:
        gsea_score = None
        gsea_pass = True

    report = {
        "top_n": args.top_n,
        "thresholds": {
            "ora_overlap": args.ora_threshold,
            "gsea_direction_consistency": args.gsea_threshold,
        },
        "run_sources": {
            "ora": ora_engine,
            "gsea": gsea_engine,
        },
        "ora": {
            "score": ora_score,
            "pass": ora_pass,
            "details": ora_details,
        },
        "gsea": {
            "applicable": gsea_applicable,
            "score": gsea_score,
            "pass": gsea_pass,
            "details": gsea_details,
        },
        "overall_pass": ora_pass and gsea_pass,
    }

    write_json(outdir / "consistency_report.json", report)

    if benchmark_mode and not report["overall_pass"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
