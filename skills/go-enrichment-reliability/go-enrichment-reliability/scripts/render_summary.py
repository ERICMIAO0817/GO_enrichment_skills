#!/usr/bin/env python3
"""Render markdown summary for a GO/GSEA run in zh/en."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import List


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Render markdown run summary")
    p.add_argument("--outdir", required=True)
    p.add_argument("--language", choices=["zh", "en"], default="en")
    p.add_argument("--fallback-note", default="")
    p.add_argument("--warning", action="append", default=[])
    return p


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _fmt_pct(v) -> str:
    if v is None:
        return "N/A"
    try:
        return f"{float(v) * 100:.1f}%"
    except Exception:
        return "N/A"


def _read_csv_rows(path: Path) -> List[dict]:
    if not path.exists():
        return []
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    except Exception:
        return []


def _top_celltypes(outdir: Path, topk: int = 3) -> List[dict]:
    rows = _read_csv_rows(outdir / "celltype_association.csv")
    if not rows:
        return []
    rows.sort(key=lambda r: float(r.get("combined_score") or 0.0), reverse=True)
    return rows[:topk]


def _basis_zh(v: str) -> str:
    mapping = {
        "mixed_marker_functional": "marker+功能双证据",
        "marker_driven": "marker主导",
        "functional_driven": "功能术语主导",
        "weak_signal": "弱信号",
    }
    return mapping.get(v or "", v or "unknown")


def _call_zh(v: str) -> str:
    mapping = {
        "associated": "相关",
        "likely_associated": "较相关",
        "possible_associated": "可能相关",
        "functional_hint_only": "功能提示相关",
        "weak_or_uncertain": "弱相关/不确定",
        "not_supported": "暂不支持",
    }
    return mapping.get(v or "", v or "未知")


def _call_en(v: str) -> str:
    mapping = {
        "associated": "associated",
        "likely_associated": "likely associated",
        "possible_associated": "possibly associated",
        "functional_hint_only": "functional hint only",
        "weak_or_uncertain": "weak/uncertain",
        "not_supported": "not supported",
    }
    return mapping.get(v or "", v or "unknown")


def render_zh(outdir: Path, meta: dict, consistency: dict, warnings: List[str], fallback_note: str) -> str:
    ora = consistency.get("ora", {})
    gsea = consistency.get("gsea", {})

    lines = [
        "# GO/GSEA 富集结果摘要",
        "",
        "## 运行概览",
        f"- 输出目录: `{outdir}`",
        f"- 物种: `{meta.get('species', 'unknown')}`",
        f"- ontology: `{','.join(meta.get('ontology', [])) if isinstance(meta.get('ontology'), list) else meta.get('ontology', 'ALL')}`",
        f"- 显著性阈值(FDR): `{meta.get('sig_threshold', 'N/A')}`",
        f"- 主执行引擎(ORA/GSEA): `{meta.get('run_sources', {}).get('ora', meta.get('engine', 'unknown'))}` / `{meta.get('run_sources', {}).get('gsea', meta.get('engine', 'unknown'))}`",
        "",
        "## 一致性评估",
        f"- ORA Top20 术语重叠率: `{_fmt_pct(ora.get('score'))}` (阈值: 60%)",
        f"- GSEA Top20 方向一致率: `{_fmt_pct(gsea.get('score'))}` (阈值: 60%)",
        f"- 综合判定: `{'通过' if consistency.get('overall_pass') else '未通过'}`",
        "",
    ]

    top_celltypes = _top_celltypes(outdir)
    lines.append("## 细胞类型关联（端到端）")
    if top_celltypes:
        for row in top_celltypes:
            basis = _basis_zh(str(row.get("evidence_basis", "")))
            qcall = _call_zh(str(row.get("qualitative_call", "")))
            lines.append(
                f"- Top{row.get('rank', '?')}: `{row.get('cell_type_name', 'N/A')}` -> 定性判断: `{qcall}` "
                f"(basis={basis}, marker={row.get('marker_hits', '0')}/{row.get('marker_total', '0')}, "
                f"keyword={row.get('keyword_hits', '0')}/{row.get('keyword_total', '0')}, ORA/GSEA={row.get('ora_term_hits', '0')}/{row.get('gsea_term_hits', '0')})"
            )
            if str(row.get("interpretation", "")).strip():
                lines.append(f"  解读: {row.get('interpretation')}")
            if str(row.get("validation_suggestion", "")).strip():
                lines.append(f"  建议: {row.get('validation_suggestion')}")
    else:
        lines.append("- 未检出稳定的细胞类型关联信号。")
    lines.append("- 详见: `celltype_association.csv`、`celltype_association.json`、`celltype_interpretation.md`、`analysis_brief.md`")
    lines.append("")

    if fallback_note:
        lines.extend(["## 降级说明", f"- {fallback_note}", ""])

    if warnings:
        lines.append("## 告警")
        lines.extend([f"- {w}" for w in warnings])
        lines.append("")

    lines.extend(
        [
            "## 输出文件",
            "- `summary.md`",
            "- `normalized_input.tsv`",
            "- `ora_results_python.csv`",
            "- `ora_results_r.csv`",
            "- `gsea_results_python.csv`",
            "- `gsea_results_r.csv`",
            "- `baseline_enrichr.csv`",
            "- `baseline_gprofiler.csv`",
            "- `baseline_clusterprofiler.csv`",
            "- `baseline_broad_gsea.csv`",
            "- `consistency_report.json`",
            "- `celltype_association.csv`",
            "- `celltype_association.json`",
            "- `celltype_interpretation.md`",
            "- `analysis_brief.md`",
            "- `analysis_brief.json`",
            "- `plots_dotplot.png`",
            "- `plots_barplot.png`",
            "",
        ]
    )

    return "\n".join(lines).rstrip() + "\n"


def render_en(outdir: Path, meta: dict, consistency: dict, warnings: List[str], fallback_note: str) -> str:
    ora = consistency.get("ora", {})
    gsea = consistency.get("gsea", {})

    lines = [
        "# GO/GSEA Enrichment Summary",
        "",
        "## Run Overview",
        f"- Output directory: `{outdir}`",
        f"- Species: `{meta.get('species', 'unknown')}`",
        f"- Ontology: `{','.join(meta.get('ontology', [])) if isinstance(meta.get('ontology'), list) else meta.get('ontology', 'ALL')}`",
        f"- Significance threshold (FDR): `{meta.get('sig_threshold', 'N/A')}`",
        "",
        "## Consistency Metrics",
        f"- ORA Top20 overlap: `{_fmt_pct(ora.get('score'))}` (threshold: 60%)",
        f"- GSEA Top20 direction consistency: `{_fmt_pct(gsea.get('score'))}` (threshold: 60%)",
        f"- Overall gate: `{'PASS' if consistency.get('overall_pass') else 'FAIL'}`",
        "",
    ]

    top_celltypes = _top_celltypes(outdir)
    lines.append("## Cell-Type Association (End-to-End)")
    if top_celltypes:
        for row in top_celltypes:
            qcall = _call_en(str(row.get("qualitative_call", "")))
            lines.append(
                f"- Top{row.get('rank', '?')}: `{row.get('cell_type_name', 'N/A')}` -> qualitative call: `{qcall}` "
                f"(basis={row.get('evidence_basis', 'unknown')}, marker={row.get('marker_hits', '0')}/{row.get('marker_total', '0')}, "
                f"keyword={row.get('keyword_hits', '0')}/{row.get('keyword_total', '0')}, ORA/GSEA={row.get('ora_term_hits', '0')}/{row.get('gsea_term_hits', '0')})"
            )
            if str(row.get("interpretation", "")).strip():
                lines.append(f"  Interpretation: {row.get('interpretation')}")
            if str(row.get("validation_suggestion", "")).strip():
                lines.append(f"  Suggestion: {row.get('validation_suggestion')}")
    else:
        lines.append("- No stable cell-type signal detected from current genes/terms.")
    lines.append("- See `celltype_association.csv`, `celltype_association.json`, `celltype_interpretation.md`, and `analysis_brief.md`.")
    lines.append("")

    if fallback_note:
        lines.extend(["## Fallback Note", f"- {fallback_note}", ""])

    if warnings:
        lines.append("## Warnings")
        lines.extend([f"- {w}" for w in warnings])
        lines.append("")

    lines.extend(
        [
            "## Output Files",
            "- `summary.md`",
            "- `normalized_input.tsv`",
            "- `ora_results_python.csv`",
            "- `ora_results_r.csv`",
            "- `gsea_results_python.csv`",
            "- `gsea_results_r.csv`",
            "- `baseline_enrichr.csv`",
            "- `baseline_gprofiler.csv`",
            "- `baseline_clusterprofiler.csv`",
            "- `baseline_broad_gsea.csv`",
            "- `consistency_report.json`",
            "- `celltype_association.csv`",
            "- `celltype_association.json`",
            "- `celltype_interpretation.md`",
            "- `analysis_brief.md`",
            "- `analysis_brief.json`",
            "- `plots_dotplot.png`",
            "- `plots_barplot.png`",
            "",
        ]
    )

    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir).expanduser().resolve()

    meta = _read_json(outdir / "run_meta.json")
    consistency = _read_json(outdir / "consistency_report.json")

    warnings = list(args.warning)
    if isinstance(meta.get("warnings"), list):
        warnings.extend([str(x) for x in meta["warnings"]])
    # De-duplicate while preserving order.
    deduped = []
    seen = set()
    for w in warnings:
        if w not in seen:
            seen.add(w)
            deduped.append(w)
    warnings = deduped

    text = (
        render_zh(outdir, meta, consistency, warnings, args.fallback_note)
        if args.language == "zh"
        else render_en(outdir, meta, consistency, warnings, args.fallback_note)
    )
    (outdir / "summary.md").write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
