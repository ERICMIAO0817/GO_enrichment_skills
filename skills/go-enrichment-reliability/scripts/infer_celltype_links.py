#!/usr/bin/env python3
"""Infer likely cell-type associations from normalized genes and enrichment terms."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Sequence

from enrichment_common import read_csv, resolve_language, write_csv, write_json

CELLTYPE_FIELDS = [
    "rank",
    "cell_type_id",
    "cell_type_name",
    "qualitative_call",
    "combined_score",
    "evidence_basis",
    "marker_hits",
    "marker_total",
    "marker_hit_ratio",
    "keyword_hits",
    "keyword_total",
    "keyword_hit_ratio",
    "ora_term_hits",
    "gsea_term_hits",
    "gsea_up_hits",
    "gsea_down_hits",
    "confidence",
    "evidence_genes",
    "missing_markers",
    "evidence_keywords",
    "top_ora_terms",
    "top_gsea_terms",
    "evidence_terms",
    "interpretation",
    "validation_suggestion",
]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Infer cell-type association from enrichment outputs")
    p.add_argument("--outdir", required=True)
    p.add_argument("--species", choices=["human", "mouse"], default="human")
    p.add_argument("--candidates", default="")
    p.add_argument("--topk", type=int, default=5)
    p.add_argument("--language", choices=["zh", "en", "auto"], default="auto")
    p.add_argument("--context-text", default="")
    return p


def _norm_text(s: str) -> str:
    return re.sub(r"[\W_]+", " ", (s or "").lower()).strip()


def _safe_float(v, default: float = 1.0) -> float:
    try:
        return float(v)
    except Exception:
        return default


def _fmt_p(v) -> str:
    x = _safe_float(v, float("nan"))
    if math.isnan(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def _pick_primary(py_rows: List[dict], r_rows: List[dict]) -> List[dict]:
    if py_rows:
        return py_rows
    if r_rows:
        return r_rows
    return []


def _top_terms_ora(rows: List[dict], top_n: int) -> List[dict]:
    sorted_rows = sorted(rows, key=lambda r: (_safe_float(r.get("padj"), 1.0), _safe_float(r.get("pvalue"), 1.0)))
    out = []
    seen = set()
    for r in sorted_rows:
        key = (r.get("term_id") or "").strip() or (r.get("term_name") or "").strip().lower()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(r)
        if len(out) >= top_n:
            break
    return out


def _top_terms_gsea(rows: List[dict], top_n: int) -> List[dict]:
    sorted_rows = sorted(rows, key=lambda r: (_safe_float(r.get("padj"), 1.0), -abs(_safe_float(r.get("nes"), 0.0))))
    out = []
    seen = set()
    for r in sorted_rows:
        key = (r.get("term_id") or "").strip() or (r.get("term_name") or "").strip().lower()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(r)
        if len(out) >= top_n:
            break
    return out


def _load_reference(script_dir: Path, species: str) -> Dict[str, dict]:
    ref_path = script_dir.parent / "references" / "celltype_markers.json"
    payload = json.loads(ref_path.read_text(encoding="utf-8"))
    species_map = payload.get(species, {})
    return species_map if isinstance(species_map, dict) else {}


def _parse_candidates(raw: str) -> List[str]:
    if not raw:
        return []
    return [_norm_text(x) for x in re.split(r"[,;，]+", raw) if _norm_text(x)]


def _select_celltypes(celltypes: Dict[str, dict], candidates: Sequence[str]) -> Dict[str, dict]:
    if not candidates:
        return celltypes

    def _match(candidate: str, name: str) -> bool:
        if not candidate or not name:
            return False
        if candidate == name:
            return True
        c_words = candidate.split()
        n_words = name.split()
        if not c_words or not n_words:
            return False
        c_set = set(c_words)
        n_set = set(n_words)
        return c_set.issubset(n_set) or n_set.issubset(c_set)

    selected: Dict[str, dict] = {}
    for cell_id, spec in celltypes.items():
        names = [cell_id, str(spec.get("display_name", ""))]
        names.extend(spec.get("aliases", []))
        norm_names = [_norm_text(n) for n in names if n]
        for c in candidates:
            if any(_match(c, n) for n in norm_names):
                selected[cell_id] = spec
                break
    return selected


def _extract_gene_universe(norm_rows: List[dict]) -> List[str]:
    genes = []
    for r in norm_rows:
        symbol = str(r.get("symbol", "")).strip()
        fallback = str(r.get("input_id", "")).strip()
        g = symbol or fallback
        if g:
            genes.append(g.upper())
    dedup = []
    seen = set()
    for g in genes:
        if g not in seen:
            seen.add(g)
            dedup.append(g)
    return dedup


def _extract_term_records(ora_rows: List[dict], gsea_rows: List[dict]) -> List[dict]:
    out: List[dict] = []
    seen = set()
    for r in _top_terms_ora(ora_rows, 80):
        term_name = str(r.get("term_name", "")).strip()
        term_id = str(r.get("term_id", "")).strip()
        key = ("ora", term_id, term_name.lower())
        if not term_name or key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "kind": "ora",
                "term_name": term_name,
                "term_norm": _norm_text(term_name),
                "padj": _safe_float(r.get("padj"), 1.0),
                "pvalue": _safe_float(r.get("pvalue"), 1.0),
                "nes": 0.0,
                "direction": "",
            }
        )
    for r in _top_terms_gsea(gsea_rows, 80):
        term_name = str(r.get("term_name", "")).strip()
        term_id = str(r.get("term_id", "")).strip()
        key = ("gsea", term_id, term_name.lower())
        if not term_name or key in seen:
            continue
        seen.add(key)
        nes = _safe_float(r.get("nes"), 0.0)
        out.append(
            {
                "kind": "gsea",
                "term_name": term_name,
                "term_norm": _norm_text(term_name),
                "padj": _safe_float(r.get("padj"), 1.0),
                "pvalue": _safe_float(r.get("pvalue"), 1.0),
                "nes": nes,
                "direction": "up" if nes >= 0 else "down",
            }
        )
    return out


def _confidence(score: float, marker_hits: int, keyword_hits: int) -> str:
    signal_hits = marker_hits + keyword_hits
    if score >= 0.55 and marker_hits >= 3 and signal_hits >= 6:
        return "high"
    if score >= 0.40 and (marker_hits >= 2 or keyword_hits >= 3):
        return "medium"
    if score >= 0.22:
        return "low"
    return "weak"


def _qualitative_call(marker_hits: int, keyword_hits: int, ora_term_hits: int, gsea_term_hits: int) -> str:
    term_hits = ora_term_hits + gsea_term_hits
    if marker_hits >= 3 and term_hits >= 1:
        return "associated"
    if marker_hits >= 2 and term_hits >= 0:
        return "likely_associated"
    if marker_hits >= 1 and term_hits >= 1:
        return "possible_associated"
    if marker_hits == 0 and keyword_hits >= 2 and term_hits >= 2:
        return "functional_hint_only"
    if marker_hits == 0 and keyword_hits == 0 and term_hits == 0:
        return "not_supported"
    return "weak_or_uncertain"


def _call_label(language: str, call: str) -> str:
    zh = {
        "associated": "相关",
        "likely_associated": "较相关",
        "possible_associated": "可能相关",
        "functional_hint_only": "功能提示相关",
        "weak_or_uncertain": "弱相关/不确定",
        "not_supported": "暂不支持",
    }
    en = {
        "associated": "associated",
        "likely_associated": "likely associated",
        "possible_associated": "possibly associated",
        "functional_hint_only": "functional hint only",
        "weak_or_uncertain": "weak/uncertain",
        "not_supported": "not supported",
    }
    mapping = zh if language == "zh" else en
    return mapping.get(call, call)


def _evidence_basis(marker_ratio: float, keyword_ratio: float, term_hits: int) -> str:
    if marker_ratio >= 0.3 and keyword_ratio >= 0.3:
        return "mixed_marker_functional"
    if marker_ratio >= 0.3:
        return "marker_driven"
    if keyword_ratio >= 0.3 or term_hits >= 4:
        return "functional_driven"
    return "weak_signal"


def _top_term_text(records: List[dict], max_items: int = 4) -> str:
    if not records:
        return ""
    ranked = sorted(
        records,
        key=lambda x: (
            _safe_float(x.get("padj"), 1.0),
            -abs(_safe_float(x.get("nes"), 0.0)),
        ),
    )
    snippets = []
    for rec in ranked[:max_items]:
        if rec.get("kind") == "gsea":
            snippets.append(
                f"{rec['term_name']} (NES={_safe_float(rec.get('nes'), 0.0):.2f}, FDR={_fmt_p(rec.get('padj'))}, {rec.get('direction', '')})"
            )
        else:
            snippets.append(f"{rec['term_name']} (FDR={_fmt_p(rec.get('padj'))})")
    return "; ".join(snippets)


def _interpretation_text(
    language: str,
    cell_type_name: str,
    basis: str,
    marker_hits: int,
    marker_total: int,
    keyword_hits: int,
    keyword_total: int,
    ora_hits: int,
    gsea_hits: int,
    qualitative_call: str,
) -> str:
    call_label = _call_label(language, qualitative_call)
    if language == "zh":
        if basis == "marker_driven":
            return (
                f"定性判断：{cell_type_name} `{call_label}`。推断以标记基因为主（{marker_hits}/{marker_total}），"
                f"功能术语为辅（关键词 {keyword_hits}/{keyword_total}，ORA {ora_hits}，GSEA {gsea_hits}）。"
            )
        if basis == "mixed_marker_functional":
            return (
                f"定性判断：{cell_type_name} `{call_label}`。同时具有标记基因和功能富集双重证据"
                f"（marker {marker_hits}/{marker_total}，关键词 {keyword_hits}/{keyword_total}，ORA {ora_hits}，GSEA {gsea_hits}）。"
            )
        if basis == "functional_driven":
            return (
                f"定性判断：{cell_type_name} `{call_label}`。主要由功能富集驱动（关键词 {keyword_hits}/{keyword_total}，ORA {ora_hits}，GSEA {gsea_hits}），"
                f"当前输入基因中经典 marker 覆盖较低（{marker_hits}/{marker_total}）。"
            )
        return (
            f"定性判断：{cell_type_name} `{call_label}`。当前证据较弱，marker 与功能术语都不充分"
            f"（marker {marker_hits}/{marker_total}，关键词 {keyword_hits}/{keyword_total}）。"
        )

    if basis == "marker_driven":
        return (
            f"Qualitative call: {cell_type_name} is {call_label}. Inference is marker-driven (markers {marker_hits}/{marker_total}) "
            f"with additional functional support (keywords {keyword_hits}/{keyword_total}, ORA {ora_hits}, GSEA {gsea_hits})."
        )
    if basis == "mixed_marker_functional":
        return (
            f"Qualitative call: {cell_type_name} is {call_label}. Supported by both marker and functional evidence "
            f"(markers {marker_hits}/{marker_total}, keywords {keyword_hits}/{keyword_total}, ORA {ora_hits}, GSEA {gsea_hits})."
        )
    if basis == "functional_driven":
        return (
            f"Qualitative call: {cell_type_name} is {call_label}. Mainly function-driven (keywords {keyword_hits}/{keyword_total}, ORA {ora_hits}, GSEA {gsea_hits}) "
            f"while canonical marker coverage is limited (markers {marker_hits}/{marker_total})."
        )
    return (
        f"Qualitative call: {cell_type_name} is {call_label}. Signal is weak with limited marker and pathway evidence "
        f"(markers {marker_hits}/{marker_total}, keywords {keyword_hits}/{keyword_total})."
    )


def _validation_tip(language: str, missing_markers: List[str], confidence: str) -> str:
    marker_text = ",".join(missing_markers[:5]) if missing_markers else ""
    if language == "zh":
        if confidence in {"weak", "low"}:
            return (
                f"建议在原始表达矩阵中优先验证缺失关键 marker：{marker_text or 'N/A'}，"
                "并结合 cluster-level DEGs 复核。"
            )
        return (
            f"建议进一步核查 marker 表达一致性（可优先看 {marker_text or 'N/A'}）并结合已知文献注释。"
        )
    if confidence in {"weak", "low"}:
        return (
            f"Validate in expression matrix first (priority markers: {marker_text or 'N/A'}) "
            "and re-check with cluster-level DEGs."
        )
    return (
        f"Verify marker coherence (priority markers: {marker_text or 'N/A'}) and compare with known references."
    )


def _infer(
    genes: List[str],
    term_records: List[dict],
    selected: Dict[str, dict],
    topk: int,
    language: str,
) -> List[dict]:
    rows: List[dict] = []
    gene_set = set(genes)

    for cell_id, spec in selected.items():
        display_name = str(spec.get("display_name", cell_id))
        markers = sorted({str(x).strip().upper() for x in spec.get("markers", []) if str(x).strip()})
        keywords = sorted({_norm_text(str(x)) for x in spec.get("keywords", []) if _norm_text(str(x))})
        marker_total = len(markers)
        keyword_total = len(keywords)

        marker_hits = sorted(gene_set & set(markers))
        missing_markers = [m for m in markers if m not in gene_set]
        marker_hit_ratio = (len(marker_hits) / marker_total) if marker_total else 0.0

        keyword_hits_set = set()
        matched_ora_terms: List[dict] = []
        matched_gsea_terms: List[dict] = []
        matched_term_names: List[str] = []
        for rec in term_records:
            rec_norm = rec.get("term_norm", "")
            if not rec_norm:
                continue
            hit_keywords = [kw for kw in keywords if kw in rec_norm]
            if not hit_keywords:
                continue
            keyword_hits_set.update(hit_keywords)
            matched_term_names.append(str(rec.get("term_name", "")))
            if rec.get("kind") == "gsea":
                matched_gsea_terms.append(rec)
            else:
                matched_ora_terms.append(rec)

        keyword_hits = sorted(keyword_hits_set)
        keyword_hit_ratio = (len(keyword_hits) / keyword_total) if keyword_total else 0.0
        ora_term_hits = len({x.get("term_name", "") for x in matched_ora_terms if x.get("term_name")})
        gsea_term_hits = len({x.get("term_name", "") for x in matched_gsea_terms if x.get("term_name")})
        gsea_up_hits = sum(1 for x in matched_gsea_terms if x.get("direction") == "up")
        gsea_down_hits = sum(1 for x in matched_gsea_terms if x.get("direction") == "down")

        # Marker evidence still takes highest weight, but functional support is explicitly decomposed.
        ora_support = min(1.0, ora_term_hits / 6.0)
        gsea_support = min(1.0, gsea_term_hits / 6.0)
        combined = 0.50 * marker_hit_ratio + 0.22 * keyword_hit_ratio + 0.18 * ora_support + 0.10 * gsea_support
        basis = _evidence_basis(marker_hit_ratio, keyword_hit_ratio, ora_term_hits + gsea_term_hits)
        qualitative_call = _qualitative_call(
            marker_hits=len(marker_hits),
            keyword_hits=len(keyword_hits),
            ora_term_hits=ora_term_hits,
            gsea_term_hits=gsea_term_hits,
        )
        confidence = _confidence(combined, len(marker_hits), len(keyword_hits))
        interpretation = _interpretation_text(
            language=language,
            cell_type_name=display_name,
            basis=basis,
            marker_hits=len(marker_hits),
            marker_total=marker_total,
            keyword_hits=len(keyword_hits),
            keyword_total=keyword_total,
            ora_hits=ora_term_hits,
            gsea_hits=gsea_term_hits,
            qualitative_call=qualitative_call,
        )
        validation_tip = _validation_tip(language, missing_markers, confidence)

        rows.append(
            {
                "rank": 0,
                "cell_type_id": cell_id,
                "cell_type_name": display_name,
                "qualitative_call": qualitative_call,
                "combined_score": round(combined, 4),
                "evidence_basis": basis,
                "marker_hits": len(marker_hits),
                "marker_total": marker_total,
                "marker_hit_ratio": round(marker_hit_ratio, 4),
                "keyword_hits": len(keyword_hits),
                "keyword_total": keyword_total,
                "keyword_hit_ratio": round(keyword_hit_ratio, 4),
                "ora_term_hits": ora_term_hits,
                "gsea_term_hits": gsea_term_hits,
                "gsea_up_hits": gsea_up_hits,
                "gsea_down_hits": gsea_down_hits,
                "confidence": confidence,
                "evidence_genes": ";".join(marker_hits[:12]),
                "missing_markers": ";".join(missing_markers[:12]),
                "evidence_keywords": ";".join(keyword_hits[:12]),
                "top_ora_terms": _top_term_text(matched_ora_terms, max_items=4),
                "top_gsea_terms": _top_term_text(matched_gsea_terms, max_items=4),
                "evidence_terms": ";".join(matched_term_names[:10]),
                "interpretation": interpretation,
                "validation_suggestion": validation_tip,
            }
        )

    rows.sort(
        key=lambda r: (
            float(r["combined_score"]),
            int(r["marker_hits"]),
            int(r["keyword_hits"]),
            int(r["ora_term_hits"]) + int(r["gsea_term_hits"]),
        ),
        reverse=True,
    )
    for idx, row in enumerate(rows[:topk], start=1):
        row["rank"] = idx
    return rows[:topk]


def _write_interpretation_md(
    outdir: Path,
    rows: List[dict],
    candidates_raw: str,
    species: str,
    genes_n: int,
    term_n: int,
    language: str,
) -> None:
    md_path = outdir / "celltype_interpretation.md"
    if language == "zh":
        lines = [
            "# 细胞类型关联深度解读",
            "",
            f"- 物种: `{species}`",
            f"- 输入基因数: `{genes_n}`",
            f"- 参与解释的富集术语数: `{term_n}`",
            f"- 候选过滤: `{candidates_raw or 'all'}`",
            "",
        ]
    else:
        lines = [
            "# Cell Type Association Deep-Dive",
            "",
            f"- Species: `{species}`",
            f"- Input genes used: `{genes_n}`",
            f"- Enrichment terms used for inference: `{term_n}`",
            f"- Candidate filter: `{candidates_raw or 'all'}`",
            "",
        ]

    if not rows:
        if language == "zh":
            lines.extend(["当前基因与术语证据不足，未推断出稳定的细胞类型关联。", ""])
        else:
            lines.extend(["No stable cell-type signal was inferred from current genes/terms.", ""])
        md_path.write_text("\n".join(lines), encoding="utf-8")
        return

    top_gap = None
    if len(rows) >= 2:
        top_gap = float(rows[0].get("combined_score") or 0.0) - float(rows[1].get("combined_score") or 0.0)

    if language == "zh":
        lines.append("## 结论概览")
        top_call = _call_label(language, str(rows[0].get("qualitative_call", "")))
        lines.append(f"- Top1: `{rows[0]['cell_type_name']}`，定性判断 `{top_call}`，置信度 `{rows[0]['confidence']}`。")
        if top_gap is not None:
            lines.append(f"- Top1 与 Top2 分差: `{top_gap:.3f}`（分差越小，类型越可能混淆）。")
        lines.append("")
    else:
        lines.append("## Executive Summary")
        lines.append(
            f"- Top1: `{rows[0]['cell_type_name']}` with qualitative call `{_call_label(language, str(rows[0].get('qualitative_call', '')))}' and confidence `{rows[0]['confidence']}`."
        )
        if top_gap is not None:
            lines.append(f"- Top1 vs Top2 margin: `{top_gap:.3f}` (smaller margin implies stronger ambiguity).")
        lines.append("")

    lines.extend(
        [
            "| Rank | Cell Type | Qualitative Call | Score | Basis | Marker Hits | Keyword Hits | ORA Terms | GSEA Terms | Confidence |",
            "|---:|---|---|---:|---|---:|---:|---:|---:|---|",
        ]
    )
    for r in rows:
        call_label = _call_label(language, str(r.get("qualitative_call", "")))
        lines.append(
            f"| {r['rank']} | {r['cell_type_name']} | {call_label} | {float(r['combined_score']):.3f} | {r['evidence_basis']} | "
            f"{r['marker_hits']}/{r['marker_total']} | {r['keyword_hits']}/{r['keyword_total']} | "
            f"{r['ora_term_hits']} | {r['gsea_term_hits']} | {r['confidence']} |"
        )
    lines.append("")

    if language == "zh":
        lines.append("## 逐项证据解读")
    else:
        lines.append("## Candidate-by-Candidate Evidence")

    for r in rows:
        lines.append(f"### {r['rank']}. {r['cell_type_name']}")
        lines.append(f"- qualitative call: `{_call_label(language, str(r.get('qualitative_call', '')) )}`")
        lines.append(
            f"- score/confidence: `{r['combined_score']}` / `{r['confidence']}`; basis: `{r['evidence_basis']}`"
        )
        lines.append(
            f"- marker evidence: hit `{r['marker_hits']}/{r['marker_total']}` -> `{r.get('evidence_genes') or 'N/A'}`"
        )
        lines.append(f"- missing key markers: `{r.get('missing_markers') or 'N/A'}`")
        lines.append(
            f"- functional keywords: `{r['keyword_hits']}/{r['keyword_total']}` -> `{r.get('evidence_keywords') or 'N/A'}`"
        )
        lines.append(
            f"- ORA support ({r['ora_term_hits']}): `{r.get('top_ora_terms') or 'N/A'}`"
        )
        lines.append(
            f"- GSEA support ({r['gsea_term_hits']}, up/down={r['gsea_up_hits']}/{r['gsea_down_hits']}): `{r.get('top_gsea_terms') or 'N/A'}`"
        )
        lines.append(f"- interpretation: {r.get('interpretation') or ''}")
        lines.append(f"- suggested validation: {r.get('validation_suggestion') or ''}")
        lines.append("")

    if language == "zh":
        lines.append("## 使用建议")
        lines.append("- 若 Top1/Top2 分差 < 0.08，建议视为并列候选并结合原始表达矩阵复核。")
        lines.append("- 若 marker 命中低而功能术语高，当前结论应定性为“功能关联提示”而非最终注释。")
    else:
        lines.append("## Practical Notes")
        lines.append("- If Top1/Top2 margin < 0.08, treat them as competing candidates.")
        lines.append("- If marker hits are low but functional terms are strong, report as functional hint rather than final annotation.")
    lines.append("")
    md_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir).expanduser().resolve()
    script_dir = Path(__file__).resolve().parent
    language = resolve_language(args.language, args.context_text)

    norm_rows = read_csv(outdir / "normalized_input.tsv")
    ora_rows = _pick_primary(
        read_csv(outdir / "ora_results_python.csv"),
        read_csv(outdir / "ora_results_r.csv"),
    )
    gsea_rows = _pick_primary(
        read_csv(outdir / "gsea_results_python.csv"),
        read_csv(outdir / "gsea_results_r.csv"),
    )

    genes = _extract_gene_universe(norm_rows)
    term_records = _extract_term_records(ora_rows, gsea_rows)
    refs = _load_reference(script_dir, args.species)
    candidates = _parse_candidates(args.candidates)
    selected = _select_celltypes(refs, candidates)

    if candidates and not selected:
        write_csv(outdir / "celltype_association.csv", CELLTYPE_FIELDS, [])
        nohit_text = (
            f"未在参考库中识别到候选细胞类型: {args.candidates}"
            if language == "zh"
            else f"No candidate cell type found in reference: {args.candidates}"
        )
        _write_interpretation_md(
            outdir=outdir,
            rows=[],
            candidates_raw=args.candidates,
            species=args.species,
            genes_n=len(genes),
            term_n=len(term_records),
            language=language,
        )
        payload = {
            "species": args.species,
            "language": language,
            "candidate_filter": args.candidates or "all",
            "gene_count": len(genes),
            "term_count": len(term_records),
            "result_count": 0,
            "top_hit": {},
            "note": nohit_text,
            "available_celltypes": sorted(refs.keys()),
        }
        write_json(outdir / "celltype_association.json", payload)
        return 0

    rows = _infer(genes, term_records, selected, args.topk, language)
    write_csv(outdir / "celltype_association.csv", CELLTYPE_FIELDS, rows)
    _write_interpretation_md(
        outdir=outdir,
        rows=rows,
        candidates_raw=args.candidates,
        species=args.species,
        genes_n=len(genes),
        term_n=len(term_records),
        language=language,
    )

    payload = {
        "species": args.species,
        "language": language,
        "candidate_filter": args.candidates or "all",
        "gene_count": len(genes),
        "term_count": len(term_records),
        "result_count": len(rows),
        "score_margin_top1_top2": (
            (float(rows[0].get("combined_score") or 0.0) - float(rows[1].get("combined_score") or 0.0))
            if len(rows) >= 2
            else None
        ),
        "top_hit": rows[0] if rows else {},
        "results": rows,
    }
    write_json(outdir / "celltype_association.json", payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
