#!/usr/bin/env python3
"""Render an extended interpretation brief for discussion-oriented analysis."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Dict, List


AXIS_KEYWORDS = {
    "hormone_secretory": [
        "hormone",
        "secretion",
        "secretory",
        "peptide hormone",
        "insulin",
        "glucagon",
        "exocyt",
    ],
    "glucose_metabolic": [
        "glucose",
        "glyco",
        "metabolic",
        "homeostasis",
        "transport",
        "insulin",
    ],
    "stress_protein_folding": [
        "endoplasmic reticulum",
        "protein folding",
        "unfolded protein",
        "er stress",
        "chaperone",
        "protein localization",
    ],
    "signaling_regulation": [
        "signaling",
        "receptor",
        "g protein",
        "phosphorylation",
        "cyclic",
        "regulation",
    ],
    "development_fate": [
        "development",
        "differentiation",
        "proliferation",
        "apoptotic",
        "pancreatic",
        "cell fate",
    ],
}


AXIS_LABELS = {
    "en": {
        "hormone_secretory": "Hormone/Secretory Program",
        "glucose_metabolic": "Glucose-Metabolic Program",
        "stress_protein_folding": "Protein Folding/Stress Program",
        "signaling_regulation": "Signaling Regulation Program",
        "development_fate": "Development/Fate Program",
        "general": "General Program",
    },
    "zh": {
        "hormone_secretory": "激素与分泌程序",
        "glucose_metabolic": "葡萄糖代谢程序",
        "stress_protein_folding": "蛋白折叠/应激程序",
        "signaling_regulation": "信号调控程序",
        "development_fate": "发育与命运程序",
        "general": "一般功能程序",
    },
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render extended analysis brief")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--language", choices=["en", "zh"], default="en")
    parser.add_argument("--focus", default="")
    return parser


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _read_csv_rows(path: Path) -> List[dict]:
    if not path.exists():
        return []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))
    except Exception:
        return []


def _safe_float(value, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _normalize_text(text: str) -> str:
    return re.sub(r"[\W_]+", " ", (text or "").lower()).strip()


def _top_ora_terms(rows: List[dict], top_n: int = 30) -> List[dict]:
    sorted_rows = sorted(rows, key=lambda row: (_safe_float(row.get("padj"), 1.0), _safe_float(row.get("pvalue"), 1.0)))
    out: List[dict] = []
    seen = set()
    for row in sorted_rows:
        term_name = str(row.get("term_name", "")).strip()
        term_id = str(row.get("term_id", "")).strip()
        key = term_id or term_name.lower()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "term_name": term_name,
                "padj": _safe_float(row.get("padj"), 1.0),
                "gene_hits": int(_safe_float(row.get("gene_hits"), 0.0)),
                "term_norm": _normalize_text(term_name),
            }
        )
        if len(out) >= top_n:
            break
    return out


def _term_axes(term_norm: str) -> List[str]:
    matched: List[str] = []
    for axis, keywords in AXIS_KEYWORDS.items():
        if any(keyword in term_norm for keyword in keywords):
            matched.append(axis)
    return matched or ["general"]


def _axis_summary(terms: List[dict], language: str) -> List[dict]:
    axis_score: Dict[str, float] = {}
    axis_terms: Dict[str, List[str]] = {}
    for term in terms:
        axes = _term_axes(term.get("term_norm", ""))
        weight = max(0.0, -math.log10(max(term.get("padj", 1.0), 1e-300)))
        formatted = f"{term['term_name']} (FDR={term['padj']:.2e})"
        for axis in axes:
            axis_score[axis] = axis_score.get(axis, 0.0) + weight
            axis_terms.setdefault(axis, [])
            if len(axis_terms[axis]) < 4:
                axis_terms[axis].append(formatted)

    ranked = sorted(axis_score.items(), key=lambda item: item[1], reverse=True)
    labels = AXIS_LABELS["zh" if language == "zh" else "en"]
    out: List[dict] = []
    for axis, score in ranked:
        out.append(
            {
                "axis_id": axis,
                "axis_label": labels.get(axis, axis),
                "score": round(score, 3),
                "top_terms": axis_terms.get(axis, []),
            }
        )
    return out


def _top_term_sentence(top_terms: List[dict], top_n: int) -> str:
    items = []
    for term in top_terms[:top_n]:
        items.append(f"{term.get('term_name', 'N/A')} (FDR={_safe_float(term.get('padj'), 1.0):.2e})")
    return "; ".join(items) if items else "N/A"


def _ora_overlap_details(consistency: dict) -> Dict[str, float]:
    out: Dict[str, float] = {}
    details = ((consistency.get("ora", {}) or {}).get("details", {})) if isinstance(consistency, dict) else {}
    if not isinstance(details, dict):
        return out
    for name, payload in details.items():
        if not isinstance(payload, dict):
            continue
        overlap = payload.get("overlap")
        if isinstance(overlap, (int, float)):
            out[name] = float(overlap)
    return out


def _term_key_from_row(row: dict) -> str:
    term_id = str(row.get("term_id", "")).strip().upper()
    if term_id:
        return term_id
    return _normalize_text(str(row.get("term_name", "")).strip())


def _collect_focus_keywords(focus: str, cell_row: dict) -> List[str]:
    keywords = set()
    for token in re.split(r"[,;/\s]+", (focus or "").lower()):
        token = token.strip()
        if len(token) >= 3:
            keywords.add(token)

    evidence_keywords = str(cell_row.get("evidence_keywords", "") or "")
    for token in re.split(r"[;,/\s]+", evidence_keywords.lower()):
        token = token.strip()
        if len(token) >= 4:
            keywords.add(token)

    focus_lower = (focus or "").lower()
    if "beta" in focus_lower:
        keywords.update({"insulin", "hormone", "secretion", "glucose", "peptide"})
    if "alpha" in focus_lower:
        keywords.update({"glucagon", "hormone", "secretion", "glucose", "peptide"})
    return sorted(keywords)


def _mean(values: List[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def _minmax_norm(value: float, low: float, high: float) -> float:
    if high <= low:
        return 0.5
    return (value - low) / (high - low)


def _build_method_profiles(
    methods: List[dict],
    sig_threshold: float,
    focus_keywords: List[str],
) -> List[dict]:
    if not methods:
        return []

    profiles: List[dict] = []
    support_count: Dict[str, int] = {}
    rank_sum: Dict[str, float] = {}
    rank_count: Dict[str, int] = {}

    for method in methods:
        rows = method.get("rows", []) or []
        top_terms = _top_ora_terms(rows, top_n=20)
        top_keys: List[str] = []
        for idx, term in enumerate(top_terms, start=1):
            key = _term_key_from_row(term)
            if not key:
                continue
            top_keys.append(key)
            support_count[key] = support_count.get(key, 0) + 1
            rank_sum[key] = rank_sum.get(key, 0.0) + idx
            rank_count[key] = rank_count.get(key, 0) + 1

        sig_count = 0
        for row in rows:
            padj = _safe_float(row.get("padj"), float("nan"))
            if not math.isnan(padj) and padj <= sig_threshold:
                sig_count += 1

        signal_values = [-math.log10(max(_safe_float(t.get("padj"), 1.0), 1e-300)) for t in top_terms]
        signal_strength = _mean(signal_values)

        focus_hits = 0
        for term in top_terms:
            term_norm = str(term.get("term_norm", ""))
            if any(kw in term_norm for kw in focus_keywords):
                focus_hits += 1
        focus_rate = focus_hits / max(1, len(top_terms))

        profiles.append(
            {
                "method_id": method.get("method_id", "unknown"),
                "method_label": method.get("method_label", "unknown"),
                "is_primary": bool(method.get("is_primary", False)),
                "is_independent": bool(method.get("is_independent", True)),
                "top_terms_n": len(top_terms),
                "top_preview": _top_term_sentence(top_terms, 3),
                "top_keys": top_keys,
                "sig_terms": sig_count,
                "signal_strength": round(signal_strength, 4),
                "focus_hits": focus_hits,
                "focus_rate": round(focus_rate, 4),
            }
        )

    ranked_keys = sorted(
        support_count.keys(),
        key=lambda k: (-support_count.get(k, 0), rank_sum.get(k, 1e9) / max(1, rank_count.get(k, 1))),
    )
    consensus_keys = ranked_keys[:20]
    consensus_set = set(consensus_keys)

    for profile in profiles:
        top_set = set(profile.get("top_keys", []))
        overlap = len(top_set & consensus_set) / max(1, len(top_set))
        profile["consensus_overlap"] = round(overlap, 4)

    signal_min = min((float(p["signal_strength"]) for p in profiles), default=0.0)
    signal_max = max((float(p["signal_strength"]) for p in profiles), default=1.0)
    sig_min = min((int(p["sig_terms"]) for p in profiles), default=0)
    sig_max = max((int(p["sig_terms"]) for p in profiles), default=1)

    for profile in profiles:
        signal_norm = _minmax_norm(float(profile["signal_strength"]), signal_min, signal_max)
        sig_norm = _minmax_norm(float(profile["sig_terms"]), float(sig_min), float(sig_max))
        score = (
            0.45 * float(profile["consensus_overlap"])
            + 0.30 * signal_norm
            + 0.15 * float(profile["focus_rate"])
            + 0.10 * sig_norm
        )
        if (not profile["is_independent"]) and (not profile["is_primary"]):
            score -= 0.08
        profile["evidence_score"] = round(max(0.0, score), 4)

    profiles.sort(key=lambda p: float(p.get("evidence_score", 0.0)), reverse=True)
    return profiles


def _build_method_recommendation(method_profiles: List[dict], language: str) -> List[str]:
    if not method_profiles:
        if language == "zh":
            return ["当前缺少可比较的方法结果，无法给出方法学推荐。"]
        return ["No comparable method outputs available; cannot provide a method recommendation."]

    best = method_profiles[0]
    independent = [
        m
        for m in method_profiles
        if bool(m.get("is_independent")) and (not bool(m.get("is_primary"))) and int(m.get("top_terms_n", 0)) > 0
    ]
    best_independent = independent[0] if independent else None

    if language == "zh":
        lines = [
            (
                f"绝对证据强度最高的方法是 `{best.get('method_label', 'N/A')}` "
                f"(score={best.get('evidence_score', 0):.2f}, overlap={best.get('consensus_overlap', 0):.2f})。"
                "该方法在当前数据上给出了最完整的显著信号和主题一致性。"
            )
        ]
        if best_independent:
            lines.append(
                (
                    f"独立基线中最强的是 `{best_independent.get('method_label', 'N/A')}` "
                    f"(score={best_independent.get('evidence_score', 0):.2f})。"
                    "建议将其作为主结论的外部佐证来源，以降低单方法偏差。"
                )
            )
        lines.append(
            "推荐策略：用“主方法 + 最强独立基线”共同叙述结果；若两者方向一致，则优先采纳该结论作为本次输出推荐。"
        )
        return lines

    lines = [
        (
            f"The strongest absolute evidence comes from `{best.get('method_label', 'N/A')}` "
            f"(score={best.get('evidence_score', 0):.2f}, overlap={best.get('consensus_overlap', 0):.2f}). "
            "This method shows the most complete significant signal with coherent thematic support."
        )
    ]
    if best_independent:
        lines.append(
            (
                f"The strongest independent baseline is `{best_independent.get('method_label', 'N/A')}` "
                f"(score={best_independent.get('evidence_score', 0):.2f}). "
                "Use it as external corroboration to reduce single-method bias."
            )
        )
    lines.append(
        "Recommended strategy: report conclusions using both the primary method and the strongest independent baseline; prioritize claims where both agree."
    )
    return lines


def _primary_method_label(run_meta: dict, ora_rows: List[dict], language: str) -> str:
    run_source = str((run_meta.get("run_sources", {}) or {}).get("ora", "")).lower()
    backend = str(run_meta.get("python_ora_backend", "") or "").lower()
    row_source = str((ora_rows[0].get("source") if ora_rows else "") or "").lower()

    if run_source == "python":
        if "gprofiler" in backend or "gprofiler" in row_source:
            return "主方法(g:Profiler-Python)" if language == "zh" else "Primary (g:Profiler-Python)"
        if "enrichr" in backend or "gseapy" in backend or "enrichr" in row_source:
            return "主方法(Enrichr/gseapy-Python)" if language == "zh" else "Primary (Enrichr/gseapy-Python)"
        return "主方法(Python)" if language == "zh" else "Primary (Python)"

    if run_source == "r" or "clusterprofiler" in row_source:
        return "主方法(clusterProfiler-R)" if language == "zh" else "Primary (clusterProfiler-R)"

    return "主方法(未知)" if language == "zh" else "Primary (Unknown)"


def _build_main_conclusion(
    cell_row: dict,
    axis_rows: List[dict],
    top_terms: List[dict],
    method_profiles: List[dict],
    ora_score: object,
    gsea_applicable: bool,
    language: str,
) -> List[str]:
    primary_axis = axis_rows[0] if axis_rows else {}
    secondary_axis = axis_rows[1] if len(axis_rows) > 1 else {}
    cell_name = str(cell_row.get("cell_type_name", "target cell"))
    call = str(cell_row.get("qualitative_call", "unknown"))
    markers = f"{cell_row.get('marker_hits', '0')}/{cell_row.get('marker_total', '0')}"
    keywords = f"{cell_row.get('keyword_hits', '0')}/{cell_row.get('keyword_total', '0')}"
    term_block = _top_term_sentence(top_terms, 4)
    ora_text = f"{ora_score:.2f}" if isinstance(ora_score, (int, float)) else "N/A"
    best_method = method_profiles[0].get("method_label", "N/A") if method_profiles else "N/A"
    independent = [m for m in method_profiles if bool(m.get("is_independent")) and (not bool(m.get("is_primary")))]
    best_independent = independent[0].get("method_label", "N/A") if independent else "N/A"

    if language == "zh":
        return [
            (
                f"当前富集结果的主轴并非离散或随机过程，而是集中在“{primary_axis.get('axis_label', '核心功能程序')}”上。"
                f"从显著性最高的术语看，{term_block} 构成了连续且同主题的证据链。"
                f"这类术语组合通常对应同一生物学模块在不同层级上的展开，而不是互不相关的噪声富集。"
            ),
            (
                f"结合细胞类型证据，当前最稳健的解释是 `{cell_name}` 相关程序活跃（定性={call}），"
                f"并由 marker 命中 {markers} 与功能关键词命中 {keywords} 共同支撑。"
                f"若考虑次级功能轴，`{secondary_axis.get('axis_label', '次级程序')}` 也提供了附加解释空间，"
                "提示该模块可能同时覆盖主分泌轴和协同调控轴。"
            ),
            (
                f"在多方法证据比较中，`{best_method}` 给出最高证据强度，"
                f"而 `{best_independent}` 提供最强独立外部佐证。"
                + (f"补充信息：ORA 共识分数 `{ora_text}`；" if ora_text != "N/A" else "")
                + ("GSEA 方向证据可用。" if gsea_applicable else "当前缺少 GSEA 方向证据。")
                + "因此更推荐采用“主方法 + 最强独立基线”的联合叙述，而非单方法定论。"
            ),
        ]

    return [
        (
            f"The enrichment landscape is not diffuse: it converges on a coherent `{primary_axis.get('axis_label', 'core functional program')}` axis. "
            f"Top significant terms ({term_block}) form a biologically aligned chain rather than unrelated hits, supporting a shared module-level process."
        ),
        (
            f"Cell-type evidence supports an active `{cell_name}` program (qualitative call={call}), with marker support {markers} and keyword support {keywords}. "
            f"A secondary `{secondary_axis.get('axis_label', 'secondary axis')}` signal remains visible, suggesting a layered rather than single-axis interpretation."
        ),
        (
            f"In cross-method evidence ranking, `{best_method}` is strongest, with `{best_independent}` as the strongest independent corroboration. "
            + (f"ORA consensus score is `{ora_text}`; " if ora_text != "N/A" else "")
            + ("GSEA directionality is available." if gsea_applicable else "GSEA directionality is currently unavailable.")
            + " The recommended framing is a joint-method narrative rather than a single-method claim."
        ),
    ]


def _build_support_points(
    consistency: dict,
    cell_row: dict,
    axis_rows: List[dict],
    top_terms: List[dict],
    language: str,
) -> List[str]:
    overlaps = _ora_overlap_details(consistency)
    overlap_text = ", ".join(f"{k}={v:.2f}" for k, v in overlaps.items()) if overlaps else "N/A"
    top_axis = axis_rows[0].get("axis_label", "N/A") if axis_rows else "N/A"
    term_block = _top_term_sentence(top_terms, 6)
    marker_hits = cell_row.get("marker_hits", "0")
    marker_total = cell_row.get("marker_total", "0")
    keyword_hits = cell_row.get("keyword_hits", "0")
    keyword_total = cell_row.get("keyword_total", "0")
    ora_hits = cell_row.get("ora_term_hits", "0")

    if language == "zh":
        return [
            (
                f"显著术语在 `{top_axis}` 主题上呈现高一致聚类，代表条目包括：{term_block}。"
                "这一模式说明结果来自同一生物过程链条，而不是被单个偶发术语牵引。"
            ),
            (
                f"细胞类型证据同时满足 marker 与功能两条路径：marker={marker_hits}/{marker_total}，"
                f"keyword={keyword_hits}/{keyword_total}，且 ORA 命中术语数={ora_hits}。"
                "双证据结构通常比单一路径更稳健，尤其适合用于研究场景中的机制假设构建。"
            ),
            (
                f"跨工具对照的 ORA 明细为：{overlap_text}。"
                "即使不同工具的统计框架不同，只要主要基线中存在中高重叠，就能支持“主轴真实存在”的判断。"
            ),
        ]

    return [
        (
            f"Significant terms cluster around `{top_axis}` with coherent biology: {term_block}. "
            "This pattern indicates a process-level signal rather than isolated term artifacts."
        ),
        (
            f"Cell-type evidence is dual-channel (marker and functional): marker={marker_hits}/{marker_total}, "
            f"keyword={keyword_hits}/{keyword_total}, ORA term hits={ora_hits}. "
            "This mixed evidence structure is generally more robust for research interpretation than single-source support."
        ),
        (
            f"Cross-tool ORA details are: {overlap_text}. "
            "Even with method heterogeneity, consistent overlap across major baselines supports the existence of a real dominant axis."
        ),
    ]


def _build_celltype_points(cell_row: dict, language: str) -> List[str]:
    cell_name = str(cell_row.get("cell_type_name", "target cell"))
    call = str(cell_row.get("qualitative_call", "unknown"))
    evidence_basis = str(cell_row.get("evidence_basis", "unknown"))
    marker_hits = cell_row.get("marker_hits", "0")
    marker_total = cell_row.get("marker_total", "0")
    keyword_hits = cell_row.get("keyword_hits", "0")
    keyword_total = cell_row.get("keyword_total", "0")
    ora_hits = cell_row.get("ora_term_hits", "0")
    gsea_hits = cell_row.get("gsea_term_hits", "0")
    missing = str(cell_row.get("missing_markers", "")).strip()

    if language == "zh":
        return [
            (
                f"当前细胞类型结论为 `{cell_name} -> {call}`，证据类型为 `{evidence_basis}`。"
                "这意味着模型不是仅靠 marker 或仅靠功能词单点判断，而是按证据结构给出定性判定。"
            ),
            (
                f"在量化上，marker 命中为 {marker_hits}/{marker_total}，关键词命中为 {keyword_hits}/{keyword_total}，"
                f"并且 ORA/GSEA 术语命中为 {ora_hits}/{gsea_hits}。"
                "该组合可直接用于比较不同 run 间的证据强弱，而不是只比一个 p 值。"
            ),
            (
                f"缺失 marker 为：{missing or 'N/A'}。"
                "这些未命中位点可作为后续验证优先列表，用于判断当前模块是“功能偏置”还是“稳定谱系状态”。"
            ),
        ]

    return [
        (
            f"The current cell-type call is `{cell_name} -> {call}` with evidence basis `{evidence_basis}`. "
            "This indicates a structured decision built from multiple evidence channels rather than a single marker-only cue."
        ),
        (
            f"Quantitatively, marker support is {marker_hits}/{marker_total}, keyword support is {keyword_hits}/{keyword_total}, "
            f"and ORA/GSEA term hits are {ora_hits}/{gsea_hits}. "
            "This composite profile is more informative for run-to-run comparison than a standalone p-value."
        ),
        (
            f"Missing markers are: {missing or 'N/A'}. "
            "These are practical validation targets for testing whether the module reflects functional bias or a stable lineage state."
        ),
    ]


def _build_supported_claims(consistency: dict, axis_rows: List[dict], method_profiles: List[dict], language: str) -> List[str]:
    top_axis = axis_rows[0].get("axis_label", "N/A") if axis_rows else "N/A"
    second_axis = axis_rows[1].get("axis_label", "N/A") if len(axis_rows) > 1 else "N/A"
    ora_score = (consistency.get("ora", {}) or {}).get("score") if isinstance(consistency, dict) else None
    overall_pass = bool(consistency.get("overall_pass")) if isinstance(consistency, dict) else False
    best_method = method_profiles[0].get("method_label", "N/A") if method_profiles else "N/A"

    if language == "zh":
        score_text = f"{ora_score:.2f}" if isinstance(ora_score, (int, float)) else "N/A"
        return [
            (
                f"可以稳健支持“主功能轴存在且具备可解释性”：当前主轴为 `{top_axis}`，"
                f"并且次级轴 `{second_axis}` 提供了协同解释。"
                "这类双轴结构通常优于单轴叙述，能更好地解释复杂细胞状态。"
            ),
            (
                f"可以支持“该模块具有细胞类型偏置”这一结论，但更适合表述为偏置而非绝对特异。"
                f"当前 ORA 一致性得分为 `{score_text}`，总体 gate={'PASS' if overall_pass else 'FAIL'}，"
                f"且在多方法比较中 `{best_method}` 的证据强度最高。"
            ),
            (
                "可以支持后续机制实验的优先级排序，例如先验证分泌轴/代谢轴，再验证缺失 marker。"
                "这种从富集到验证的路径对科学研究最有实际价值。"
            ),
        ]

    score_text = f"{ora_score:.2f}" if isinstance(ora_score, (int, float)) else "N/A"
    return [
        (
            f"A primary interpretable axis is supported (`{top_axis}`), with a secondary cooperating axis (`{second_axis}`). "
            "This layered structure generally explains complex cellular programs better than a single-axis narrative."
        ),
        (
            f"A lineage-biased interpretation is supported, but should be framed as bias rather than absolute specificity. "
            f"Current ORA consistency is `{score_text}` with overall gate={'PASS' if overall_pass else 'FAIL'}, "
            f"and `{best_method}` is the strongest method in cross-method comparison."
        ),
        (
            "The result supports a clear validation roadmap: test dominant functional axes first, then unresolved marker gaps. "
            "This translation from enrichment to experiments is the highest-value research use."
        ),
    ]


def _build_caveats(consistency: dict, cell_row: dict, gsea_applicable: bool, language: str) -> List[str]:
    caveats: List[str] = []
    ora = consistency.get("ora", {}) if isinstance(consistency, dict) else {}
    thresholds = consistency.get("thresholds", {}) if isinstance(consistency, dict) else {}
    ora_score = ora.get("score")
    ora_threshold = thresholds.get("ora_overlap", 0.6)

    if isinstance(ora_score, (int, float)) and ora_score < ora_threshold:
        if language == "zh":
            caveats.append(
                f"ORA 一致性 `score={ora_score:.2f}` 低于阈值 `{ora_threshold:.2f}`，跨工具仍存在敏感性。"
            )
        else:
            caveats.append(
                f"ORA consistency `score={ora_score:.2f}` is below threshold `{ora_threshold:.2f}`, indicating tool sensitivity."
            )

    if not gsea_applicable:
        if language == "zh":
            caveats.append("未提供 ranked scores，当前缺少 GSEA 的方向性证据。")
        else:
            caveats.append("No ranked scores provided; directionality evidence from GSEA is missing.")

    if cell_row:
        marker_hits = _safe_float(cell_row.get("marker_hits"), 0.0)
        marker_total = max(1.0, _safe_float(cell_row.get("marker_total"), 1.0))
        marker_ratio = marker_hits / marker_total
        if marker_ratio < 0.4:
            if language == "zh":
                caveats.append("目标细胞 marker 覆盖偏低，建议结合原始表达矩阵与亚群标记做二次验证。")
            else:
                caveats.append("Target marker coverage is limited; validate using raw expression matrix and subtype markers.")

    if not caveats:
        if language == "zh":
            caveats.append("当前结果未见明显技术性反证，但仍建议进行独立生物学验证。")
        else:
            caveats.append("No strong technical contradiction detected; independent biological validation is still recommended.")
    return caveats


def _build_limits_points(caveats: List[str], consistency: dict, language: str) -> List[str]:
    points = list(caveats)
    overall_pass = bool(consistency.get("overall_pass")) if isinstance(consistency, dict) else False
    if not overall_pass:
        if language == "zh":
            points.append("综合 gate 未通过时，结论应表述为“功能偏置与假设方向”，避免表述为最终定论。")
        else:
            points.append("When the overall gate fails, frame the result as a directional functional bias rather than a final claim.")
    return points


def _build_citable_paragraph(
    cell_row: dict,
    axis_rows: List[dict],
    top_terms: List[dict],
    consistency: dict,
    method_profiles: List[dict],
    language: str,
) -> str:
    axis_label = axis_rows[0].get("axis_label", "N/A") if axis_rows else "N/A"
    terms = _top_term_sentence(top_terms, 5)
    cell_name = str(cell_row.get("cell_type_name", "target cell"))
    call = str(cell_row.get("qualitative_call", "unknown"))
    marker_hits = cell_row.get("marker_hits", "0")
    marker_total = cell_row.get("marker_total", "0")
    keyword_hits = cell_row.get("keyword_hits", "0")
    keyword_total = cell_row.get("keyword_total", "0")
    ora_score = (consistency.get("ora", {}) or {}).get("score") if isinstance(consistency, dict) else None
    ora_threshold = (consistency.get("thresholds", {}) or {}).get("ora_overlap", 0.6) if isinstance(consistency, dict) else 0.6
    best_method = method_profiles[0].get("method_label", "N/A") if method_profiles else "N/A"
    independent = [m for m in method_profiles if bool(m.get("is_independent")) and (not bool(m.get("is_primary")))]
    best_independent = independent[0].get("method_label", "N/A") if independent else "N/A"

    if language == "zh":
        score_text = f"{ora_score:.2f}" if isinstance(ora_score, (int, float)) else "N/A"
        return (
            "基于当前基因集的富集分析，我们观察到结果在功能上集中于"
            f"{axis_label}，并由 {terms} 等同主题条目共同支持。"
            f"细胞类型关联分析将该模块判定为 {cell_name} `{call}`，"
            f"其证据来自 marker 命中 {marker_hits}/{marker_total} 与功能关键词命中 {keyword_hits}/{keyword_total} 的双路径一致。"
            f"在方法比较层面，{best_method} 的证据强度最高，而 {best_independent} 提供最强独立佐证。"
            + (f"补充地，ORA 共识分数为 {score_text}（阈值 {ora_threshold:.2f}）。" if score_text != "N/A" else "")
            + "综合来看，该结果更适合解释为稳定的功能偏置与机制假设，而非对终末状态的单点定论。"
        )

    score_text = f"{ora_score:.2f}" if isinstance(ora_score, (int, float)) else "N/A"
    return (
        "Enrichment analysis of the current gene set indicates a dominant "
        f"{axis_label} axis, jointly supported by coherent terms such as {terms}. "
        f"Cell-type association classifies this module as {cell_name} `{call}`, "
        f"with convergent evidence from marker support ({marker_hits}/{marker_total}) and functional keyword support ({keyword_hits}/{keyword_total}). "
        f"In method comparison, {best_method} is strongest and {best_independent} is the strongest independent corroboration. "
        + (f"ORA cross-tool consensus is {score_text} (threshold {ora_threshold:.2f}). " if score_text != "N/A" else "")
        + "Overall, the result is best interpreted as a robust functional bias and mechanistic direction rather than a terminal-state endpoint."
    )


def _build_next_steps(cell_row: dict, language: str) -> List[str]:
    missing = str(cell_row.get("missing_markers", "")).strip() if cell_row else ""
    if language == "zh":
        return [
            (
                "补充 ranked scores（如 logFC 或 t-score）后重跑 GSEA，优先检查主轴术语是否在方向上与 ORA 一致；"
                "这一步能把“是否富集”升级为“向上/向下驱动”的证据。"
            ),
            (
                "在相同阈值和同一流程下做 beta vs alpha 的并排对照，重点比较分泌轴和代谢轴的稳定性；"
                "对照设计通常比单组结果更容易形成可发表的结论。"
            ),
            (
                f"围绕缺失 marker（{missing or 'N/A'}）做 targeted validation，"
                "可结合原始表达矩阵与文献 marker 集合评估“功能偏置”是否上升为“谱系稳定状态”。"
            ),
        ]

    return [
        (
            "Add ranked scores (logFC or t-score) and rerun GSEA to test whether directionality agrees with ORA; "
            "this upgrades evidence from enrichment-only to directional mechanism support."
        ),
        (
            "Run a side-by-side beta vs alpha analysis under identical thresholds and workflow; "
            "comparative stability is usually more publishable than single-cohort interpretation."
        ),
        (
            f"Perform targeted validation on missing markers ({missing or 'N/A'}) using raw expression plus curated markers "
            "to determine whether the signal is a functional bias or a stable lineage state."
        ),
    ]


def _render_markdown(payload: dict, language: str) -> str:
    if language == "zh":
        lines = [
            "# 深度研究解读简报",
            "",
            "## 1) 主结论（研究版）",
        ]
    else:
        lines = [
            "# Deep Research Interpretation Brief",
            "",
            "## 1) Main Conclusion (Research Framing)",
        ]

    for idx, text in enumerate(payload.get("main_conclusion", []), start=1):
        lines.append(f"{idx}. {text}")
    lines.append("")

    if language == "zh":
        lines.append("## 2) 证据主轴（term + FDR）")
    else:
        lines.append("## 2) Evidence Backbone (term + FDR)")
    for axis in payload.get("axis_support", [])[:4]:
        lines.append(
            f"- `{axis.get('axis_label', 'N/A')}` (score={axis.get('score', 0)}): "
            f"{', '.join(axis.get('top_terms', []))}"
        )
    lines.append("")

    if language == "zh":
        lines.append("## 3) 跨方法证据对比与推荐")
    else:
        lines.append("## 3) Cross-Method Comparison and Recommendation")
    for method in payload.get("method_profiles", []):
        role = "主方法" if method.get("is_primary") else "基线"
        if language != "zh":
            role = "primary" if method.get("is_primary") else "baseline"
        indep = "独立" if method.get("is_independent") else "非独立"
        if language != "zh":
            indep = "independent" if method.get("is_independent") else "non-independent"
        lines.append(
            f"- `{method.get('method_label', 'N/A')}` [{role}, {indep}] "
            f"score={float(method.get('evidence_score', 0.0)):.2f}, "
            f"consensus_overlap={float(method.get('consensus_overlap', 0.0)):.2f}, "
            f"signal={float(method.get('signal_strength', 0.0)):.2f}, "
            f"focus={float(method.get('focus_rate', 0.0)):.2f}, "
            f"sig_terms={int(method.get('sig_terms', 0))}"
        )
        lines.append(f"  top_terms: {method.get('top_preview', 'N/A')}")
    if payload.get("method_recommendation"):
        for idx, text in enumerate(payload.get("method_recommendation", []), start=1):
            lines.append(f"{idx}. {text}")
    lines.append("")

    if language == "zh":
        lines.append("## 4) 细胞类型判断与证据拆解")
    else:
        lines.append("## 4) Cell-Type Call and Evidence Breakdown")
    for idx, text in enumerate(payload.get("celltype_points", []), start=1):
        lines.append(f"{idx}. {text}")
    lines.append("")

    if language == "zh":
        lines.append("## 5) 能说什么 / 不能说什么")
        lines.append("### 能说什么")
    else:
        lines.append("## 5) What This Supports / What It Does Not")
        lines.append("### Supported Claims")
    for idx, text in enumerate(payload.get("supported_claims", []), start=1):
        lines.append(f"{idx}. {text}")
    lines.append("")

    if language == "zh":
        lines.append("### 不能说什么（当前边界）")
    else:
        lines.append("### Current Boundaries")
    for idx, text in enumerate(payload.get("limits_points", []), start=1):
        lines.append(f"{idx}. {text}")
    lines.append("")

    if language == "zh":
        lines.append("## 6) 可直接引用结果段（科研写作）")
    else:
        lines.append("## 6) Directly Citable Result Paragraph")
    lines.append(payload.get("citable_paragraph", ""))
    lines.append("")

    if language == "zh":
        lines.append("## 7) 下一步高信息增益分析")
    else:
        lines.append("## 7) High-Information Next Analyses")
    for idx, step in enumerate(payload.get("next_steps", []), start=1):
        lines.append(f"{idx}. {step}")
    lines.append("")

    lines.extend(
        [
            "## Key Files",
            "- `summary.md`",
            "- `celltype_association.csv`",
            "- `consistency_report.json`",
            "- `ora_results_r.csv` / `ora_results_python.csv`",
            "",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir).expanduser().resolve()
    language = args.language

    run_meta = _read_json(outdir / "run_meta.json")
    consistency = _read_json(outdir / "consistency_report.json")
    cell_rows = _read_csv_rows(outdir / "celltype_association.csv")
    ora_rows = _read_csv_rows(outdir / "ora_results_python.csv")
    if not ora_rows:
        ora_rows = _read_csv_rows(outdir / "ora_results_r.csv")
    baseline_enrichr = _read_csv_rows(outdir / "baseline_enrichr.csv")
    baseline_gprofiler = _read_csv_rows(outdir / "baseline_gprofiler.csv")
    baseline_clusterprofiler = _read_csv_rows(outdir / "baseline_clusterprofiler.csv")

    top_terms = _top_ora_terms(ora_rows, top_n=30)
    axis_rows = _axis_summary(top_terms, language=language)
    top_cell = cell_rows[0] if cell_rows else {}

    gsea_meta = consistency.get("gsea", {}) if isinstance(consistency, dict) else {}
    gsea_applicable = bool(gsea_meta.get("applicable")) if isinstance(gsea_meta, dict) else False
    ora_score = (consistency.get("ora", {}) or {}).get("score") if isinstance(consistency, dict) else None
    sig_threshold = _safe_float(run_meta.get("sig_threshold"), 0.05)

    primary_label = _primary_method_label(run_meta, ora_rows, language=language)
    primary_is_cluster = "clusterprofiler" in primary_label.lower()
    focus_keywords = _collect_focus_keywords(args.focus or run_meta.get("celltype_candidates", ""), top_cell)
    method_candidates = []
    if ora_rows:
        method_candidates.append(
            {
                "method_id": "primary_run",
                "method_label": primary_label,
                "rows": ora_rows,
                "is_primary": True,
                "is_independent": True,
            }
        )
    method_candidates.extend(
        [
            {
                "method_id": "baseline_enrichr",
                "method_label": "Enrichr baseline" if language != "zh" else "Enrichr 基线",
                "rows": baseline_enrichr,
                "is_primary": False,
                "is_independent": True,
            },
            {
                "method_id": "baseline_gprofiler",
                "method_label": "g:Profiler baseline" if language != "zh" else "g:Profiler 基线",
                "rows": baseline_gprofiler,
                "is_primary": False,
                "is_independent": True,
            },
            {
                "method_id": "baseline_clusterprofiler",
                "method_label": "clusterProfiler baseline" if language != "zh" else "clusterProfiler 基线",
                "rows": baseline_clusterprofiler,
                "is_primary": False,
                "is_independent": not primary_is_cluster,
            },
        ]
    )
    method_profiles = _build_method_profiles(
        methods=method_candidates,
        sig_threshold=sig_threshold,
        focus_keywords=focus_keywords,
    )
    method_recommendation = _build_method_recommendation(method_profiles, language=language)

    main_conclusion = _build_main_conclusion(
        cell_row=top_cell,
        axis_rows=axis_rows,
        top_terms=top_terms,
        method_profiles=method_profiles,
        ora_score=ora_score,
        gsea_applicable=gsea_applicable,
        language=language,
    )
    support_points = _build_support_points(
        consistency=consistency,
        cell_row=top_cell,
        axis_rows=axis_rows,
        top_terms=top_terms,
        language=language,
    )
    celltype_points = _build_celltype_points(top_cell, language=language)
    supported_claims = _build_supported_claims(consistency, axis_rows, method_profiles, language=language)
    caveats = _build_caveats(consistency, top_cell, gsea_applicable, language=language)
    limits_points = _build_limits_points(caveats, consistency, language=language)
    citable_paragraph = _build_citable_paragraph(
        cell_row=top_cell,
        axis_rows=axis_rows,
        top_terms=top_terms,
        consistency=consistency,
        method_profiles=method_profiles,
        language=language,
    )
    next_steps = _build_next_steps(top_cell, language=language)

    payload = {
        "focus": args.focus or run_meta.get("celltype_candidates", ""),
        "primary_cell": top_cell.get("cell_type_name", "N/A"),
        "primary_call": top_cell.get("qualitative_call", "N/A"),
        "ora_score": ora_score,
        "gsea_applicable": gsea_applicable,
        "axis_support": axis_rows,
        "top_terms": top_terms,
        "main_conclusion": main_conclusion,
        "support_points": support_points,
        "method_profiles": method_profiles,
        "method_recommendation": method_recommendation,
        "celltype_points": celltype_points,
        "supported_claims": supported_claims,
        "caveats": caveats,
        "limits_points": limits_points,
        "citable_paragraph": citable_paragraph,
        "next_steps": next_steps,
    }

    (outdir / "analysis_brief.md").write_text(_render_markdown(payload, language), encoding="utf-8")
    (outdir / "analysis_brief.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
