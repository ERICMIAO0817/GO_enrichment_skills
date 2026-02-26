#!/usr/bin/env python3
"""Parse natural-language GO enrichment requests into normalized parameters."""

from __future__ import annotations

import argparse
import ast
import json
import re
import shlex
import sys
from typing import List, Optional


REQUIRED_FIELDS = ("genes", "sig_threshold")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Parse natural-language GO/GSEA request text")
    p.add_argument("--request", required=True, help="Natural-language request text")
    p.add_argument("--format", choices=["json", "shell"], default="json")
    return p


def _dedupe_keep_order(items: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for item in items:
        key = item.upper()
        if key in seen:
            continue
        seen.add(key)
        out.append(item)
    return out


def _clean_gene_token(token: str) -> str:
    token = token.strip().strip("'\"")
    token = token.replace("，", ",")
    return token


def _parse_genes_from_brackets(text: str) -> List[str]:
    candidates = re.findall(r"\[[^\[\]]+\]", text, flags=re.DOTALL)
    for candidate in candidates:
        try:
            parsed = ast.literal_eval(candidate)
        except Exception:
            parsed = None

        if isinstance(parsed, list):
            genes = [_clean_gene_token(str(item)) for item in parsed if str(item).strip()]
            genes = [g for g in genes if re.match(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", g)]
            if genes:
                return _dedupe_keep_order(genes)

        raw = candidate.strip().strip("[]")
        tokens = re.split(r"[,\s]+", raw)
        genes = [_clean_gene_token(t) for t in tokens if t.strip()]
        genes = [g for g in genes if re.match(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", g)]
        if genes:
            return _dedupe_keep_order(genes)
    return []


def _parse_genes_from_labeled_segment(text: str) -> List[str]:
    patterns = [
        r"(?:genes?|gene\s*list|基因(?:列表)?)[^:：=]{0,20}[:：=]\s*([^\n;；。]+)",
        r"(?:跑这组基因|这组基因)\s*[:：]?\s*([^\n;；。]+)",
    ]
    for pattern in patterns:
        m = re.search(pattern, text, flags=re.IGNORECASE)
        if not m:
            continue
        segment = m.group(1).strip()
        segment = re.split(r"(?:sig[_-]?threshold|species|ontology|focus|并重点看|并关注)", segment, maxsplit=1, flags=re.IGNORECASE)[0]
        tokens = re.split(r"[,\s]+", segment.replace("，", ","))
        genes = [_clean_gene_token(t) for t in tokens if t.strip()]
        genes = [g for g in genes if re.match(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", g)]
        if genes:
            return _dedupe_keep_order(genes)
    return []


def _parse_genes(text: str) -> List[str]:
    genes = _parse_genes_from_brackets(text)
    if genes:
        return genes
    return _parse_genes_from_labeled_segment(text)


def _parse_sig_threshold(text: str) -> Optional[float]:
    patterns = [
        r"sig[_\-\s]?threshold\s*[:=]\s*([0-9]*\.?[0-9]+(?:e-?[0-9]+)?)",
        r"fdr\s*[:=]\s*([0-9]*\.?[0-9]+(?:e-?[0-9]+)?)",
        r"阈值\s*[:：=]?\s*([0-9]*\.?[0-9]+(?:e-?[0-9]+)?)",
    ]
    for pattern in patterns:
        m = re.search(pattern, text, flags=re.IGNORECASE)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                return None
    return None


def _parse_species(text: str) -> str:
    lowered = text.lower()
    if re.search(r"\bmouse\b", lowered) or ("小鼠" in text) or ("鼠" in text and "人" not in text):
        return "mouse"
    if re.search(r"\bhuman\b", lowered) or ("人类" in text) or ("人源" in text):
        return "human"

    m = re.search(r"(?:species|物种)\s*[:：=]\s*([A-Za-z\u4e00-\u9fa5]+)", text, flags=re.IGNORECASE)
    if m:
        value = m.group(1).strip().lower()
        if value in {"mouse", "mm", "小鼠", "鼠"}:
            return "mouse"
        if value in {"human", "hs", "人", "人类"}:
            return "human"
    return "human"


def _parse_ontology(text: str) -> str:
    m = re.search(r"ontology\s*[:=]\s*([A-Za-z,\s]+)", text, flags=re.IGNORECASE)
    if m:
        value = m.group(1).strip().upper().replace(" ", "")
        if value:
            if value == "ALL":
                return "ALL"
            terms = [x for x in value.split(",") if x in {"BP", "MF", "CC"}]
            if terms:
                return ",".join(terms)
    if re.search(r"\bALL\b", text, flags=re.IGNORECASE):
        return "ALL"
    return "ALL"


def _parse_focus(text: str) -> str:
    patterns = [
        r"重点看\s*([^\n，。,;；]+?)(?:的关联|关联|相关|$)",
        r"关注\s*([^\n，。,;；]+?)(?:的关联|关联|相关|$)",
        r"focus(?:\s+on)?\s*[:=]?\s*([^\n,.;]+)",
    ]
    for pattern in patterns:
        m = re.search(pattern, text, flags=re.IGNORECASE)
        if not m:
            continue
        value = m.group(1).strip()
        value = re.sub(r"\s+and\s+", ",", value, flags=re.IGNORECASE)
        value = value.replace("，", ",")
        value = re.sub(r"\s*,\s*", ",", value)
        value = value.strip(", ")
        if value:
            return value
    return ""


def parse_request(text: str) -> dict:
    genes = _parse_genes(text)
    sig_threshold = _parse_sig_threshold(text)
    species = _parse_species(text)
    ontology = _parse_ontology(text)
    focus = _parse_focus(text)

    payload = {
        "genes": genes,
        "genes_inline": ",".join(genes),
        "sig_threshold": sig_threshold,
        "species": species,
        "ontology": ontology,
        "focus": focus,
    }

    missing = []
    if not genes:
        missing.append("genes")
    if sig_threshold is None:
        missing.append("sig_threshold")
    payload["missing_required"] = missing
    payload["required_fields"] = list(REQUIRED_FIELDS)
    return payload


def _as_shell(payload: dict) -> str:
    def q(value: object) -> str:
        if value is None:
            value = ""
        return shlex.quote(str(value))

    lines = [
        f"PARSED_GENES={q(payload.get('genes_inline', ''))}",
        f"PARSED_SIG_THRESHOLD={q(payload.get('sig_threshold', ''))}",
        f"PARSED_SPECIES={q(payload.get('species', 'human'))}",
        f"PARSED_ONTOLOGY={q(payload.get('ontology', 'ALL'))}",
        f"PARSED_FOCUS={q(payload.get('focus', ''))}",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    args = build_parser().parse_args()
    payload = parse_request(args.request)

    missing = payload.get("missing_required", [])
    if args.format == "shell":
        sys.stdout.write(_as_shell(payload))
    else:
        sys.stdout.write(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")

    if missing:
        sys.stderr.write(f"Missing required fields: {', '.join(missing)}\n")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
