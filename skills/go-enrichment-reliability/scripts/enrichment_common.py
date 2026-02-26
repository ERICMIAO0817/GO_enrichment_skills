#!/usr/bin/env python3
"""Common helpers for GO/GSEA enrichment scripts."""

from __future__ import annotations

import csv
import datetime as dt
import json
import os
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

ORA_FIELDS = [
    "term_id",
    "term_name",
    "ontology",
    "source",
    "pvalue",
    "padj",
    "gene_hits",
    "gene_ratio",
    "engine",
]

GSEA_FIELDS = [
    "term_id",
    "term_name",
    "ontology",
    "source",
    "nes",
    "pvalue",
    "padj",
    "leading_edge",
    "direction",
    "engine",
]


INLINE_PREFIX_RE = re.compile(r"^(genes?|symbols?|gene_list)\s*[:=]\s*", flags=re.IGNORECASE)
INLINE_WRAP_RE = re.compile(r"^[\[\]\(\)\{\}'\"“”‘’`]+|[\[\]\(\)\{\}'\"“”‘’`]+$")


def timestamp_tag() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_output_dir(outdir: str | None, runs_root: str = "runs") -> Path:
    """Create output directory. Defaults to runs/<timestamp>."""
    if outdir:
        target = Path(outdir).expanduser().resolve()
    else:
        target = Path.cwd() / runs_root / timestamp_tag()
    target.mkdir(parents=True, exist_ok=True)
    return target


def _clean_inline_token(token: str) -> str:
    t = (token or "").strip()
    if not t:
        return ""
    t = INLINE_PREFIX_RE.sub("", t)
    prev = None
    while prev != t:
        prev = t
        t = INLINE_WRAP_RE.sub("", t).strip()
    if t.lower() in {"genes", "gene", "symbols", "gene_list"}:
        return ""
    return t


def _split_inline_values(raw: str) -> List[str]:
    cleaned_raw = INLINE_PREFIX_RE.sub("", (raw or "").strip())
    vals = []
    for v in re.split(r"[\s,;]+", cleaned_raw):
        x = _clean_inline_token(v)
        if x:
            vals.append(x)
    return vals


def parse_gene_input(raw: str) -> List[str]:
    """Parse genes from file path or inline list."""
    p = Path(raw).expanduser()
    if p.exists() and p.is_file():
        values: List[str] = []
        with p.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = re.split(r"[\t,\s]+", line)
                if parts and parts[0].strip():
                    token = _clean_inline_token(parts[0].strip())
                    if token:
                        values.append(token)
        return dedupe_keep_order(values)
    return dedupe_keep_order(_split_inline_values(raw))


def parse_score_input(raw: str | None) -> List[Tuple[str, float]]:
    """Parse ranked list from file path or inline pairs.

    Inline format supports:
      - "GENE1:2.1,GENE2:-1.9"
      - "GENE1=2.1 GENE2=-1.9"
    File format supports two-column delimited records: gene<tab|,|space>score
    """
    if not raw:
        return []

    p = Path(raw).expanduser()
    records: List[Tuple[str, float]] = []

    if p.exists() and p.is_file():
        with p.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = re.split(r"[\t,\s]+", line)
                if len(parts) < 2:
                    continue
                try:
                    gene = _clean_inline_token(parts[0].strip())
                    if not gene:
                        continue
                    records.append((gene, float(parts[1])))
                except ValueError:
                    continue
        return dedupe_ranked(records)

    for item in re.split(r"[,;\s]+", raw.strip()):
        if not item:
            continue
        pair = re.split(r"[:=]", item, maxsplit=1)
        if len(pair) != 2:
            continue
        gene = _clean_inline_token(pair[0].strip())
        try:
            score = float(pair[1])
        except ValueError:
            continue
        if gene:
            records.append((gene, score))
    return dedupe_ranked(records)


def parse_ontology(value: str | None) -> List[str]:
    if not value or value.strip().upper() == "ALL":
        return ["BP", "MF", "CC"]
    allowed = {"BP", "MF", "CC"}
    parsed = []
    for x in re.split(r"[,\s]+", value.strip().upper()):
        if x in allowed and x not in parsed:
            parsed.append(x)
    return parsed or ["BP", "MF", "CC"]


def dedupe_keep_order(values: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def dedupe_ranked(values: Sequence[Tuple[str, float]]) -> List[Tuple[str, float]]:
    seen = set()
    out: List[Tuple[str, float]] = []
    for gene, score in values:
        if gene not in seen:
            seen.add(gene)
            out.append((gene, score))
    return out


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        reader = csv.DictReader(f, delimiter=delimiter)
        rows = list(reader)
        # Backward compatibility: older normalized_input.tsv files were comma-delimited.
        if delimiter == "\t" and reader.fieldnames and len(reader.fieldnames) == 1 and "," in reader.fieldnames[0]:
            f.seek(0)
            rows = list(csv.DictReader(f, delimiter=","))
        return rows


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def contains_cjk(text: str) -> bool:
    return bool(re.search(r"[\u4e00-\u9fff]", text))


def resolve_language(language: str, context: str = "") -> str:
    language = (language or "auto").lower()
    if language in {"zh", "en"}:
        return language
    if contains_cjk(context):
        return "zh"
    env_lang = os.environ.get("LANG", "")
    if env_lang.lower().startswith("zh"):
        return "zh"
    return "en"


def write_placeholder_png(path: Path) -> None:
    """Write a tiny valid PNG placeholder when plotting stack is unavailable."""
    png_bytes = bytes(
        [
            137,
            80,
            78,
            71,
            13,
            10,
            26,
            10,
            0,
            0,
            0,
            13,
            73,
            72,
            68,
            82,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            1,
            8,
            6,
            0,
            0,
            0,
            31,
            21,
            196,
            137,
            0,
            0,
            0,
            13,
            73,
            68,
            65,
            84,
            120,
            156,
            99,
            248,
            207,
            192,
            0,
            0,
            3,
            1,
            1,
            0,
            24,
            221,
            141,
            205,
            0,
            0,
            0,
            0,
            73,
            69,
            78,
            68,
            174,
            66,
            96,
            130,
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(png_bytes)
