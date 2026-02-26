#!/usr/bin/env python3
"""Normalize gene identifiers (SYMBOL/ENTREZ/ENSEMBL) for human and mouse."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

from enrichment_common import parse_gene_input, write_csv, write_json

ID_RE_ENSEMBL = re.compile(r"^(ENSG\d+|ENSMUSG\d+)(?:\.\d+)?$", re.IGNORECASE)
ID_RE_DIGITS = re.compile(r"^\d+$")


def detect_input_id_type(values: List[str]) -> str:
    if not values:
        return "UNKNOWN"
    total = len(values)
    ensembl_hits = sum(1 for v in values if ID_RE_ENSEMBL.match(v))
    digit_hits = sum(1 for v in values if ID_RE_DIGITS.match(v))

    if ensembl_hits / total >= 0.7:
        return "ENSEMBL"
    if digit_hits / total >= 0.7:
        return "ENTREZ"
    return "SYMBOL"


def _extract_ensembl(value):
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, list) and value:
        head = value[0]
        if isinstance(head, dict):
            return str(head.get("gene", ""))
        return str(head)
    if isinstance(value, dict):
        return str(value.get("gene", ""))
    return ""


def _normalize_with_mygene(genes: List[str], species: str, detected_input: str) -> List[dict]:
    # Delayed import to keep script functional even without optional dependency.
    import mygene  # type: ignore

    mg = mygene.MyGeneInfo()
    scopes = {
        "SYMBOL": "symbol",
        "ENTREZ": "entrezgene",
        "ENSEMBL": "ensembl.gene",
    }.get(detected_input, "symbol,entrezgene,ensembl.gene")

    query_res = mg.querymany(
        genes,
        scopes=scopes,
        species="human" if species == "human" else "mouse",
        fields="symbol,entrezgene,ensembl.gene",
        as_dataframe=False,
        returnall=False,
        verbose=False,
    )

    by_query: Dict[str, dict] = {}
    for item in query_res:
        q = str(item.get("query", ""))
        if q and q not in by_query:
            by_query[q] = item

    rows = []
    for g in genes:
        payload = by_query.get(g, {})
        symbol = str(payload.get("symbol", "") or "")
        entrez = str(payload.get("entrezgene", "") or "")
        ensembl = _extract_ensembl(payload.get("ensembl"))
        status = "mapped" if (symbol or entrez or ensembl) else "unmapped"

        # Fill a usable fallback for partially-mapped or missing records.
        if not symbol and detected_input == "SYMBOL":
            symbol = g
        if not entrez and detected_input == "ENTREZ":
            entrez = g
        if not ensembl and detected_input == "ENSEMBL":
            ensembl = g

        rows.append(
            {
                "input_id": g,
                "detected_input_type": detected_input,
                "symbol": symbol,
                "entrez": entrez,
                "ensembl": ensembl,
                "status": status,
                "note": "",
            }
        )
    return rows


def _normalize_fallback(genes: List[str], detected_input: str) -> List[dict]:
    rows = []
    note = "mygene package unavailable; used heuristic pass-through mapping"
    for g in genes:
        symbol = g if detected_input == "SYMBOL" else ""
        entrez = g if detected_input == "ENTREZ" else ""
        ensembl = g if detected_input == "ENSEMBL" else ""
        rows.append(
            {
                "input_id": g,
                "detected_input_type": detected_input,
                "symbol": symbol,
                "entrez": entrez,
                "ensembl": ensembl,
                "status": "mapped" if (symbol or entrez or ensembl) else "unmapped",
                "note": note,
            }
        )
    return rows


def normalize_genes(genes: List[str], species: str) -> dict:
    detected = detect_input_id_type(genes)
    try:
        rows = _normalize_with_mygene(genes, species, detected)
        method = "mygene"
    except Exception as exc:  # pragma: no cover - dependency/network dependent
        rows = _normalize_fallback(genes, detected)
        method = f"fallback:{type(exc).__name__}"

    mapped = sum(1 for r in rows if r["status"] == "mapped")
    return {
        "rows": rows,
        "stats": {
            "total": len(rows),
            "mapped": mapped,
            "unmapped": len(rows) - mapped,
            "detected_input_type": detected,
            "method": method,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Normalize gene IDs for enrichment runs")
    p.add_argument("--genes", required=True, help="Path or inline list of genes")
    p.add_argument("--species", choices=["human", "mouse"], default="human")
    p.add_argument("--out", required=True, help="Output TSV path")
    p.add_argument("--meta-out", help="Optional JSON metadata path")
    return p


def main() -> int:
    args = build_parser().parse_args()
    genes = parse_gene_input(args.genes)
    result = normalize_genes(genes, args.species)

    out_path = Path(args.out).expanduser().resolve()
    write_csv(
        out_path,
        ["input_id", "detected_input_type", "symbol", "entrez", "ensembl", "status", "note"],
        result["rows"],
    )

    if args.meta_out:
        write_json(Path(args.meta_out).expanduser().resolve(), result["stats"])

    print(f"Normalized {result['stats']['total']} genes ({result['stats']['mapped']} mapped).")
    if result["stats"]["unmapped"]:
        print("Hint: provide SYMBOL (TP53), ENTREZ (7157), or ENSEMBL (ENSG...) IDs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
