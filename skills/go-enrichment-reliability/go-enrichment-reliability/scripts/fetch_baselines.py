#!/usr/bin/env python3
"""Fetch baseline enrichment outputs from online and local reference tools."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

from enrichment_common import GSEA_FIELDS, ORA_FIELDS, parse_gene_input, parse_ontology, parse_score_input, write_csv
from normalize_gene_ids import normalize_genes

ONTOLOGY_SOURCE_MAP = {
    "BP": "GO:BP",
    "MF": "GO:MF",
    "CC": "GO:CC",
}

ONTOLOGY_LIBRARY_MAP = {
    "BP": "GO_Biological_Process_2023",
    "MF": "GO_Molecular_Function_2023",
    "CC": "GO_Cellular_Component_2023",
}


def _as_scalar(v):
    try:
        if hasattr(v, "iloc"):
            return v.iloc[-1] if len(v) else ""
    except Exception:
        pass
    return v


def _parse_term_value(raw) -> tuple[str, str]:
    text = str(_as_scalar(raw) or "").strip()
    if "dtype:" in text and "GO:" in text:
        m = re.search(r"([^\n]+\(GO:\d+\))", text)
        if m:
            text = m.group(1).strip()

    term_id = ""
    term_name = text
    m = re.search(r"\((GO:\d+)\)\s*$", text, flags=re.IGNORECASE)
    if m:
        term_id = m.group(1).upper()
        term_name = text[: m.start()].strip()
    return term_id, term_name


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Fetch baseline results for consistency checks")
    p.add_argument("--genes", required=True)
    p.add_argument("--scores")
    p.add_argument("--species", choices=["human", "mouse"], default="human")
    p.add_argument("--ontology", default="ALL")
    p.add_argument("--sig-threshold", type=float, required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--benchmark-mode", action="store_true")
    return p


def _save_empty_outputs(outdir: Path) -> None:
    write_csv(outdir / "baseline_enrichr.csv", ORA_FIELDS, [])
    write_csv(outdir / "baseline_gprofiler.csv", ORA_FIELDS, [])
    write_csv(outdir / "baseline_clusterprofiler.csv", ORA_FIELDS, [])
    write_csv(outdir / "baseline_broad_gsea.csv", GSEA_FIELDS, [])
    write_csv(outdir / "baseline_clusterprofiler_gsea.csv", GSEA_FIELDS, [])


def fetch_enrichr(genes: List[str], ontologies: List[str], sig_threshold: float) -> Tuple[List[dict], str | None]:
    try:
        import requests  # type: ignore

        add_resp = requests.post(
            "https://maayanlab.cloud/Enrichr/addList",
            files={"list": (None, "\n".join(genes)), "description": (None, "go-enrichment-skill")},
            timeout=30,
        )
        add_resp.raise_for_status()
        user_list_id = add_resp.json().get("userListId")
        if not user_list_id:
            return [], "Enrichr addList returned no userListId"

        rows = []
        for ont in ontologies:
            library = ONTOLOGY_LIBRARY_MAP[ont]
            enrich_resp = requests.get(
                "https://maayanlab.cloud/Enrichr/enrich",
                params={"userListId": user_list_id, "backgroundType": library},
                timeout=30,
            )
            enrich_resp.raise_for_status()
            payload = enrich_resp.json().get(library, [])
            for item in payload:
                # Expected format: [rank, term_name, p, z, combined, genes, adj_p, old_p, old_adj_p]
                term_name = str(item[1]) if len(item) > 1 else ""
                term_id = ""
                m = re.search(r"\((GO:\d+)\)\s*$", term_name, flags=re.IGNORECASE)
                if m:
                    term_id = m.group(1).upper()
                    term_name = term_name[: m.start()].strip()
                pval = item[2] if len(item) > 2 else ""
                padj = item[6] if len(item) > 6 else ""
                gene_hits = len(str(item[5]).split(";")) if len(item) > 5 and str(item[5]).strip() else 0
                rows.append(
                    {
                        "term_id": term_id,
                        "term_name": term_name,
                        "ontology": ont,
                        "source": "enrichr",
                        "pvalue": pval,
                        "padj": padj,
                        "gene_hits": gene_hits,
                        "gene_ratio": "",
                        "engine": "baseline",
                    }
                )

        # Filter if padj is numeric.
        filt_rows = []
        for r in rows:
            try:
                if float(r["padj"]) <= sig_threshold:
                    filt_rows.append(r)
            except Exception:
                filt_rows.append(r)
        return filt_rows, None
    except Exception as exc:  # pragma: no cover - network-dependent
        return [], f"Enrichr baseline failed: {exc}"


def fetch_gprofiler(genes: List[str], species: str, ontologies: List[str], sig_threshold: float) -> Tuple[List[dict], str | None]:
    organism = "hsapiens" if species == "human" else "mmusculus"
    sources = [ONTOLOGY_SOURCE_MAP[o] for o in ontologies]

    # Prefer direct HTTP API to avoid hard dependency on python wrapper.
    try:
        import requests  # type: ignore

        payload = {
            "organism": organism,
            "query": genes,
            "sources": sources,
            "user_threshold": sig_threshold,
            "significance_threshold_method": "fdr",
        }
        resp = requests.post(
            "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
            data=json.dumps(payload),
            headers={"Content-Type": "application/json"},
            timeout=40,
        )
        resp.raise_for_status()
        result = resp.json().get("result", [])

        rows = []
        for r in result:
            ont = str(r.get("source", "")).replace("GO:", "")
            if ont not in ontologies:
                continue
            intersection_size = int(r.get("intersection_size", 0) or 0)
            query_size = int(r.get("query_size", len(genes)) or len(genes))
            rows.append(
                {
                    "term_id": str(r.get("native", "")),
                    "term_name": str(r.get("name", "")),
                    "ontology": ont,
                    "source": "gprofiler",
                    "pvalue": r.get("p_value", ""),
                    "padj": r.get("p_value", ""),
                    "gene_hits": intersection_size,
                    "gene_ratio": f"{intersection_size}/{max(query_size, 1)}",
                    "engine": "baseline",
                }
            )
        return rows, None
    except Exception as exc:  # pragma: no cover - network-dependent
        return [], f"g:Profiler baseline failed: {exc}"


def ensure_normalized_input(outdir: Path, genes: List[str], species: str) -> Path:
    normalized_input = outdir / "normalized_input.tsv"
    if normalized_input.exists():
        return normalized_input

    norm = normalize_genes(genes, species)
    write_csv(
        normalized_input,
        ["input_id", "detected_input_type", "symbol", "entrez", "ensembl", "status", "note"],
        norm["rows"],
    )
    return normalized_input


def ensure_scores_file(outdir: Path, parsed_scores: List[tuple[str, float]]) -> Path | None:
    if not parsed_scores:
        return None
    path = outdir / "scores_input.tsv"
    if path.exists():
        return path
    with path.open("w", encoding="utf-8") as f:
        for gene, score in parsed_scores:
            f.write(f"{gene}\t{score}\n")
    return path


def fetch_clusterprofiler(
    outdir: Path,
    normalized_input: Path,
    scores_path: Path | None,
    species: str,
    ontology: str,
    sig_threshold: float,
) -> Tuple[List[dict], List[dict], str | None]:
    script_dir = Path(__file__).resolve().parent
    r_script = script_dir / "run_go_enrichment_r.R"
    ora_tmp = outdir / ".tmp_clusterprofiler_ora.csv"
    gsea_tmp = outdir / ".tmp_clusterprofiler_gsea.csv"

    cmd = [
        "Rscript",
        str(r_script),
        "--normalized-input",
        str(normalized_input),
        "--species",
        species,
        "--ontology",
        ontology,
        "--sig-threshold",
        str(sig_threshold),
        "--mode",
        "gsea" if scores_path else "ora",
        "--ora-out",
        str(ora_tmp),
        "--gsea-out",
        str(gsea_tmp),
    ]
    if scores_path:
        cmd.extend(["--scores", str(scores_path)])

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            err = proc.stderr.strip() or proc.stdout.strip() or f"clusterProfiler failed ({proc.returncode})"
            return [], [], err

        ora_rows = []
        gsea_rows = []
        if ora_tmp.exists():
            import csv

            with ora_tmp.open("r", encoding="utf-8", newline="") as f:
                ora_rows = list(csv.DictReader(f))
        if gsea_tmp.exists():
            import csv

            with gsea_tmp.open("r", encoding="utf-8", newline="") as f:
                gsea_rows = list(csv.DictReader(f))

        return ora_rows, gsea_rows, None
    finally:
        if ora_tmp.exists():
            ora_tmp.unlink()
        if gsea_tmp.exists():
            gsea_tmp.unlink()


def fetch_broad_gsea(scores_path: Path | None, ontologies: List[str], sig_threshold: float) -> Tuple[List[dict], str | None]:
    if not scores_path:
        return [], None

    # Optional integration: if user sets BROAD_GSEA_JAR and BROAD_GSEA_GMT, run CLI.
    gsea_jar = os.environ.get("BROAD_GSEA_JAR", "").strip()
    gmt_path = os.environ.get("BROAD_GSEA_GMT", "").strip()
    if not gsea_jar or not gmt_path:
        return [], "Broad GSEA CLI not configured (set BROAD_GSEA_JAR and BROAD_GSEA_GMT)"

    # Lightweight fallback: run gseapy prerank as a proxy if CLI is unavailable.
    # If CLI is configured, users can patch this section for full Broad compatibility.
    try:
        import pandas as pd  # type: ignore
        import gseapy as gp  # type: ignore
        perm_num = max(10, int(os.environ.get("BROAD_PROXY_PERM_NUM", "100")))
        threads = max(1, int(os.environ.get("BROAD_PROXY_THREADS", "1")))

        rank_df = pd.read_csv(scores_path, sep="\t", header=None, names=["gene", "score"])
        rank_df = rank_df.drop_duplicates(subset=["gene"]).dropna(subset=["score"])
        rows = []

        for ont in ontologies:
            pre = gp.prerank(
                rnk=rank_df,
                gene_sets=gmt_path,
                min_size=5,
                max_size=500,
                permutation_num=perm_num,
                threads=threads,
                outdir=None,
                seed=7,
                no_plot=True,
            )
            res = getattr(pre, "res2d", None)
            if res is None or len(res) == 0:
                continue
            if "Term" in getattr(res, "columns", []):
                res = res.reset_index(drop=True)
            else:
                res = res.reset_index().rename(columns={"index": "Term"})
            for _, r in res.iterrows():
                term_id, term_name = _parse_term_value(r.get("Term", ""))
                nes = float(r.get("NES", 0.0) or 0.0)
                row = {
                    "term_id": term_id,
                    "term_name": term_name,
                    "ontology": ont,
                    "source": "broad_gsea_proxy",
                    "nes": nes,
                    "pvalue": _as_scalar(r.get("NOM p-val", "")),
                    "padj": _as_scalar(r.get("FDR q-val", "")),
                    "leading_edge": str(_as_scalar(r.get("Lead_genes", "")) or ""),
                    "direction": "up" if nes >= 0 else "down",
                    "engine": "baseline",
                }
                rows.append(row)

        filtered = []
        for r in rows:
            try:
                if float(r["padj"]) <= sig_threshold:
                    filtered.append(r)
            except Exception:
                filtered.append(r)
        return filtered, None
    except Exception as exc:
        return [], f"Broad GSEA baseline unavailable: {exc}"


def main() -> int:
    args = build_parser().parse_args()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    _save_empty_outputs(outdir)

    genes = parse_gene_input(args.genes)
    ontologies = parse_ontology(args.ontology)
    scores = parse_score_input(args.scores)

    errors: List[str] = []

    enrichr_rows, enrichr_err = fetch_enrichr(genes, ontologies, args.sig_threshold)
    write_csv(outdir / "baseline_enrichr.csv", ORA_FIELDS, enrichr_rows)
    if enrichr_err:
        errors.append(enrichr_err)

    gprof_rows, gprof_err = fetch_gprofiler(genes, args.species, ontologies, args.sig_threshold)
    write_csv(outdir / "baseline_gprofiler.csv", ORA_FIELDS, gprof_rows)
    if gprof_err:
        errors.append(gprof_err)

    normalized_input = ensure_normalized_input(outdir, genes, args.species)
    scores_path = ensure_scores_file(outdir, scores)
    cp_ora_rows, cp_gsea_rows, cp_err = fetch_clusterprofiler(
        outdir,
        normalized_input,
        scores_path,
        args.species,
        args.ontology,
        args.sig_threshold,
    )
    write_csv(outdir / "baseline_clusterprofiler.csv", ORA_FIELDS, cp_ora_rows)
    write_csv(outdir / "baseline_clusterprofiler_gsea.csv", GSEA_FIELDS, cp_gsea_rows)
    if cp_err:
        errors.append(cp_err)

    broad_rows, broad_err = fetch_broad_gsea(scores_path, ontologies, args.sig_threshold)
    write_csv(outdir / "baseline_broad_gsea.csv", GSEA_FIELDS, broad_rows)
    if broad_err:
        errors.append(broad_err)

    if errors:
        for e in errors:
            print(e)
        if args.benchmark_mode:
            return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
