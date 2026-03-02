#!/usr/bin/env python3
"""Primary GO/GSEA entrypoint with Python-first and R fallback strategy."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

from enrichment_common import (
    GSEA_FIELDS,
    ORA_FIELDS,
    ensure_output_dir,
    parse_gene_input,
    parse_ontology,
    parse_score_input,
    read_csv,
    resolve_language,
    write_csv,
    write_json,
    write_placeholder_png,
)
from normalize_gene_ids import normalize_genes


ONTOLOGY_LIBRARY_MAP = {
    "BP": "GO_Biological_Process_2023",
    "MF": "GO_Molecular_Function_2023",
    "CC": "GO_Cellular_Component_2023",
}

ONTOLOGY_SOURCE_MAP = {
    "BP": "GO:BP",
    "MF": "GO:MF",
    "CC": "GO:CC",
}

ONTOLOGY_COLORS = {
    "BP": "#4C78A8",
    "CC": "#54A24B",
    "MF": "#F58518",
    "UNKNOWN": "#9E9E9E",
}

THREAD_ENV_VARS = [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run GO ORA/GSEA with reliability outputs")
    p.add_argument("--task", default="go_enrichment")
    p.add_argument("--genes", required=True, help="Path or inline gene list")
    p.add_argument("--sig-threshold", type=float, required=True)
    p.add_argument("--species", choices=["human", "mouse"], default="human")
    p.add_argument("--ontology", default="ALL")
    p.add_argument("--scores", help="Path or inline ranked gene scores")
    p.add_argument("--mode", choices=["auto", "ora", "gsea"], default="auto")
    p.add_argument("--engine", choices=["auto", "python", "r"], default="auto")
    p.add_argument("--background", help="Optional background gene file path")
    p.add_argument("--language", choices=["zh", "en", "auto"], default="auto")
    p.add_argument("--outdir", help="Output directory; default runs/<timestamp>")
    p.add_argument("--benchmark-mode", action="store_true", help="Fail on baseline collection errors")
    p.add_argument("--context-text", default="", help="Optional text for auto language detection")
    p.add_argument("--skip-baselines", action="store_true")
    p.add_argument("--celltype-candidates", default="", help="Optional cell-type focus list, comma-separated")
    p.add_argument("--celltype-topk", type=int, default=5)
    p.add_argument("--skip-celltype-link", action="store_true")
    return p


def configure_runtime_threads() -> int:
    raw = os.environ.get("GO_ENRICH_THREADS", "").strip() or os.environ.get("GSEAPY_THREADS", "").strip() or "1"
    try:
        threads = max(1, int(raw))
    except ValueError:
        threads = 1

    for key in THREAD_ENV_VARS:
        os.environ.setdefault(key, str(threads))
    os.environ.setdefault("GSEAPY_THREADS", str(threads))
    return threads


def load_background(path: Optional[str]) -> List[str]:
    if not path:
        return []
    p = Path(path).expanduser()
    if not p.exists():
        return []
    return parse_gene_input(str(p))


def normalize_and_write_input(genes: List[str], species: str, outdir: Path) -> Dict:
    norm = normalize_genes(genes, species)
    normalized_path = outdir / "normalized_input.tsv"
    write_csv(
        normalized_path,
        ["input_id", "detected_input_type", "symbol", "entrez", "ensembl", "status", "note"],
        norm["rows"],
    )
    return norm


def _symbol_genes(normalized_rows: List[dict], fallback_genes: List[str]) -> List[str]:
    symbols = [r.get("symbol", "").strip() for r in normalized_rows if r.get("symbol", "").strip()]
    if symbols:
        # Keep order while deduping.
        seen = set()
        out = []
        for s in symbols:
            if s not in seen:
                seen.add(s)
                out.append(s)
        return out
    return fallback_genes


def _as_scalar(v):
    try:
        if hasattr(v, "iloc"):
            return v.iloc[-1] if len(v) else ""
    except Exception:
        pass
    return v


def _parse_term_value(raw) -> tuple[str, str]:
    text = str(_as_scalar(raw) or "").strip()
    # gseapy can emit stringified pandas Series when duplicate "Term" columns exist.
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


def run_python_ora(
    genes: List[str],
    species: str,
    ontologies: List[str],
    sig_threshold: float,
    background: List[str],
) -> tuple[List[dict], List[dict], str]:
    rows: List[dict] = []
    plot_rows: List[dict] = []
    organism = "hsapiens" if species == "human" else "mmusculus"
    backend = ""

    def _gprofiler_to_rows(df) -> List[dict]:
        parsed = []
        for _, r in df.iterrows():
            ont = str(r.get("source", "")).replace("GO:", "")
            if ont not in ontologies:
                continue
            intersection_size = int(r.get("intersection_size", 0) or 0)
            query_size = int(r.get("query_size", len(genes)) or len(genes))
            parsed.append(
                {
                    "term_id": str(r.get("native", "")),
                    "term_name": str(r.get("name", "")),
                    "ontology": ont,
                    "source": "gprofiler",
                    "pvalue": r.get("p_value", ""),
                    "padj": r.get("p_value", ""),
                    "gene_hits": intersection_size,
                    "gene_ratio": f"{intersection_size}/{max(query_size, 1)}",
                    "engine": "python",
                }
            )
        return parsed

    try:
        from gprofiler import GProfiler  # type: ignore

        gp = GProfiler(return_dataframe=True)
        sources = [ONTOLOGY_SOURCE_MAP[o] for o in ontologies]
        kwargs = {
            "organism": organism,
            "query": genes,
            "sources": sources,
            "user_threshold": sig_threshold,
            "significance_threshold_method": "fdr",
        }
        if background:
            kwargs["background"] = background
        df = gp.profile(**kwargs)

        if df is not None:
            backend = "gprofiler"
        if df is not None and len(df) > 0:
            rows = _gprofiler_to_rows(df)
            plot_rows = list(rows)
            return rows, plot_rows, backend

        # If thresholded result is empty, fetch nominal top terms for plotting context.
        kwargs_nominal = dict(kwargs)
        kwargs_nominal["user_threshold"] = 1.0
        df_nominal = gp.profile(**kwargs_nominal)
        if df_nominal is not None and len(df_nominal) > 0:
            plot_rows = _gprofiler_to_rows(df_nominal)[:30]
    except Exception:
        pass

    # Secondary fallback path for ORA using Enrichr libraries through gseapy.
    try:
        import gseapy as gp  # type: ignore
        backend = backend or "gseapy-enrichr"

        for ont in ontologies:
            library = ONTOLOGY_LIBRARY_MAP[ont]
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=library,
                outdir=None,
                cutoff=sig_threshold,
                background=background if background else None,
                no_plot=True,
            )
            result_df = getattr(enr, "results", None)
            if result_df is None or len(result_df) == 0:
                continue
            for _, r in result_df.iterrows():
                term = str(r.get("Term", ""))
                term_id = ""
                term_name = term
                if "(" in term and term.endswith(")"):
                    maybe_name, maybe_id = term.rsplit("(", 1)
                    term_name = maybe_name.strip()
                    term_id = maybe_id[:-1].strip()
                overlap = str(r.get("Overlap", "0/1"))
                numerator = overlap.split("/")[0] if "/" in overlap else overlap
                try:
                    gene_hits = int(numerator)
                except ValueError:
                    gene_hits = 0
                rows.append(
                    {
                        "term_id": term_id,
                        "term_name": term_name,
                        "ontology": ont,
                        "source": "enrichr",
                        "pvalue": r.get("P-value", ""),
                        "padj": r.get("Adjusted P-value", ""),
                        "gene_hits": gene_hits,
                        "gene_ratio": overlap,
                        "engine": "python",
                    }
                )

        if rows:
            plot_rows = list(rows)
        else:
            # Try nominal terms for plotting only.
            for ont in ontologies:
                library = ONTOLOGY_LIBRARY_MAP[ont]
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=library,
                    outdir=None,
                    cutoff=1.0,
                    background=background if background else None,
                    no_plot=True,
                )
                result_df = getattr(enr, "results", None)
                if result_df is None or len(result_df) == 0:
                    continue
                for _, r in result_df.head(30).iterrows():
                    term = str(r.get("Term", ""))
                    term_id = ""
                    term_name = term
                    if "(" in term and term.endswith(")"):
                        maybe_name, maybe_id = term.rsplit("(", 1)
                        term_name = maybe_name.strip()
                        term_id = maybe_id[:-1].strip()
                    overlap = str(r.get("Overlap", "0/1"))
                    numerator = overlap.split("/")[0] if "/" in overlap else overlap
                    try:
                        gene_hits = int(numerator)
                    except ValueError:
                        gene_hits = 0
                    plot_rows.append(
                        {
                            "term_id": term_id,
                            "term_name": term_name,
                            "ontology": ont,
                            "source": "enrichr",
                            "pvalue": r.get("P-value", ""),
                            "padj": r.get("Adjusted P-value", ""),
                            "gene_hits": gene_hits,
                            "gene_ratio": overlap,
                            "engine": "python",
                        }
                    )
    except Exception:
        pass

    if not backend:
        raise RuntimeError("No Python ORA backend available (install gprofiler-official or gseapy)")
    return rows, plot_rows, backend


def run_python_gsea(
    scores: List[tuple[str, float]],
    ontologies: List[str],
    sig_threshold: float,
) -> tuple[List[dict], str]:
    rows: List[dict] = []
    if not scores:
        return rows, "skipped"
    backend = ""

    try:
        import pandas as pd  # type: ignore
        import gseapy as gp  # type: ignore
        backend = "gseapy-prerank"
        perm_num = max(10, int(os.environ.get("GSEAPY_PERM_NUM", "100")))
        threads = max(1, int(os.environ.get("GSEAPY_THREADS", "1")))

        rank_df = pd.DataFrame(scores, columns=["gene", "score"]).drop_duplicates(subset=["gene"])

        for ont in ontologies:
            library = ONTOLOGY_LIBRARY_MAP[ont]
            prerank_res = gp.prerank(
                rnk=rank_df,
                gene_sets=library,
                min_size=5,
                max_size=500,
                permutation_num=perm_num,
                threads=threads,
                outdir=None,
                seed=7,
                no_plot=True,
            )
            res_df = getattr(prerank_res, "res2d", None)
            if res_df is None or len(res_df) == 0:
                continue
            if "Term" in getattr(res_df, "columns", []):
                res_df = res_df.reset_index(drop=True)
            else:
                res_df = res_df.reset_index().rename(columns={"index": "Term"})
            for _, r in res_df.iterrows():
                term_id, term_name = _parse_term_value(r.get("Term", ""))
                nes = float(r.get("NES", 0.0) or 0.0)
                fdr = _as_scalar(r.get("FDR q-val", ""))
                pval = _as_scalar(r.get("NOM p-val", ""))
                lead = str(_as_scalar(r.get("Lead_genes", "")) or "")
                rows.append(
                    {
                        "term_id": term_id,
                        "term_name": term_name,
                        "ontology": ont,
                        "source": "gseapy",
                        "nes": nes,
                        "pvalue": pval,
                        "padj": fdr,
                        "leading_edge": lead,
                        "direction": "up" if nes >= 0 else "down",
                        "engine": "python",
                    }
                )
    except Exception:
        backend = ""

    # Filter by threshold if possible.
    filtered = []
    for r in rows:
        try:
            padj = float(r.get("padj", "nan"))
            if math.isnan(padj) or padj <= sig_threshold:
                filtered.append(r)
        except Exception:
            filtered.append(r)
    if not backend:
        raise RuntimeError("No Python GSEA backend available (install pandas and gseapy)")
    return filtered, backend


def create_plots(ora_rows: List[dict], outdir: Path, plot_fallback_rows: Optional[List[dict]] = None) -> tuple[bool, str]:
    dot_path = outdir / "plots_dotplot.png"
    bar_path = outdir / "plots_barplot.png"

    try:
        import matplotlib.pyplot as plt  # type: ignore
        from matplotlib.lines import Line2D  # type: ignore
        from matplotlib.patches import Patch  # type: ignore

        source_rows = ora_rows if ora_rows else (plot_fallback_rows or [])
        top = sorted(
            source_rows,
            key=lambda r: float(r.get("padj") or 1.0),
        )[:15]

        if not top:
            write_placeholder_png(dot_path)
            write_placeholder_png(bar_path)
            return False, "No ORA terms available for plotting"

        names = [r.get("term_name", "")[:48] for r in top]
        x_vals = []
        y_vals = []
        sizes = []
        ontologies = []
        for row in top:
            padj = float(row.get("padj") or 1.0)
            hits = float(row.get("gene_hits") or 1)
            ont = str(row.get("ontology") or "").strip().upper()
            if ont.startswith("GO:"):
                ont = ont[3:]
            if ont not in {"BP", "CC", "MF"}:
                ont = "UNKNOWN"
            x_vals.append(-math.log10(max(padj, 1e-300)))
            y_vals.append(hits)
            sizes.append(max(hits * 12.0, 20.0))
            ontologies.append(ont)

        point_colors = [ONTOLOGY_COLORS.get(ont, ONTOLOGY_COLORS["UNKNOWN"]) for ont in ontologies]
        ontology_order = ["BP", "CC", "MF", "UNKNOWN"]
        present_ontologies = [ont for ont in ontology_order if ont in set(ontologies)]
        label_map = {"BP": "BP", "CC": "CC", "MF": "MF", "UNKNOWN": "Other"}

        plt.figure(figsize=(8, 5))
        plt.scatter(x_vals, y_vals, s=sizes, alpha=0.78, c=point_colors, edgecolors="#1F2D3D", linewidths=0.6)
        plt.xlabel("-log10(FDR)")
        plt.ylabel("Gene hits")
        plt.title("GO Enrichment Dotplot")
        if present_ontologies:
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    markerfacecolor=ONTOLOGY_COLORS[ont],
                    markeredgecolor="#1F2D3D",
                    markersize=8,
                    label=label_map[ont],
                )
                for ont in present_ontologies
            ]
            plt.legend(handles=handles, title="Ontology", frameon=False, loc="best")
        plt.tight_layout()
        plt.savefig(dot_path, dpi=180)
        plt.close()

        bar_colors = [ONTOLOGY_COLORS.get(ont, ONTOLOGY_COLORS["UNKNOWN"]) for ont in ontologies[::-1]]
        plt.figure(figsize=(9, 5))
        plt.barh(names[::-1], x_vals[::-1], color=bar_colors)
        plt.xlabel("-log10(FDR)")
        plt.title("GO Enrichment Top Terms")
        if present_ontologies:
            handles = [Patch(facecolor=ONTOLOGY_COLORS[ont], label=label_map[ont]) for ont in present_ontologies]
            plt.legend(handles=handles, title="Ontology", frameon=False, loc="best")
        plt.tight_layout()
        plt.savefig(bar_path, dpi=180)
        plt.close()
    except Exception as exc:
        write_placeholder_png(dot_path)
        write_placeholder_png(bar_path)
        return False, f"Plot rendering fallback used: {exc}"

    if (not ora_rows) and source_rows:
        return True, "No significant ORA terms; plots use nominal top terms for context."
    return True, ""


def _safe_float(value: object, default: float) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _ora_rows_for_plot_from_csv(path: Path, sig_threshold: float) -> tuple[List[dict], List[dict]]:
    all_rows: List[dict] = []
    sig_rows: List[dict] = []
    for row in read_csv(path):
        term_name = str(
            row.get("term_name")
            or row.get("Description")
            or row.get("Term")
            or row.get("ID")
            or ""
        ).strip()
        if not term_name:
            continue

        padj = row.get("padj") or row.get("p.adjust") or row.get("Adjusted P-value") or row.get("FDR") or 1.0
        gene_hits = row.get("gene_hits") or row.get("Count")
        if not gene_hits:
            ratio = str(row.get("gene_ratio") or row.get("GeneRatio") or "")
            if "/" in ratio:
                gene_hits = ratio.split("/", 1)[0]
        plot_row = {
            "term_name": term_name,
            "padj": _safe_float(padj, 1.0),
            "gene_hits": _safe_float(gene_hits, 1.0),
            "ontology": str(row.get("ontology") or row.get("ONTOLOGY") or "").strip(),
        }
        all_rows.append(plot_row)
        if plot_row["padj"] <= sig_threshold:
            sig_rows.append(plot_row)
    return sig_rows, all_rows


def run_r_script(
    script_dir: Path,
    normalized_input: Path,
    scores_path: Optional[Path],
    species: str,
    ontology: str,
    sig_threshold: float,
    mode: str,
    background: Optional[str],
    ora_out: Path,
    gsea_out: Path,
    dotplot_out: Path,
    barplot_out: Path,
) -> tuple[bool, str]:
    r_script = script_dir / "run_go_enrichment_r.R"
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
        mode,
        "--ora-out",
        str(ora_out),
        "--gsea-out",
        str(gsea_out),
        "--dotplot-out",
        str(dotplot_out),
        "--barplot-out",
        str(barplot_out),
    ]
    if scores_path:
        cmd.extend(["--scores", str(scores_path)])
    if background:
        cmd.extend(["--background", background])

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        return True, proc.stdout.strip() or "R engine completed"

    msg = proc.stderr.strip() or proc.stdout.strip() or f"R engine failed with code {proc.returncode}"
    return False, msg


def write_scores_file(scores: List[tuple[str, float]], outdir: Path) -> Optional[Path]:
    if not scores:
        return None
    p = outdir / "scores_input.tsv"
    with p.open("w", encoding="utf-8") as f:
        for gene, score in scores:
            f.write(f"{gene}\t{score}\n")
    return p


def run_baseline_fetch(script_dir: Path, args: argparse.Namespace, outdir: Path, scores_path: Optional[Path]) -> tuple[bool, str]:
    cmd = [
        sys.executable,
        str(script_dir / "fetch_baselines.py"),
        "--genes",
        args.genes,
        "--species",
        args.species,
        "--ontology",
        args.ontology,
        "--sig-threshold",
        str(args.sig_threshold),
        "--outdir",
        str(outdir),
    ]
    if scores_path:
        cmd.extend(["--scores", str(scores_path)])
    if args.benchmark_mode:
        cmd.append("--benchmark-mode")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    ok = proc.returncode == 0
    msg = (proc.stdout + "\n" + proc.stderr).strip()
    return ok, msg


def run_consistency_compare(script_dir: Path, args: argparse.Namespace, outdir: Path) -> tuple[bool, str]:
    cmd = [
        sys.executable,
        str(script_dir / "compare_results.py"),
        "--outdir",
        str(outdir),
        "--ora-threshold",
        "0.6",
        "--gsea-threshold",
        "0.6",
        "--top-n",
        "20",
        "--benchmark-mode",
        "true" if args.benchmark_mode else "false",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    ok = proc.returncode == 0
    msg = (proc.stdout + "\n" + proc.stderr).strip()
    return ok, msg


def run_summary_render(script_dir: Path, outdir: Path, language: str, fallback_note: str, warnings: List[str]) -> None:
    cmd = [
        sys.executable,
        str(script_dir / "render_summary.py"),
        "--outdir",
        str(outdir),
        "--language",
        language,
    ]
    if fallback_note:
        cmd.extend(["--fallback-note", fallback_note])
    for w in warnings:
        cmd.extend(["--warning", w])
    subprocess.run(cmd, capture_output=True, text=True)


def run_analysis_brief(script_dir: Path, outdir: Path, language: str, focus: str) -> tuple[bool, str]:
    cmd = [
        sys.executable,
        str(script_dir / "render_analysis_brief.py"),
        "--outdir",
        str(outdir),
        "--language",
        language,
    ]
    if focus.strip():
        cmd.extend(["--focus", focus.strip()])

    proc = subprocess.run(cmd, capture_output=True, text=True)
    ok = proc.returncode == 0
    msg = (proc.stdout + "\n" + proc.stderr).strip()
    return ok, msg


def run_celltype_link(script_dir: Path, args: argparse.Namespace, outdir: Path, language: str) -> tuple[bool, str]:
    cmd = [
        sys.executable,
        str(script_dir / "infer_celltype_links.py"),
        "--outdir",
        str(outdir),
        "--species",
        args.species,
        "--topk",
        str(max(1, int(args.celltype_topk))),
        "--language",
        language,
        "--context-text",
        args.context_text or "",
    ]
    if args.celltype_candidates.strip():
        cmd.extend(["--candidates", args.celltype_candidates.strip()])

    proc = subprocess.run(cmd, capture_output=True, text=True)
    ok = proc.returncode == 0
    msg = (proc.stdout + "\n" + proc.stderr).strip()
    return ok, msg


def main() -> int:
    args = build_parser().parse_args()
    configure_runtime_threads()
    if args.sig_threshold is None:
        print("sig_threshold is required")
        return 2

    genes = parse_gene_input(args.genes)
    if not genes:
        print("No genes found from --genes input")
        return 2

    scores = parse_score_input(args.scores)
    ontologies = parse_ontology(args.ontology)
    language = resolve_language(args.language, args.context_text)
    outdir = ensure_output_dir(args.outdir)
    script_dir = Path(__file__).resolve().parent

    fallback_note = ""
    requested_gsea = args.mode in {"auto", "gsea"}
    run_gsea = requested_gsea and bool(scores)
    if requested_gsea and not scores:
        fallback_note = "No scores were provided. GSEA was skipped and pipeline downgraded to ORA."

    # Write stable CSV shells so downstream logic always has the expected files.
    py_ora_path = outdir / "ora_results_python.csv"
    py_gsea_path = outdir / "gsea_results_python.csv"
    r_ora_path = outdir / "ora_results_r.csv"
    r_gsea_path = outdir / "gsea_results_r.csv"
    write_csv(py_ora_path, ORA_FIELDS, [])
    write_csv(py_gsea_path, GSEA_FIELDS, [])
    write_csv(r_ora_path, ORA_FIELDS, [])
    write_csv(r_gsea_path, GSEA_FIELDS, [])
    write_csv(outdir / "baseline_enrichr.csv", ORA_FIELDS, [])
    write_csv(outdir / "baseline_gprofiler.csv", ORA_FIELDS, [])
    write_csv(outdir / "baseline_clusterprofiler.csv", ORA_FIELDS, [])
    write_csv(outdir / "baseline_broad_gsea.csv", GSEA_FIELDS, [])
    write_csv(
        outdir / "celltype_association.csv",
        [
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
        ],
        [],
    )
    write_json(outdir / "celltype_association.json", {})
    (outdir / "celltype_interpretation.md").write_text("", encoding="utf-8")
    (outdir / "analysis_brief.md").write_text("", encoding="utf-8")
    write_json(outdir / "analysis_brief.json", {})

    warnings: List[str] = []

    norm = normalize_and_write_input(genes, args.species, outdir)
    normalized_rows = norm["rows"]
    symbol_genes = _symbol_genes(normalized_rows, genes)

    scores_path = write_scores_file(scores, outdir)
    background = load_background(args.background)

    python_ok = False
    python_ora_backend = ""
    python_gsea_backend = ""
    ora_rows_for_plot: List[dict] = []
    if args.engine in {"auto", "python"}:
        try:
            if args.mode in {"auto", "ora", "gsea"}:
                ora_rows, ora_rows_for_plot, python_ora_backend = run_python_ora(
                    genes=symbol_genes,
                    species=args.species,
                    ontologies=ontologies,
                    sig_threshold=args.sig_threshold,
                    background=background,
                )
                write_csv(py_ora_path, ORA_FIELDS, ora_rows)
            else:
                ora_rows = []

            gsea_rows = []
            if run_gsea:
                gsea_rows, python_gsea_backend = run_python_gsea(
                    scores=scores,
                    ontologies=ontologies,
                    sig_threshold=args.sig_threshold,
                )
                write_csv(py_gsea_path, GSEA_FIELDS, gsea_rows)

            plot_ok, plot_note = create_plots(ora_rows, outdir, ora_rows_for_plot)
            if plot_note:
                warnings.append(plot_note)
            if args.mode in {"auto", "ora", "gsea"} and not ora_rows and ora_rows_for_plot:
                warnings.append(
                    "No significant ORA terms under current threshold; plots show nominal top terms for biological context."
                )
            python_ok = True
        except Exception as exc:
            warnings.append(f"Python engine failed: {exc}")

    r_ok = False
    if args.engine == "r" or (args.engine == "auto" and not python_ok):
        dotplot_path = outdir / "plots_dotplot.png"
        barplot_path = outdir / "plots_barplot.png"
        mode_for_r = args.mode
        if args.mode in {"auto", "gsea"} and not scores:
            mode_for_r = "ora"
        r_ok, r_msg = run_r_script(
            script_dir=script_dir,
            normalized_input=outdir / "normalized_input.tsv",
            scores_path=scores_path if run_gsea else None,
            species=args.species,
            ontology=args.ontology,
            sig_threshold=args.sig_threshold,
            mode=mode_for_r,
            background=args.background,
            ora_out=r_ora_path,
            gsea_out=r_gsea_path,
            dotplot_out=dotplot_path,
            barplot_out=barplot_path,
        )
        if not r_ok:
            warnings.append(f"R engine failed: {r_msg}")

        if r_ok and (not dotplot_path.exists() or not barplot_path.exists()):
            r_sig_rows, r_all_rows = _ora_rows_for_plot_from_csv(r_ora_path, args.sig_threshold)
            if r_all_rows:
                plot_ok, plot_note = create_plots(r_sig_rows, outdir, r_all_rows)
                if plot_note:
                    warnings.append(plot_note)

        # Ensure plots exist even if both plotting paths were unavailable.
        if not (outdir / "plots_dotplot.png").exists():
            write_placeholder_png(outdir / "plots_dotplot.png")
        if not (outdir / "plots_barplot.png").exists():
            write_placeholder_png(outdir / "plots_barplot.png")

    no_engine_success = False
    if args.engine == "python" and not python_ok:
        return 3
    if args.engine == "r" and not r_ok:
        return 4
    if args.engine == "auto" and not python_ok and not r_ok:
        warnings.append("Both Python and R engines failed; no enrichment results generated.")
        no_engine_success = True

    baseline_ok = True
    if not args.skip_baselines:
        baseline_ok, baseline_msg = run_baseline_fetch(script_dir, args, outdir, scores_path if run_gsea else None)
        if not baseline_ok:
            if baseline_msg.strip():
                warnings.append(f"Baseline fetch issues: {baseline_msg}")
            else:
                warnings.append("Baseline fetch issues: unknown error")

    if args.skip_baselines:
        consistency_ok = True
        write_json(
            outdir / "consistency_report.json",
            {
                "top_n": 20,
                "thresholds": {
                    "ora_overlap": 0.6,
                    "gsea_direction_consistency": 0.6,
                },
                "run_sources": {
                    "ora": "python" if python_ok else ("r" if r_ok else "none"),
                    "gsea": "python" if (python_ok and run_gsea) else ("r" if (r_ok and run_gsea) else "none"),
                },
                "ora": {"score": None, "pass": True, "details": {}},
                "gsea": {"applicable": bool(run_gsea), "score": None, "pass": True, "details": {}},
                "overall_pass": True,
                "note": "Baselines skipped by --skip-baselines; consistency gate not evaluated.",
            },
        )
    else:
        consistency_ok, consistency_msg = run_consistency_compare(script_dir, args, outdir)
        if not consistency_ok:
            if consistency_msg.strip():
                warnings.append(f"Consistency comparison issues: {consistency_msg}")
            else:
                warnings.append("Consistency comparison issues: threshold not met")

    celltype_ok = True
    if not args.skip_celltype_link:
        celltype_ok, celltype_msg = run_celltype_link(script_dir, args, outdir, language)
        if not celltype_ok:
            if celltype_msg.strip():
                warnings.append(f"Cell-type linking issues: {celltype_msg}")
            else:
                warnings.append("Cell-type linking issues: unknown error")

    consistency_payload = {}
    consistency_path = outdir / "consistency_report.json"
    if consistency_path.exists():
        try:
            consistency_payload = json.loads(consistency_path.read_text(encoding="utf-8"))
        except Exception:
            consistency_payload = {}

    meta = {
        "task": args.task,
        "species": args.species,
        "ontology": ontologies,
        "sig_threshold": args.sig_threshold,
        "mode": args.mode,
        "engine": args.engine,
        "language": language,
        "requested_gsea": requested_gsea,
        "executed_gsea": run_gsea,
        "fallback_note": fallback_note,
        "normalization": norm["stats"],
        "python_ora_backend": python_ora_backend,
        "python_gsea_backend": python_gsea_backend,
        "python_engine_ok": python_ok,
        "r_engine_ok": r_ok,
        "baseline_ok": baseline_ok,
        "consistency_ok": consistency_ok,
        "celltype_link_ok": celltype_ok,
        "celltype_candidates": args.celltype_candidates,
        "run_sources": consistency_payload.get("run_sources", {}),
        "consistency_scores": {
            "ora": (consistency_payload.get("ora", {}) or {}).get("score"),
            "gsea": (consistency_payload.get("gsea", {}) or {}).get("score"),
            "overall_pass": consistency_payload.get("overall_pass"),
        },
        "analysis_brief_ok": False,
        "warnings": warnings,
    }
    write_json(outdir / "run_meta.json", meta)

    analysis_brief_ok, analysis_brief_msg = run_analysis_brief(
        script_dir=script_dir,
        outdir=outdir,
        language=language,
        focus=args.celltype_candidates,
    )
    if not analysis_brief_ok:
        if analysis_brief_msg.strip():
            warnings.append(f"Analysis brief issues: {analysis_brief_msg}")
        else:
            warnings.append("Analysis brief issues: unknown error")
    meta["analysis_brief_ok"] = analysis_brief_ok
    meta["warnings"] = warnings
    write_json(outdir / "run_meta.json", meta)

    run_summary_render(script_dir, outdir, language, fallback_note, warnings)

    if args.benchmark_mode and (not baseline_ok or not consistency_ok):
        return 5
    if no_engine_success:
        return 4

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
