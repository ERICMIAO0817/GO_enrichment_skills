---
name: go-enrichment-reliability
description: Natural-language-first GO/GSEA copilot that auto-runs reliability-checked enrichment (Python-first, R fallback), benchmarks against Enrichr/g:Profiler/clusterProfiler/Broad proxy, and returns interpretable cell-type-focused conclusions with follow-up discussion prompts.
---

# GO Enrichment Reliability

## Positioning

This is a **conversation-first analysis skill**, not a bare script wrapper.

Goal:
- User describes the task in natural language.
- Agent executes reliable GO/GSEA analysis automatically.
- Agent returns interpretable conclusions and supports follow-up discussion.

## Trigger Rules

Use this skill whenever user intent matches any of:
- "Run GO enrichment / GSEA"
- "Analyze this gene list"
- "Focus on beta/alpha/immune/cell-type association"
- "Compare reliability / consistency across tools"
- "Summarize enrichment and discuss biological meaning"

Natural language input is the default entrypoint. CLI-like input is also accepted.

## Input Understanding Contract

Accepted user input forms:
- Free text (Chinese/English)
- Inline list (e.g., `["TP53", "EGFR"]` or `TP53,EGFR`)
- File path for genes
- Optional ranked scores file/list
- Optional CLI-like key-value block

Required semantic fields:
- `genes`
- `sig_threshold`

Defaults if missing:
- `species=human`
- `ontology=ALL`
- `mode=auto`
- `engine=auto`
- `language=auto`

Strict rule:
- If `sig_threshold` is missing, ask for it or block run.

## Agent Execution Policy

When triggered, agent should:
1. Parse user intent and normalize parameters.
2. Run pipeline automatically (do not require user to run commands manually).
3. Prefer Python engine, fallback to R when needed.
4. Run baseline comparisons unless explicitly skipped.
5. Compute consistency metrics.
6. Infer cell-type links with focus candidates if user asked.
7. Generate extended interpretation brief (`analysis_brief.md/json`) with primary and alternative hypotheses.
8. Return biological interpretation, not just raw files.

Only surface commands when user explicitly asks for CLI reproduction.

## Analysis Pipeline

1. Parse genes and optional scores (`run_go_enrichment.py`).
2. Normalize IDs (`normalize_gene_ids.py`).
3. Main run: Python first, R fallback (`run_go_enrichment_r.R`).
4. Fetch baselines (`fetch_baselines.py`).
5. Compare consistency (`compare_results.py`).
6. Infer cell-type association (`infer_celltype_links.py`).
7. Render analysis brief (`render_analysis_brief.py`).
8. Render summary (`render_summary.py`).

## Reliability Definition

- ORA consistency: overlap of run top terms vs baseline consensus top terms.
- GSEA consistency: direction agreement vs baseline consensus direction.
- Baseline comparison is also used for **method-strength ranking** (not only pass/fail gates).
- Default pass thresholds:
  - ORA overlap >= 0.60
  - GSEA direction consistency >= 0.60

`benchmark-mode`:
- Baseline failures and gate failures are hard failures.

Normal mode:
- Baseline issues are warnings; analysis continues.

## Output Contract (for conversation reply)

Audience is fixed to **scientific research**.  
Default depth is fixed to **deep** (not short summary mode).

Agent response should use this narrative template:
1. **Main conclusion (research framing)**: at least 3 full sentences that state the dominant biological axis and what the module most likely represents.
2. **Evidence backbone (term + FDR)**: at least 4-6 key terms grouped by shared biological theme, not a flat term dump.
3. **Cross-method comparison and recommendation**: show primary + baseline methods side-by-side and recommend the strongest evidence source (plus strongest independent corroboration).
4. **Cell-type call and evidence breakdown**: qualitative call + marker hits + keyword hits + ORA/GSEA term hits, with interpretation of what this supports biologically.
5. **What this supports / what it does not**: explicit supported claims and explicit boundaries (consistency limits, missing GSEA directionality, marker coverage limits).
6. **Directly citable paragraph**: one polished paragraph that can be pasted into result slides/manuscript draft.
7. **High-information next analyses**: 2-3 concrete next actions with rationale.

Avoid流水账-style file listing in the main narrative unless user asks for file inventory.

## Stable Files

Each run outputs to `runs/<timestamp>/` (or `--outdir`):

- `summary.md`
- `run_meta.json`
- `normalized_input.tsv`
- `ora_results_python.csv`
- `ora_results_r.csv`
- `gsea_results_python.csv` (if scores)
- `gsea_results_r.csv` (if scores)
- `baseline_enrichr.csv`
- `baseline_gprofiler.csv`
- `baseline_clusterprofiler.csv`
- `baseline_broad_gsea.csv` (GSEA scenario)
- `consistency_report.json`
- `celltype_association.csv`
- `celltype_association.json`
- `celltype_interpretation.md`
- `analysis_brief.md`
- `analysis_brief.json`
- `plots_dotplot.png`
- `plots_barplot.png`

## CLI Reproduction (optional)

If user asks for command-line reproduction, use:

```bash
bash scripts/go-enrich.sh \
  --genes "TP53,EGFR,BRCA1,STAT3,MYC" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL \
  --focus "monocyte,macrophage"
```

Natural-language CLI style is also supported:

```bash
bash scripts/go-enrich.sh \
  --request "用这组基因:['Ins2','Iapp','Pdx1']，sig_threshold=0.05，species=mouse，ontology=ALL，并重点看 beta cell 的关联。"
```

## Troubleshooting

- Missing Python ORA backend: install `gprofiler-official` or `gseapy`.
- R fallback missing deps: install `clusterProfiler`, `AnnotationDbi`, `org.Hs.eg.db` or `org.Mm.eg.db`.
- Broad GSEA baseline: set `BROAD_GSEA_JAR` and `BROAD_GSEA_GMT`.
- Empty enrichment: check species/ID namespace, increase gene count, or relax strictness.
