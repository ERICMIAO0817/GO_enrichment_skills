# GO/GSEA Skill Usage Guide (English)

This guide documents the `go-enrichment-reliability` skill in a publish-ready, conversation-first format.

The skill is designed for scientific research workflows where users want:
1. End-to-end GO ORA/GSEA execution from natural language.
2. Cross-method comparison against baselines.
3. Cell-type-focused biological interpretation with explicit evidence boundaries.

## 0. One-command install (recommended)

From repository root:

```bash
bash install.sh
```

Reinstall/update:

```bash
bash install.sh --force
```

Uninstall:

```bash
bash uninstall.sh
```

Default install location:
- `$CODEX_HOME/skills/go-enrichment-reliability` (if `CODEX_HOME` is set)
- `~/.codex/skills/go-enrichment-reliability` (fallback)

## 1. Trigger from a conversation (recommended)

Use natural language directly in your AI assistant:

```text
Run $go-enrichment-reliability on this gene set:
['TP53','EGFR','BRCA1','STAT3','MYC'],
sig_threshold=0.05, species=human, ontology=ALL,
and focus on monocyte/macrophage association.
```

Recommended fields to include each time:
1. `genes` (inline list or file path)
2. `sig_threshold`
3. `species`
4. `ontology`
5. `scores` (only if you want GSEA directionality)

## 2. CLI reproduction (single entrypoint)

Use the wrapper script instead of calling raw Python scripts:

```bash
bash scripts/go-enrich.sh \
  --genes "TP53,EGFR,BRCA1,STAT3,MYC" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL \
  --focus "monocyte,macrophage"
```

Natural-language CLI input is also supported:

```bash
bash scripts/go-enrich.sh \
  --request "Run this list:['Ins2','Hmgn3','Iapp','Pdx1'], sig_threshold=0.05, species=mouse, ontology=ALL, focus on beta cell association."
```

## 3. Parameter contract

- `--genes`: required, file path or inline list.
- `--sig-threshold`: required.
- `--species`: optional, `human|mouse`, default `human`.
- `--ontology`: optional, `ALL` or `BP,MF,CC`.
- `--scores`: optional ranked file/list for GSEA.
- `--mode`: `auto|ora|gsea`.
- `--engine`: `auto|python|r`.
- `--focus`: optional focus targets for cell-type interpretation.
- `--outdir`: optional output directory.

If scores are missing in `auto` mode, the pipeline degrades to ORA and records that downgrade in `run_meta.json` and `summary.md`.

## 4. Output files to read first

Each run writes to `runs/<timestamp>/` (or `--outdir`):

1. `summary.md`: one-page operational summary.
2. `analysis_brief.md`: deep research narrative.
3. `celltype_association.csv`: structured cell-type evidence.
4. `consistency_report.json`: overlap/consensus metrics.

Then inspect:
- `ora_results_r.csv` / `ora_results_python.csv`
- `gsea_results_*.csv` (if scores were provided and GSEA executed)
- baseline files (`baseline_enrichr.csv`, `baseline_gprofiler.csv`, `baseline_clusterprofiler.csv`, `baseline_broad_gsea.csv`)

## 5. Recommended deep narrative format

For scientific audiences, keep the response structure fixed:

1. Main conclusion (research framing): describe the dominant biological axis and why the result is coherent.
2. Evidence backbone (term + FDR): list 4-8 terms grouped by shared theme, not a flat dump.
3. Cross-method comparison and recommendation: compare primary method and baselines side by side.
4. Cell-type call and evidence breakdown: qualitative call, marker hits, keyword hits, ORA/GSEA term support.
5. What this supports / what it does not: explicit claims and explicit boundaries.
6. Directly citable paragraph: manuscript-ready prose.
7. High-information next analyses: 2-3 concrete follow-up actions with rationale.

## 6. Reliability interpretation

Reliability includes two distinct layers:

1. Consistency gates:
- ORA consensus overlap threshold (default 0.60).
- GSEA direction consistency threshold (default 0.60, only when applicable).

2. Method-strength ranking:
- The framework compares primary and baseline methods and recommends:
  - the strongest evidence source;
  - the strongest independent corroboration source.

A gate failure does not mean “biologically meaningless.” It means cross-tool sensitivity remains and conclusions should be framed as functional bias/hypothesis, not final state claims.

## 7. Common issues

1. Empty `gsea_results_*.csv`:
- Usually no scores were provided, or ranking input is insufficient.

2. Empty Python ORA file but non-empty R ORA file:
- Python backend unavailable, R fallback succeeded.

3. Baseline files sparse or empty:
- External API/network constraints or tool-specific fallback conditions.

4. Weak cell-type call:
- Often caused by low marker coverage and a dominant generic process axis (for example ribosome/translation terms).

## 8. Publish checklist

Before publishing to GitHub:

1. Ensure bilingual docs are linked in `README.md`.
2. Include at least one positive and one boundary-case tested example.
3. Keep examples reproducible with explicit gene lists and run dates.
4. State clearly what can be claimed and what cannot.
5. Avoid script-internal paths in user-facing docs; keep commands at wrapper level (`go-enrich.sh`).
