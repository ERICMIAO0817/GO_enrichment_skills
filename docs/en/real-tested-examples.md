# GO Enrichment Reliability: Real Tested Examples (English)

This page provides 2 real run snapshots for scientific communication using a fixed deep template:
1. Main conclusion (research framing)
2. Evidence backbone (term + FDR)
3. Cross-method comparison and recommendation
4. Cell-type call and evidence breakdown
5. What this supports / what it does not
6. Directly citable result paragraph
7. High-information next analyses

Snapshot date: 2026-02-26.  
Note: all values below come from real local run artifacts in this repository.

## Case 1: Mouse beta-cell positive case (clear secretion axis, gate not passed)

### Input request

```text
Run $go-enrichment-reliability on this gene set:
['Ins2','Hmgn3','Iapp','Fos','Pdx1','Nnat','Gnas','Ttr','Dynll1','Tuba1a','Calm2','Calr','Rest','Pyy','Hspa5'],
sig_threshold=0.05, species=mouse, ontology=ALL,
and focus on beta-cell association.
```

### Run snapshot

1. Output directory: `runs/beta_mouse_ins2_set_20260226_v4`
2. Input genes: 15, mapped: 15/15 (`mygene`)
3. No ranked scores were provided, so GSEA was skipped and downgraded to ORA
4. Primary execution source: R (Python ORA backend unavailable in this environment)
5. ORA consistency: `0.50` (threshold `0.60`), overall gate: `FAIL`

### 1. Main conclusion (research framing)

1. The enrichment signal is not a noisy term collection; it is concentrated on a coherent hormone/secretion axis. Top terms are semantically continuous and support a single biological program interpretation.
2. At the cell-type layer, the most stable interpretation is `Beta cell` `associated`, supported by dual evidence channels rather than a single indicator.
3. In cross-method comparison, `clusterProfiler-R` provides the strongest primary evidence and `g:Profiler` provides the strongest independent corroboration.
4. Because ORA consensus is below threshold and no GSEA directionality is available, this should be framed as a stable functional bias and hypothesis direction, not a terminal state claim.

### 2. Evidence backbone (term + FDR)

1. `hormone transport` (`3.17e-07`)
2. `regulation of protein secretion` (`1.71e-06`)
3. `regulation of hormone secretion` (`1.71e-06`)
4. `hormone secretion` (`3.56e-06`)
5. `regulation of insulin secretion` (`3.56e-06`)
6. `insulin secretion` (`7.36e-06`)
7. `response to endoplasmic reticulum stress` (`2.80e-04`)

These terms cover a coherent chain: secretion regulation -> hormone release -> insulin release -> ER stress coupling.

### 3. Cross-method comparison and recommendation

1. Primary method (`clusterProfiler-R`): `score=0.82`, `consensus_overlap=1.00`, `focus=0.95`, `sig_terms=281`.
2. `g:Profiler` baseline: `score=0.81`, `consensus_overlap=0.65`, `focus=0.80`, `sig_terms=591`.
3. `clusterProfiler` baseline: `score=0.74`, same family as primary method, weaker independence.
4. `Enrichr` baseline: `score=0.15`, `consensus_overlap=0.10`, clearly weaker in this case.
5. Recommendation: report `clusterProfiler-R` as primary evidence, and `g:Profiler` as independent corroboration.

### 4. Beta-cell call and evidence breakdown

1. Call: `Beta cell -> associated`, evidence basis `mixed_marker_functional`.
2. Marker hits: `3/10` (`IAPP; INS2; PDX1`).
3. Functional keyword hits: `5/7`.
4. ORA/GSEA term support: `13/0` (no GSEA run in this case).
5. Missing markers: `INS1; ISL1; MAFA; NKX6-1; PCSK1; PCSK2; SLC2A2`.

### 5. What this supports / what it does not

1. Supports: a stable and interpretable beta-biased secretion axis with multi-method support.
2. Supports: using this module as a mechanistic hypothesis starting point for validation.
3. Does not support: a final cross-tool robust claim, because `ORA consensus=0.50 < 0.60`.
4. Does not support: pathway direction claims, because ranked scores were not provided and GSEA was not executed.

### 6. Directly citable result paragraph

Enrichment analysis of this gene set showed a coherent hormone/secretion program supported by convergent terms including `hormone transport`, `regulation of protein/hormone secretion`, `regulation of insulin secretion`, and `insulin secretion`. Cell-type association classified the module as Beta cell associated, with dual-path support from marker hits (3/10) and functional keyword hits (5/7). In cross-method comparison, clusterProfiler-R provided the strongest primary evidence, while g:Profiler provided the strongest independent corroboration. Overall, this result is best interpreted as a stable functional bias and mechanistic hypothesis direction rather than a terminal-state single-point conclusion.

### 7. High-information next analyses

1. Add logFC or t-scores and rerun GSEA to test directional agreement of the secretion axis.
2. Run matched `beta vs alpha` comparisons under identical thresholds to strengthen publishability.
3. Prioritize targeted validation on missing markers to separate functional bias from stable lineage state.

## Case 2: Human skin/senescence focus (relevant terms present, but weak cell-type call)

### Input request

```text
Run $go-enrichment-reliability on this gene set:
['MT-CO1','HLA-B','NEAT1','RPLP1','CYBA','MT2A','COL6A3','FN1','RPS12','S100A10','IGFBP7','RPL32','MT1E','MT-ATP6','TMSB4X','ATP5F1E','RPS28','CDKN1A','MT-CO2','IFITM1','ADM','RPL41','RPS14','RPS15A','RPL37A','DST','COL12A1','FTL','LY6E','FLNA'],
sig_threshold=0.05, species=human, ontology=ALL,
and focus on skin-cell and senescence association.
```

### Run snapshot

1. Output directory: `runs/skin_senescence_human_20260226_v1`
2. Input genes: 30, mapped: 30/30 (`mygene`)
3. No ranked scores were provided, so GSEA was skipped and downgraded to ORA
4. Primary execution source: R (Python ORA backend unavailable in this environment)
5. ORA consistency: `0.65` (threshold `0.60`), overall gate: `PASS`

### 1. Main conclusion (research framing)

1. The dominant signal is a broad translation/ribosome plus energy-metabolism program rather than a highly specific skin-lineage program.
2. Under the requested skin/senescence focus, secondary evidence does appear, including ECM/collagen, wound healing, and aging-related terms.
3. The cell-type inference output is `Fibroblast -> weak_or_uncertain`, indicating directional clues but insufficient strength for a high-confidence identity call.
4. This case should be framed as skin-repair/ECM and aging-related tendency, not confirmed fibroblast lineage or definitive senescence state.

### 2. Evidence backbone (term + FDR)

1. `cytoplasmic translation` (`1.10e-07`)
2. `proton transmembrane transport` (`2.32e-03`)
3. `oxidative phosphorylation` (`8.64e-03`)
4. `collagen-containing extracellular matrix` (`5.12e-04`)
5. `extracellular matrix structural constituent` (`5.55e-03`)
6. `wound healing` (`1.45e-02`)
7. `aging` (`3.17e-02`)
8. `positive regulation of fibroblast proliferation` (`3.68e-02`)

### 3. Cross-method comparison and recommendation

1. Primary method (`clusterProfiler-R`): `score=0.71`, `consensus_overlap=0.95`, `focus=0.15`, `sig_terms=135`.
2. `g:Profiler` baseline: `score=0.66`, `consensus_overlap=0.50`, `focus=0.20`, `sig_terms=165`.
3. `clusterProfiler` baseline: `score=0.63`, same-family method, weaker independence.
4. `Enrichr` baseline: `score=0.28`, `consensus_overlap=0.55`, useful as auxiliary context only.
5. Recommendation: use `clusterProfiler-R` as primary source and `g:Profiler` as independent support, while explicitly separating dominant axis from focus-axis evidence.

### 4. Skin/senescence call and evidence breakdown

1. Call: `Fibroblast -> weak_or_uncertain`, evidence basis `weak_signal`.
2. Marker hits: `0/10`.
3. Functional keyword hits: `1/5` (`wound healing`).
4. ORA/GSEA term support: `1/0` (no GSEA run in this case).
5. Missing markers: `ACTA2; COL1A1; COL1A2; COL3A1; DCN; FAP; LUM; PDGFRB; TAGLN; THY1`.

### 5. What this supports / what it does not

1. Supports: skin-repair/ECM and aging-related terms are present as secondary evidence.
2. Supports: process-level stability is acceptable (`ORA=0.65`, gate passed).
3. Does not support: a strong fibroblast-specific or definitive senescence-state claim.
4. Does not support: pathway directionality, because ranked scores were not provided and GSEA was not executed.

### 6. Directly citable result paragraph

In this human gene-set enrichment run, the strongest signal was a broad translation/ribosome-centered process axis, indicating a general cellular state component. In parallel, secondary support was observed for skin-relevant and aging-related biology, including `collagen-containing extracellular matrix`, `extracellular matrix structural constituent`, `wound healing`, `aging`, and `positive regulation of fibroblast proliferation`. Cross-method comparison identified clusterProfiler-R as the strongest primary evidence source and g:Profiler as the strongest independent corroboration source. Considering the cell-type evidence profile, the most defensible interpretation is Fibroblast weak_or_uncertain, which supports a directional hypothesis rather than a definitive cell-state assignment.

### 7. High-information next analyses

1. Add ranked scores and run GSEA to verify direction-consistent aging/ECM pathway behavior.
2. Expand marker panels with higher skin-lineage resolution and rerun cell-type inference.
3. Add matched non-skin/non-senescence controls under the same pipeline to improve publication-grade contrast.
