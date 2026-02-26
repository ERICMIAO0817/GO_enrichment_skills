# Benchmark Protocol

## Goal

Validate run reliability by comparing skill outputs against external baselines.

## Datasets

- ORA benchmark: use curated public GO-related gene sets (human and mouse).
- GSEA benchmark: use Broad tutorial-style ranked datasets.

Store benchmark inputs under repository `benchmarks/data/`.

## Required Baselines

- ORA: Enrichr, g:Profiler, clusterProfiler.
- GSEA: clusterProfiler gseGO, Broad GSEA proxy output.

## Metrics

- ORA metric: Top20 term overlap rate.
- GSEA metric: Top20 pathway direction consistency rate.

## Passing Gates

- ORA overlap >= 0.60
- GSEA direction consistency >= 0.60

## Run Modes

- Normal mode: baseline errors are warnings.
- Benchmark mode (`--benchmark-mode`): any baseline or gate failure fails the run.

## Reproducibility Requirements

- Always keep `runs/<timestamp>/` artifacts.
- Keep `consistency_report.json` for each benchmark run.
- Record dependency lock state using Conda environment file.
