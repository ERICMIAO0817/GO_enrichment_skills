# GO Enrichment Skills

A publish-ready repository for reliability-first GO/GSEA agent skills.

This project packages reproducible enrichment workflows as reusable skills, with:
- Python-first execution and R fallback
- Baseline comparisons (Enrichr, g:Profiler, clusterProfiler, Broad GSEA proxy)
- Consistency gates and structured reports
- End-to-end cell-type association interpretation

Inspired by the curation standards used in [VoltAgent/awesome-agent-skills](https://github.com/VoltAgent/awesome-agent-skills): explicit triggering, portability, security notice, and quality bar.

## Publication Pack (Bilingual)

- Usage (ZH): `docs/zh/usage.md`
- Usage (EN): `docs/en/usage.md`
- Real tested cases (ZH): `docs/zh/real-tested-examples.md`
- Real tested cases (EN): `docs/en/real-tested-examples.md`
- Full docs index: `docs/README.md`

## Table of Contents

- [Publication Pack (Bilingual)](#publication-pack-bilingual)
- [Security Notice](#security-notice)
- [Included Skills](#included-skills)
- [One-Command Install (Codex)](#one-command-install-codex)
- [Installation](#installation)
- [Skills Paths for AI Assistants](#skills-paths-for-ai-assistants)
- [Quick Start](#quick-start)
- [Publish Artifacts](#publish-artifacts)
- [Repository Layout](#repository-layout)
- [Quality Standards](#quality-standards)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## Security Notice

Skills and scripts in this repository are maintained code, not security-audited artifacts.

Before production use:
- Review scripts and external API dependencies
- Validate environment variables and data sources
- Run on non-sensitive datasets first

## Included Skills

<!-- BEGIN_SKILLS_TABLE -->
| Name | Description | Path |
|---|---|---|
| `go-enrichment-reliability` | Natural-language-first GO/GSEA copilot that auto-runs reliability-checked enrichment (Python-first, R fallback), benchmarks against Enrichr/g:Profiler/clusterProfiler/Broad proxy, and returns interpretable cell-type-focused conclusions with follow-up discussion prompts. | `skills/go-enrichment-reliability/` |
<!-- END_SKILLS_TABLE -->

## One-Command Install (Codex)

Install directly into Codex skills directory:

```bash
bash install.sh
```

Update/reinstall:

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

## Installation

### Codex (manual copy/symlink)

Copy a skill folder into your Codex skills path:

```bash
cp -R skills/go-enrichment-reliability ~/.codex/skills/
```

### Claude Code (plugin marketplace)

```bash
/plugin marketplace add <your-github-user>/GO_enrichment_skills
/plugin install go-enrichment-reliability@go-enrichment-skills
```

### Cursor (plugin manifest)

Use `.cursor-plugin/plugin.json` and the `skills/` directory in this repository.

### skills.sh ecosystem

```bash
npx skills add <your-github-user>/GO_enrichment_skills
```

## Skills Paths for AI Assistants

| Tool | Project Path |
|---|---|
| This repository | `skills/` |
| Codex | `.codex/skills/` |
| Claude Code | `.claude/skills/` |
| Cursor | `.cursor/skills/` |
| Gemini CLI | `.gemini/skills/` |
| GitHub Copilot | `.github/skills/` |

## Quick Start

0. Clean old run outputs (optional but recommended):

```bash
bash scripts/clean_runs.sh
```

1. Use natural language in your AI assistant:

```text
用 $go-enrichment-reliability 跑这组基因:[TP53,EGFR,BRCA1,STAT3,MYC]，
sig_threshold=0.05，species=human，ontology=ALL，
并重点看 monocyte/macrophage 的关联。
```

2. Or run via unified CLI wrapper (not raw python path):

```bash
bash scripts/go-enrich.sh \
  --genes "TP53,EGFR,BRCA1,STAT3,MYC" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL \
  --focus "monocyte,macrophage"
```

3. Or pass a natural-language request directly to CLI:

```bash
bash scripts/go-enrich.sh \
  --request "用这组基因:['Ins2','Iapp','Pdx1']，sig_threshold=0.05，species=mouse，ontology=ALL，并重点看 beta cell 的关联。"
```

4. Run benchmark checks:

```bash
bash benchmarks/run_reliability_benchmarks.sh
```

5. Read outputs from `runs/<timestamp>/`:
- `summary.md`
- `consistency_report.json`
- `celltype_association.csv`
- `celltype_interpretation.md`
- `analysis_brief.md`

## Publish Artifacts

Generate catalog and plugin metadata from `SKILL.md` files:

```bash
bash scripts/publish.sh
```

This regenerates:
- `skills/README.md`
- `agents/AGENTS.md`
- `.claude-plugin/marketplace.json`
- `.claude-plugin/plugin.json`
- `.cursor-plugin/plugin.json`

## Repository Layout

```text
.
├── AGENTS.md
├── CLAUDE.md
├── install.sh
├── uninstall.sh
├── .claude-plugin/
├── .cursor-plugin/
├── agents/
├── scripts/
├── skills/
│   └── go-enrichment-reliability/
│       ├── SKILL.md
│       ├── agents/openai.yaml
│       ├── scripts/
│       └── references/
├── benchmarks/
├── docs/
│   ├── en/
│   ├── zh/
│   ├── reports/
│   └── README.md
├── specs/
└── runs/                  # generated artifacts (ignored in git)
```

## Quality Standards

| Area | Guideline |
|---|---|
| Description quality | State what the skill does and when to trigger it in `SKILL.md` frontmatter. |
| Progressive disclosure | Keep core guidance in `SKILL.md`; load `references/` only when needed. |
| Portability | Avoid machine-specific absolute paths in docs and scripts. |
| Reproducibility | Preserve stable output contracts and explicit thresholds. |
| Reliability | Report overlap/consistency scores and pass/fail gates. |

## Documentation

- Docs index: `docs/README.md`
- Usage (ZH): `docs/zh/usage.md`
- Usage (EN): `docs/en/usage.md`
- Real tested examples (ZH): `docs/zh/real-tested-examples.md`
- Real tested examples (EN): `docs/en/real-tested-examples.md`
- Repo workflow standard (ZH): `docs/zh/repo-workflow.md`
- Social-media test cases (ZH): `docs/zh/social-test-cases.md`
- Tool landscape and selection notes: `docs/reports/tool-landscape.md`
- Benchmark report: `docs/reports/reliability-benchmark-report.md`
- Skills ecosystem benchmark: `docs/reports/skills-ecosystem-benchmark.md`
- Pathway phase-2 spec: `specs/pathway-enrichment-spec.md`

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution rules and quality checklist.

## License

MIT. See [LICENSE](LICENSE).
