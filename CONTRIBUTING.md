# Contributing to GO Enrichment Skills

This repository hosts local skill implementations (not just links). Contributions should improve reliability, reproducibility, and portability.

## Add or Update a Skill

Place skills under `skills/<skill-name>/` with this minimum structure:

```text
skills/<skill-name>/
├── SKILL.md
└── agents/openai.yaml
```

Optional:
- `scripts/` for deterministic execution
- `references/` for on-demand domain docs
- `assets/` for output templates/resources

## Required Quality Bar

1. `SKILL.md` frontmatter includes only:
   - `name`
   - `description` (must clearly include trigger conditions)
2. Avoid absolute machine paths (for example `/Users/...`) in docs and scripts.
3. Keep reusable scripts executable and tested with at least one representative run.
4. Keep outputs stable and documented when a skill produces files.
5. Add/update docs if behavior or interfaces changed.

## Pull Request Checklist

Before opening a PR:
- [ ] `SKILL.md` trigger description is specific and searchable
- [ ] New/changed scripts run successfully in this repo
- [ ] `bash scripts/publish.sh` has been run after skill metadata changes
- [ ] No generated run artifacts committed (`runs/`, `logs/`, temp files)
- [ ] Docs updated (`README.md`, `docs/`, `specs/` where applicable)

## Scope Note

This project prioritizes skills that are:
- Practically usable end-to-end
- Reproducible with clear input/output contracts
- Supported by validation or benchmark evidence
