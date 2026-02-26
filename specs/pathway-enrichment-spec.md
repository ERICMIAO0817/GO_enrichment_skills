# Pathway Enrichment Skill 二阶段规格（Reactome + KEGG）

## 1. 目标

在 `go-enrichment-reliability` 成熟后，新增 `pathway-enrichment` skill，支持通路 ORA/GSEA，复用当前 GO skill 的输入契约、输出结构、可靠性门禁与报告框架。

## 2. 范围

- 数据库：Reactome + KEGG（第一优先）
- 模式：ORA + GSEA
- 物种：human + mouse
- 输入方式：文件路径 + 消息内粘贴
- 引擎策略：Python 优先，R fallback

不在二阶段范围：

- 多组学联合富集
- 网络拓扑打分
- 单细胞专用差异通路模块

## 3. 输入契约（沿用 GO skill）

```yaml
task: pathway_enrichment
genes: "<file_path|inline_list>"        # required
sig_threshold: 0.05                       # required each run
species: "human|mouse"                  # default human
databases: "REACTOME,KEGG|ALL"          # default ALL
scores: "<file_path|inline_ranked>"     # optional
mode: "auto|ora|gsea"                   # default auto
engine: "auto|python|r"                 # default auto
background: "<file_path>"               # optional
language: "zh|en|auto"                  # default auto
outdir: "<path>"                        # default runs/<timestamp>
```

行为规则：

1. 缺失 `sig_threshold` 时直接阻断。
2. `mode=auto` 且有 `scores` 时执行 ORA+GSEA。
3. `mode=auto|gsea` 且无 `scores` 时降级 ORA 并写入摘要。

## 4. 输出契约（沿用 GO schema）

目录：`runs/<timestamp>/`

固定文件：

- `summary.md`
- `normalized_input.tsv`
- `ora_results_python.csv`
- `ora_results_r.csv`
- `gsea_results_python.csv`（有 scores 时）
- `gsea_results_r.csv`（有 scores 时）
- `baseline_reactome.csv`
- `baseline_kegg.csv`
- `baseline_clusterprofiler.csv`
- `baseline_broad_gsea.csv`（GSEA 场景）
- `consistency_report.json`
- `plots_dotplot.png`
- `plots_barplot.png`

ORA 字段：
`term_id, term_name, ontology, source, pvalue, padj, gene_hits, gene_ratio, engine`

GSEA 字段：
`term_id, term_name, ontology, source, nes, pvalue, padj, leading_edge, direction, engine`

## 5. 基线与验收

### ORA 基线

- Reactome Analysis Service
- KEGG（合法接口/许可策略下）
- clusterProfiler enrichPathway/enrichKEGG

### GSEA 基线

- clusterProfiler gsePathway/gseKEGG
- Broad GSEA（或兼容代理）

### 指标

- ORA：Top20 通路重叠率 >= 60%
- GSEA：Top20 通路方向一致率 >= 60%

### 失败策略

- 普通分析模式：网络基线失败仅告警。
- benchmark mode：基线失败或阈值不达标即失败。

## 6. 工程实现建议

- 新增脚本：
  - `run_pathway_enrichment.py`
  - `run_pathway_enrichment_r.R`
  - `fetch_pathway_baselines.py`
- 复用脚本：
  - `normalize_gene_ids.py`
  - `compare_results.py`
  - `render_summary.py`
- 可配置项：
  - KEGG 接口策略（REST/API key/许可提醒）
  - Reactome 版本标记

## 7. 兼容与迁移

- 保持与 GO skill 同一参数命名习惯。
- `compare_results.py` 增加数据库维度标签，不破坏现有 JSON 字段。
- 通过环境变量配置 KEGG/Reactome 基线开关，便于离线环境运行。

## 8. 里程碑

1. M1：Pathway ORA Python/R 双引擎与输出 schema 对齐。
2. M2：Pathway GSEA + baseline 对比打通。
3. M3：benchmark gate + 文档与样例数据定版。
