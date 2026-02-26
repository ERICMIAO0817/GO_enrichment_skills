# GO/GSEA Skill 使用说明（中文）

本文档对应 skill：
`go-enrichment-reliability`

适用场景：
给一组基因（可带 ranked score），直接端到端得到：
1. GO ORA/GSEA 结果
2. 与主流工具的一致性评估
3. 细胞类型关联解释（以定性生物结论为主，不需要手动在网站和大模型之间搬运结果）

该 skill 设计目标是：**用户自然语言提需求，AI 自动执行并继续解释/追问**，而不是让用户手工拼装底层脚本命令。

配套文档：
- 英文使用说明：`docs/en/usage.md`
- 中文实测案例：`docs/zh/real-tested-examples.md`
- 英文实测案例：`docs/en/real-tested-examples.md`

## 0. 一键安装（推荐）

在仓库根目录直接执行：

```bash
bash install.sh
```

覆盖更新：

```bash
bash install.sh --force
```

卸载：

```bash
bash uninstall.sh
```

默认安装位置：
- 若设置 `CODEX_HOME`：`$CODEX_HOME/skills/go-enrichment-reliability`
- 否则：`~/.codex/skills/go-enrichment-reliability`

## 1. 在 Codex 对话里直接调用

你可以在对话里直接说：

```text
用 $go-enrichment-reliability 跑这组基因，sig_threshold=0.05，species=human，ontology=ALL，
并重点看 monocyte, macrophage, t cell, nk cell 的关联。
```

建议你每次都明确给出：
- `sig_threshold`
- `species`
- `genes`（文件路径或内联列表）
- 如果要 GSEA，再给 `scores`

## 2. 在终端运行（推荐用统一入口）

推荐入口（避免直接写 `python .../run_go_enrichment.py`）：
`scripts/go-enrich.sh`

### 2.1 ORA + GSEA + 细胞类型关联

```bash
bash scripts/go-enrich.sh \
  --genes /abs/path/genes.txt \
  --scores /abs/path/ranked_scores.tsv \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL \
  --language zh \
  --focus "monocyte,macrophage,t cell,nk cell"
```

### 2.2 自然语言命令行入口（不手写参数）

```bash
bash scripts/go-enrich.sh \
  --request "用这组基因:['Ins2','Hmgn3','Iapp','Pdx1']，sig_threshold=0.05，species=mouse，ontology=ALL，并重点看 beta cell 的关联。"
```

### 2.3 只有 ORA（无 scores）

```bash
bash scripts/go-enrich.sh \
  --genes "TNF,IL6,IL1B,NFKB1,RELA,STAT3" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL \
  --language zh
```

说明：
- `mode=auto` 且没有 `scores` 时，会自动降级为 ORA（摘要中会注明）。
- `sig-threshold` 缺失会直接报错并阻断运行。

## 2.4 对话式输出建议（科研深度版）

每次运行后，建议 AI 固定按以下结构输出（默认 deep，不走简版）：
- 主结论（研究表述）：至少 3 句，说明主功能轴与模块生物学含义
- 证据主轴（term + FDR）：至少 4-6 条同主题术语，避免平铺术语清单
- 跨方法证据对比与推荐：并排展示主方法与各 baseline，给出“证据最强方法 + 最强独立佐证方法”
- 细胞类型判断：定性结论 + marker/keyword/术语命中 + 生物学解释
- 能说什么 / 不能说什么：显式写出支持点与边界条件
- 可直接引用结果段：一段可贴 PPT/论文草稿的文字
- 下一步高信息增益分析：2-3 条并说明为什么优先

## 3. 关键参数（统一入口 `go-enrich.sh`）

- `--genes`：基因文件路径或内联列表（必填）
- `--scores`：ranked genes（两列：gene + score；可选）
- `--sig-threshold`：显著性阈值（必填）
- `--species`：`human|mouse`，默认 `human`
- `--ontology`：`ALL` 或 `BP,MF,CC`
- `--mode`：`auto|ora|gsea`
- `--engine`：`auto|python|r`
- `--background`：可选背景基因文件
- `--focus`：可选，限制要重点解释的细胞类型（内部映射到 `--celltype-candidates`）
- `--outdir`：输出目录，默认 `runs/<timestamp>`
- `--benchmark-mode`：开启可靠性门禁（用于基准验证）

## 4. 输出文件怎么看

每次运行产出目录：`runs/<timestamp>/`

最常用文件：
- `summary.md`：一页总结（推荐先看）
- `celltype_association.csv`：细胞类型打分表（含证据拆解字段）
- `celltype_association.json`：结构化结果（适合程序化读取）
- `celltype_interpretation.md`：细胞类型深度解读（marker/术语/建议）
- `analysis_brief.md`：扩展分析简报（主假设 + 替代假设 + 不确定性 + 下一步讨论）
- `analysis_brief.json`：扩展分析结构化结果
- `ora_results_python.csv` / `gsea_results_python.csv`：主结果
- `consistency_report.json`：可靠性分数与 pass/fail

`celltype_association.csv` 关键新增字段：
- `qualitative_call`：定性判断（相关/较相关/可能相关/功能提示相关/暂不支持）
- `evidence_basis`：证据类型（marker 主导 / 功能主导 / 混合）
- `ora_term_hits`、`gsea_term_hits`：支持术语命中数
- `top_ora_terms`、`top_gsea_terms`：关键支持术语（含 FDR/NES）
- `missing_markers`：未命中的关键 marker
- `interpretation`：自动生成解释
- `validation_suggestion`：验证建议

定性解读优先阅读：
- `summary.md` 中 “细胞类型关联（端到端）”
- `celltype_interpretation.md` 中 “结论概览” 和 “逐项证据解读”

## 5. 可靠性基准（低占用）

脚本路径：
`benchmarks/run_reliability_benchmarks.sh`

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
GSEAPY_THREADS=1 BROAD_PROXY_THREADS=1 \
GSEAPY_PERM_NUM=20 BROAD_PROXY_PERM_NUM=20 \
BENCH_NICE_LEVEL=15 \
bash benchmarks/run_reliability_benchmarks.sh
```

如果你只想跑 ORA 基准：

```bash
RUN_GSEA_CASE=0 BENCH_NICE_LEVEL=15 \
bash benchmarks/run_reliability_benchmarks.sh
```

## 6. 常见问题

1. 为什么我机器会卡？
   - 用上面的低占用参数（线程=1、`nice`）。
   - 避免同时并行跑多个 benchmark。

2. 为什么 GSEA 可能为空？
   - 可能是 `scores` 太短或分布问题。
   - 先检查 `scores` 至少几十个基因且有正负分布。

3. Broad 基线为什么报配置问题？
   - 需要配置 `BROAD_GSEA_JAR` 和 `BROAD_GSEA_GMT`。
   - 基准脚本默认已提供本地 proxy GMT，可先跑通流程。
