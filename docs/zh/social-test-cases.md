# 社交媒体测试用例清单（对话优先）

目标：让观众看到这不是“脚本库”，而是一个能持续对话的富集分析 skill。

运行前可先清理历史结果：

```bash
bash scripts/clean_runs.sh
```

## A. 对话提示词（直接复制到 AI）

### Case 1: 小鼠 Beta 细胞正例

```text
用 $go-enrichment-reliability 跑这组基因:
['Ins2','Hmgn3','Iapp','Fos','Pdx1','Nnat','Gnas','Ttr','Dynll1','Tuba1a','Calm2','Calr','Rest','Pyy','Hspa5']，
sig_threshold=0.05，species=mouse，ontology=ALL，
并重点看 beta cell 的关联。
```

看点：
- insulin secretion / glucose homeostasis
- beta cell 定性判断与证据基因

### Case 2: 小鼠 Alpha 细胞正例

```text
用 $go-enrichment-reliability 跑这组基因:
['Hmgn3','Slc25a5','Hspa8','Gcg','Actb','Hsp90aa1','Gnas','Iapp','Slc38a5','Ttr','Dynll1','Fos','Calm2','Rps2','Cnksr3','H3f3b']，
sig_threshold=0.05，species=mouse，ontology=ALL，
并重点看 alpha cell 的关联。
```

看点：
- glucagon response / glucagon receptor binding
- alpha cell 关联结论 + 后续验证建议

### Case 3: 物种错配教育向示例

```text
用 $go-enrichment-reliability 跑这组基因:
['Ins2','Hmgn3','Iapp','Fos','Pdx1','Nnat','Gnas','Ttr','Dynll1','Tuba1a','Calm2','Calr','Rest','Pyy','Hspa5']，
sig_threshold=0.05，species=human，ontology=ALL，
并解释为什么结果可能不稳定。
```

看点：
- 同一基因表在不同 species 下稳定性差异
- AI 自动指出 ID/物种匹配问题

### Case 4: Human HGNC 修正对照

```text
用 $go-enrichment-reliability 跑这组基因:
['INS','HMGN3','IAPP','FOS','PDX1','NNAT','GNAS','TTR','DYNLL1','TUBA1A','CALM2','CALR','REST','PYY','HSPA5']，
sig_threshold=0.05，species=human，ontology=ALL，
并和上一个结果做对照解释。
```

看点：
- “错误输入 -> 修正输入”的可解释对照

### Case 5: 触发 GSEA（有 scores）

```text
用 $go-enrichment-reliability 跑：
genes=benchmarks/data/ora_human_inflammation_genes.txt，
scores=benchmarks/data/gsea_broad_tutorial_like_rank.tsv，
sig_threshold=0.05，species=human，ontology=ALL。
请重点解释 GSEA 方向一致性。
```

看点：
- `executed_gsea=true`
- GSEA 一致性和 ORA 一致性的区别

### Case 6: 无 scores 自动降级 ORA

```text
用 $go-enrichment-reliability 跑这组基因:
['TNF','IL6','IL1B','NFKB1','RELA','STAT3']，
sig_threshold=0.05，species=human，ontology=ALL。
并告诉我为什么没有 GSEA 结果。
```

看点：
- AI 主动说明“自动降级到 ORA”

### Case 7: 严格门禁模式

```text
用 $go-enrichment-reliability 以 benchmark_mode 跑：
genes=benchmarks/data/ora_human_inflammation_genes.txt，
scores=benchmarks/data/gsea_broad_tutorial_like_rank.tsv，
sig_threshold=0.05，species=human，ontology=ALL。
并解释 pass/fail 含义。
```

看点：
- baseline/一致性失败时的硬失败语义

### Case 8: 对话式追问（展示可持续探讨）

在任意 Case 1/2 之后追加：

```text
基于你刚才的结果，给我 3 个后续实验/分析建议，
并说明每个建议对应要看的证据文件。
```

看点：
- skill 不只出表，还能继续科学讨论

## B. CLI 复现实验（可选）

如果你要录“命令行版”，用统一入口而不是直接 `python ...`：

```bash
bash scripts/go-enrich.sh --help
```

示例：

```bash
bash scripts/go-enrich.sh \
  --genes "Ins2,Hmgn3,Iapp,Fos,Pdx1,Nnat,Gnas,Ttr,Dynll1,Tuba1a,Calm2,Calr,Rest,Pyy,Hspa5" \
  --sig-threshold 0.05 \
  --species mouse \
  --ontology ALL \
  --focus "beta cell,pancreatic beta cell,islet beta cell" \
  --outdir runs/case_cli_mouse_beta
```

## C. 发帖模板建议

1. 输入提示词截图（展示自然语言入口）
2. `summary.md` 结论段截图
3. `celltype_association.csv` 关键行截图
4. 一句话点评：
   - 例："这个 skill 会自动跑富集 + 基线一致性 + cell-type 解释，并支持继续追问。"

## D. 结果核对清单

- `run_meta.json -> run_sources.ora`
- `run_meta.json -> executed_gsea`
- `consistency_report.json`
- `celltype_association.csv`
