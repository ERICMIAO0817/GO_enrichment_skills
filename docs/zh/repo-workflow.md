# 仓库标准流程（2026-02）

这份文档用于统一本仓库的日常使用方式，避免 `runs/` 与实验命名混乱。

## 1. 目录约定

- `skills/go-enrichment-reliability/`：技能本体（脚本 + 参考资料）
- `docs/`：使用说明、流程文档、测试用例
- `benchmarks/`：基准输入数据与基准脚本
- `runs/`：每次分析的输出（可清空、不可手改）
- `logs/`：安装/运行日志（可清空）

## 2. 单次分析标准流程

1. 清理历史运行结果（可选但推荐）

```bash
bash scripts/clean_runs.sh
```

2. 首选对话调用（推荐）

```text
用 $go-enrichment-reliability 跑这组基因:[TP53,EGFR,BRCA1,STAT3,MYC]，
sig_threshold=0.05，species=human，ontology=ALL。
```

3. 如需命令行复现，使用统一入口

```bash
bash scripts/go-enrich.sh \
  --genes "TP53,EGFR,BRCA1,STAT3,MYC" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL
```

4. 优先检查这 4 个文件

```bash
python3 skills/go-enrichment-reliability/scripts/run_go_enrichment.py \
  --genes "TP53,EGFR,BRCA1,STAT3,MYC" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL \
  --mode auto \
  --engine auto
```

- `summary.md`
- `run_meta.json`
- `ora_results_r.csv` 或 `ora_results_python.csv`（看 `run_sources.ora`）
- `consistency_report.json`

## 3. 输出判读规则

- `gsea_results_*.csv` 空：通常是没提供 `scores`，属于预期行为。
- `ora_results_python.csv` 空但 `ora_results_r.csv` 非空：Python 主引擎失败后 R fallback，属于预期行为。
- baseline 大量为空：通常是外部接口失败或手动跳过 baseline。
- `overall_pass=false`：表示一致性门禁未通过，不代表主结果一定无生物学意义。

## 4. 命名规范（建议）

建议每次固定 `--outdir` 命名：

```text
runs/YYYYMMDD_<species>_<focus>_<tag>
```

示例：

- `runs/20260226_mouse_beta_userA`
- `runs/20260226_mouse_alpha_smoke`

## 5. 版本发布前检查

1. 清空运行目录：`bash scripts/clean_runs.sh`
2. 更新发布工件：`bash scripts/publish.sh`
3. 检查文档引用是否仍正确（README/docs）

## 6. 常用命令

预览将要删除的历史 run：

```bash
bash scripts/clean_runs.sh --dry-run
```

只跑 ORA（不带 scores）：

```bash
bash scripts/go-enrich.sh \
  --genes "TNF,IL6,IL1B,NFKB1,RELA,STAT3" \
  --sig-threshold 0.05 \
  --species human \
  --ontology ALL
```
