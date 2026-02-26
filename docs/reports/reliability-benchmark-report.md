# GO/GSEA Skill 可靠性基准与实装测试报告（更新）

日期：2026-02-25

## 1. 测试目标

验证 `go-enrichment-reliability` 是否满足以下门禁：

- ORA Top20 重叠率 >= 60%
- GSEA Top20 方向一致率 >= 60%

## 2. 测试框架

- 主入口：`skills/go-enrichment-reliability/scripts/run_go_enrichment.py`
- 基线抓取：`fetch_baselines.py`
- 一致性评估：`compare_results.py`
- 汇总报告：`render_summary.py`

基准输入目录：`benchmarks/data/`

- `ora_human_inflammation_genes.txt`
- `ora_mouse_immune_genes.txt`
- `gsea_broad_tutorial_like_rank.tsv`

## 3. 已执行用例（冒烟）

| 用例 | 命令模式 | 输出目录 | 结果 |
|---|---|---|---|
| ORA human（无基线） | `mode=ora, engine=auto, --skip-baselines` | `runs/smoke_ora3` | 失败（双引擎均不可用） |
| GSEA请求但无score（mouse） | `mode=gsea, engine=auto, --skip-baselines` | `runs/smoke_mouse_inline` | 失败（按规则降级 ORA，但双引擎不可用） |

## 4. 环境实装进展（本次）

- Python 环境已完成（`go-enrichment-reliability` Conda env）：
  - `pandas`, `matplotlib`, `requests`, `gseapy`, `gprofiler-official`, `mygene` 可导入。
- 由于 base conda 求解器版本和环境一致性问题，R/Bioc 依赖（`clusterProfiler`, `AnnotationDbi`, `org.Hs.eg.db`, `org.Mm.eg.db`）尚未完成安装。
- 已清理 Conda 缓存释放磁盘空间（约数 GB），当前可继续安装 R 依赖。

## 5. 功能实装测试（Python 引擎）

### 5.1 文件输入 + ORA/GSEA

- 命令：`engine=python, mode=auto, genes=file, scores=file`
- 输出：`runs/test_python_file_ora_gsea`
- 结果：
  - `ora_results_python.csv`: 2507 行（含表头）
  - `gsea_results_python.csv`: 58 行（含表头）
  - `python_ora_backend=gprofiler`, `python_gsea_backend=gseapy-prerank`
  - 生成完整产物（summary/csv/png/consistency json）

### 5.2 内联输入 + 无 score 自动降级

- 命令：`engine=python, mode=gsea, genes=inline, no scores`
- 输出：`runs/test_python_inline_no_scores`
- 结果：
  - `executed_gsea=false`
  - `summary.md` 包含降级说明（GSEA skipped, downgraded to ORA）
  - `ora_results_python.csv`: 1464 行（含表头）
  - `gsea_results_python.csv`: 1 行（仅表头）

### 5.3 自定义背景集

- 命令：`engine=python, mode=ora, --background benchmarks/data/background_human_test.txt`
- 输出：`runs/test_python_background`
- 结果：任务正常完成并生成全套输出。

### 5.4 阈值必填校验

- 命令：省略 `--sig-threshold`
- 结果：命令行参数校验直接失败，退出码 `2`（符合设计）。

## 6. 基准脚本执行结果（门禁）

执行命令：

```bash
benchmarks/run_reliability_benchmarks.sh
```

终端结果：

- `ora_human`: FAIL (`rc=5`)
- `ora_gsea_human`: FAIL (`rc=5`)
- 汇总：`failed=2`，脚本退出码 `1`

说明：`rc=5` 由 `run_go_enrichment.py --benchmark-mode` 返回，代表基线/一致性门禁未通过。

补充实测（Python 引擎 benchmark mode）：

- 输出：`runs/test_python_benchmark_mode`
- 退出码：`5`
- 主要告警：
  - 缺少 R 依赖导致 `clusterProfiler` 基线不可用
  - 未配置 `BROAD_GSEA_JAR/BROAD_GSEA_GMT`
  - ORA/GSEA 一致性分数未达阈值（60%）

## 7. 失败原因

当前瓶颈已从“Python 与 R 都不可运行”更新为“R 基线未就绪”：

- Python 侧：已就绪。
- R 侧：缺少 `clusterProfiler`, `AnnotationDbi`, `org.Hs.eg.db`, `org.Mm.eg.db`。
- Broad GSEA 严格基线：未配置本地 CLI 所需环境变量。

## 8. 结论

- 框架状态：已就绪（参数契约、输出 schema、对标脚本、门禁逻辑全部可运行）。
- 运行状态：Python 主流程已可执行；R 基线依赖待补齐。
- 风险状态：中（功能可用，但 benchmark 门禁受 R 依赖限制）。

## 9. 下一步（执行级）

1. 安装 R/Bioc 依赖：`clusterProfiler`, `AnnotationDbi`, `org.Hs.eg.db`, `org.Mm.eg.db`。
2. 配置 `BROAD_GSEA_JAR` 与 `BROAD_GSEA_GMT`（若需严格 Broad 基线）。
3. 重新执行 `benchmarks/run_reliability_benchmarks.sh`，并以 `consistency_report.json` 判定是否达标。
