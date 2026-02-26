# GO/Pathway 富集工具全景与选型（2026-02-25）

## 评估目标

为 `go-enrichment-reliability` skill 选择可复现、可自动化、许可风险可控、并可稳定对标的工具组合。

评分维度（1-5，越高越好）：

1. 可复现性
2. 自动化能力
3. 许可/使用约束友好度
4. 更新活跃度
5. 输出可比性（跨工具对比难度）

## 候选评分矩阵

| 工具 | 可复现 | 自动化 | 许可/约束 | 更新活跃 | 可比性 | 备注 |
|---|---:|---:|---:|---:|---:|---|
| GSEApy | 5 | 5 | 4 | 5 | 4 | Python 本地执行，适合流水线与批量任务 |
| clusterProfiler | 5 | 5 | 4 | 4 | 5 | 生信社区经典，R 生态成熟 |
| Enrichr | 3 | 4 | 4 | 4 | 4 | API 简洁，在线服务稳定性受网络影响 |
| g:Profiler | 4 | 5 | 4 | 4 | 5 | API 结构化结果好，跨物种支持好 |
| WebGestalt | 3 | 4 | 4 | 4 | 4 | 2024 API/文档完善，适合补充验证 |
| Reactome | 4 | 4 | 4 | 5 | 4 | 数据发布周期明确，通路分析价值高 |
| DAVID | 3 | 3 | 4 | 3 | 3 | 历史积累强，但自动化与现代接口一般 |
| Broad GSEA | 5 | 3 | 4 | 4 | 5 | 方法学权威，CLI 集成成本偏高 |

## 更新活跃度证据（关键点）

- GSEApy：GitHub releases 显示近期版本仍在发布（例：`v1.1.10`，2025-10-16）。
- clusterProfiler：Bioconductor 页面显示当前 release 线可用（BioC 3.22 分支）。
- Reactome：发布存档显示 `Release 94`（2025-12-10）。
- WebGestalt：官方 2024 手册/API 文档可用。
- DAVID：官方更新页面长期维护但节奏相对慢于前几者。

## 推荐组合

### 组合A（主方案，推荐）

- 执行栈：`Python(g:Profiler/GSEApy) + R(clusterProfiler fallback)`
- 基线对标：`Enrichr + g:Profiler + clusterProfiler + Broad GSEA(代理/CLI)`
- 适用：自动化运行、批量分析、可靠性门禁。
- 优点：可复现性和自动化能力最强；与现有 skill 脚本契合。
- 风险：Broad GSEA 需要额外本地配置。

### 组合B（备选方案）

- 执行栈：`R(clusterProfiler-first)`
- 基线对标：`g:Profiler + WebGestalt + Reactome`
- 适用：R 团队主导的分析流程。
- 优点：统计解释一致性强，结果与文献社区习惯匹配。
- 风险：环境搭建复杂度高于 Python-first。

## 结论

- 第一优先：采用组合A落地当前 skill（已实施）。
- 第二阶段（pathway-enrichment）：以 Reactome + KEGG 为核心通路库，并沿用同一可靠性门禁框架。

## 参考链接

- GSEApy releases: <https://github.com/zqfang/GSEApy/releases>
- clusterProfiler (Bioconductor): <https://bioconductor.org/packages/clusterProfiler/>
- Enrichr API: <https://maayanlab.cloud/Enrichr/help#api>
- g:Profiler APIs: <https://biit.cs.ut.ee/gprofiler/page/apis>
- WebGestalt API manual: <https://www.webgestalt.org/WebGestalt_2024_Manual_API.html>
- Reactome release archive: <https://reactome.org/tag/release>
- DAVID updates: <https://david.ncifcrf.gov/content.jsp?file=updates.html>
- Broad GSEA user guide: <https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm>
