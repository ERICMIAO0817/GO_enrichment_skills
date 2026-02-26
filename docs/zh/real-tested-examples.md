# GO Enrichment Reliability 实测案例（中文）

本页给出 2 个真实运行快照，面向科研场景，采用固定深度模板：
1. 主结论（研究版）
2. 证据主轴（term + FDR）
3. 跨方法证据对比与推荐
4. 细胞类型判断与证据拆解
5. 能说什么 / 不能说什么
6. 可直接引用结果段（科研写作）
7. 下一步高信息增益分析

数据快照日期：2026-02-26。  
说明：以下数值来自本仓库本地运行目录，不是手写示意。

## 案例 1：小鼠 beta 细胞正例（分泌轴清晰，但一致性门禁未通过）

### 输入请求

```text
用 $go-enrichment-reliability 跑这组基因:
['Ins2','Hmgn3','Iapp','Fos','Pdx1','Nnat','Gnas','Ttr','Dynll1','Tuba1a','Calm2','Calr','Rest','Pyy','Hspa5']，
sig_threshold=0.05，species=mouse，ontology=ALL，
并重点看 beta cell 的关联。
```

### 运行快照

1. 输出目录：`runs/beta_mouse_ins2_set_20260226_v4`
2. 输入基因：15，映射成功：15/15（`mygene`）
3. 未提供 scores，GSEA 未执行，自动降级为 ORA
4. 主执行来源：R（Python ORA backend 在本环境缺失）
5. ORA 一致性：`0.50`（阈值 `0.60`），总体 gate：`FAIL`

### 1. 主结论（研究版）

1. 该基因集的富集信号高度集中在“激素与分泌程序”，不是离散噪声条目。显著术语在语义上连续，构成同一生物过程链条，支持“单一主轴驱动”的解释。
2. 细胞类型层面，当前最稳健结论是 `Beta cell` 相关（`associated`），并且不是单点证据，而是 marker 与功能关键词的双路径一致支持。
3. 方法比较层面，`clusterProfiler-R` 给出最强主证据，`g:Profiler` 给出最强独立外部佐证。综合建议采用“主方法 + 最强独立基线”的联合表述，而不是单方法定论。
4. 由于 ORA 共识分数未过阈值且缺失 GSEA 方向信息，该结论更适合作为“稳定功能偏置 + 机制假设起点”，不应被表述为终局状态判定。

### 2. 证据主轴（term + FDR）

1. `hormone transport`（`3.17e-07`）
2. `regulation of protein secretion`（`1.71e-06`）
3. `regulation of hormone secretion`（`1.71e-06`）
4. `hormone secretion`（`3.56e-06`）
5. `regulation of insulin secretion`（`3.56e-06`）
6. `insulin secretion`（`7.36e-06`）
7. `response to endoplasmic reticulum stress`（`2.80e-04`）

这些术语覆盖“分泌调控 -> 激素释放 -> 胰岛素释放 -> ER 应激协同”链条，证据结构完整且具备可解释性。

### 3. 跨方法证据对比与推荐

1. 主方法（`clusterProfiler-R`）：`score=0.82`，`consensus_overlap=1.00`，`focus=0.95`，`sig_terms=281`。
2. `g:Profiler` 基线：`score=0.81`，`consensus_overlap=0.65`，`focus=0.80`，`sig_terms=591`。
3. `clusterProfiler` 基线：`score=0.74`，与主方法同源，独立性较弱。
4. `Enrichr` 基线：`score=0.15`，`consensus_overlap=0.10`，本例证据明显偏弱。
5. 推荐：以 `clusterProfiler-R` 作为主结论来源，以 `g:Profiler` 作为独立佐证来源。

### 4. beta cell 判断与证据拆解

1. 结论：`Beta cell -> associated`，证据类型 `mixed_marker_functional`。
2. marker 命中：`3/10`（`IAPP; INS2; PDX1`）。
3. 功能关键词命中：`5/7`。
4. ORA/GSEA 术语命中：`13/0`（本次无 GSEA）。
5. 未命中 marker：`INS1; ISL1; MAFA; NKX6-1; PCSK1; PCSK2; SLC2A2`。

### 5. 能说什么 / 不能说什么

1. 能说：存在稳定可解释的 beta-biased 分泌功能轴，且有跨方法支持。
2. 能说：该模块可作为机制假设与后续验证起点。
3. 不能说：`ORA consensus=0.50 < 0.60`，跨工具仍有敏感性，不能作为终局定论。
4. 不能说：未提供 ranked scores，GSEA 方向证据缺失，无法判断通路方向驱动。

### 6. 可直接引用结果段（科研写作）

基于该基因集的富集分析，我们观察到结果在功能上集中于激素与分泌程序，并由 `hormone transport`、`regulation of protein/hormone secretion`、`regulation of insulin secretion` 与 `insulin secretion` 等同主题条目共同支持。细胞类型关联分析将该模块判定为 Beta cell associated，其证据来自 marker 命中（3/10）与功能关键词命中（5/7）的双路径一致。在跨方法比较中，clusterProfiler-R 的证据强度最高，而 g:Profiler 提供最强独立佐证。综合来看，该结果更适合解释为稳定的功能偏置与机制假设方向，而非对终末状态的单点定论。

### 7. 下一步高信息增益分析

1. 补充 logFC 或 t-score 后重跑 GSEA，验证分泌轴方向一致性。
2. 在同阈值同流程下做 `beta vs alpha` 并排对照，提升可发表性。
3. 以缺失 marker 为优先清单做 targeted validation，区分“功能偏置”与“稳定谱系状态”。

## 案例 2：人类皮肤/衰老焦点（有相关条目，但细胞类型结论偏弱）

### 输入请求

```text
用 $go-enrichment-reliability 跑这组基因:
['MT-CO1','HLA-B','NEAT1','RPLP1','CYBA','MT2A','COL6A3','FN1','RPS12','S100A10','IGFBP7','RPL32','MT1E','MT-ATP6','TMSB4X','ATP5F1E','RPS28','CDKN1A','MT-CO2','IFITM1','ADM','RPL41','RPS14','RPS15A','RPL37A','DST','COL12A1','FTL','LY6E','FLNA']，
sig_threshold=0.05，species=human，ontology=ALL，
并重点看皮肤细胞与衰老的关联。
```

### 运行快照

1. 输出目录：`runs/skin_senescence_human_20260226_v1`
2. 输入基因：30，映射成功：30/30（`mygene`）
3. 未提供 scores，GSEA 未执行，自动降级为 ORA
4. 主执行来源：R（Python ORA backend 在本环境缺失）
5. ORA 一致性：`0.65`（阈值 `0.60`），总体 gate：`PASS`

### 1. 主结论（研究版）

1. 该基因集的主轴首先是一般性翻译/核糖体与能量代谢程序，而不是强特异的皮肤谱系程序。主显著条目集中在 translation/ribosome 相关过程，提示存在广谱细胞状态信号。
2. 在你指定的“皮肤细胞与衰老”焦点下，确实可见二级支持证据，包括 ECM/胶原、创伤修复和 aging 相关条目，但这些条目强度次于主轴。
3. 细胞类型层面，系统给出 `Fibroblast -> weak_or_uncertain`，说明“有方向性线索，但证据不足以支撑强结论”。
4. 因此该案例适合写成“存在皮肤修复/ECM 与衰老相关信号倾向”，而非“已确认某一皮肤细胞亚群或衰老状态”。

### 2. 证据主轴（term + FDR）

1. `cytoplasmic translation`（`1.10e-07`）
2. `proton transmembrane transport`（`2.32e-03`）
3. `oxidative phosphorylation`（`8.64e-03`）
4. `collagen-containing extracellular matrix`（`5.12e-04`）
5. `extracellular matrix structural constituent`（`5.55e-03`）
6. `wound healing`（`1.45e-02`）
7. `aging`（`3.17e-02`）
8. `positive regulation of fibroblast proliferation`（`3.68e-02`）

### 3. 跨方法证据对比与推荐

1. 主方法（`clusterProfiler-R`）：`score=0.71`，`consensus_overlap=0.95`，`focus=0.15`，`sig_terms=135`。
2. `g:Profiler` 基线：`score=0.66`，`consensus_overlap=0.50`，`focus=0.20`，`sig_terms=165`。
3. `clusterProfiler` 基线：`score=0.63`，与主方法同源，独立性较弱。
4. `Enrichr` 基线：`score=0.28`，`consensus_overlap=0.55`，可作辅助但不是主证据。
5. 推荐：以 `clusterProfiler-R` 作为主证据来源，同时用 `g:Profiler` 提供独立对照，强调主轴与焦点轴的层级差异。

### 4. 皮肤细胞/衰老判断与证据拆解

1. 结论：`Fibroblast -> weak_or_uncertain`，证据类型 `weak_signal`。
2. marker 命中：`0/10`。
3. 功能关键词命中：`1/5`（`wound healing`）。
4. ORA/GSEA 术语命中：`1/0`（本次无 GSEA）。
5. 缺失 marker：`ACTA2; COL1A1; COL1A2; COL3A1; DCN; FAP; LUM; PDGFRB; TAGLN; THY1`。

### 5. 能说什么 / 不能说什么

1. 能说：存在皮肤修复/ECM 与衰老相关条目，但属于二级信号，不是主轴。
2. 能说：方法层面总体一致性通过（ORA=0.65），主结论具备流程稳定性。
3. 不能说：当前证据不足以宣称“强 fibroblast 特异”或“明确 senescence 状态”。
4. 不能说：无 ranked scores，缺少 GSEA 方向证据，无法判断这些过程的上下调驱动方向。

### 6. 可直接引用结果段（科研写作）

在该人类基因集的富集分析中，主显著信号集中于翻译与核糖体相关过程，提示存在广谱细胞状态轴。与此同时，我们观察到 `collagen-containing extracellular matrix`、`extracellular matrix structural constituent`、`wound healing`、`aging` 及 `positive regulation of fibroblast proliferation` 等条目，提示皮肤修复/ECM 与衰老相关功能具有次级支持。跨方法比较显示 clusterProfiler-R 为最强主证据来源，g:Profiler 为最强独立佐证。综合细胞类型证据，当前更稳健的结论是“Fibroblast weak_or_uncertain”，因此结果应被解释为方向性线索与假设起点，而非确定性细胞状态定论。

### 7. 下一步高信息增益分析

1. 引入 ranked scores 运行 GSEA，优先验证 aging/ECM 相关通路是否方向一致。
2. 补充更具皮肤谱系分辨率的 marker 集后重跑 cell-type inference。
3. 增加对照组（非皮肤或非衰老相关基因集）做同流程并排比较，提升结论可发表性。
