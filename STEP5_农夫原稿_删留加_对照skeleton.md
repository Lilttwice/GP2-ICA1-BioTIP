# Step 5：在 skeleton 指导下组装终稿（删 / 留 / 加）

## 用哪个文件？

| 文件 | 作用 |
|------|------|
| **`GP2_ICA1_report_skeleton.md`** | **目录 + 必交结构**（照着建 Word 标题，不必在 md 里写完） |
| **`农夫原稿.docx`** | **文字仓库**（复制英文段落） |
| **`figures_for_word/*.pdf`** | **唯一要用的图**（11 张） |
| **`GP2_ICA1_终稿.docx`** | **你最后交的 Word**（推荐在这里拼，比纯 md 省事） |

**推荐流程：** 打开 skeleton **当 checklist** → 在 Word 终稿里建同样标题 → 从农夫原稿复制/删/插新图。  
不必先把 skeleton 填满再转 Word（除非你熟悉 Pandoc）。

---

## 一、整份删除（农夫原稿里）

| 删什么 | 原因 |
|--------|------|
| **Keywords** 整行 | guidance 不要求 |
| **3.1 BioTIP case-study rationale** 整节文字 | 只有空话，没有 GSE87038 图；改成真 case study |
| **旧 Figure 1**（composition 那张） | 用 `fig01_cell_composition.pdf` **替换** |
| **旧 Figure 2**（两根柱状 CTS） | 用 `fig02_biotip_cts_results.pdf` **替换**（数字已更新） |
| Results 里 **「Ce l ls reta ined / PBS / D2…」** 乱码行 | Word 导出坏字，删干净 |
| Methods 里 **「reduced analysis scope…tractable」** 整句 | 改成下面「留」里的新参数句（不自我贬低） |
| Abstract 里 **day 6 and 9** 等旧结果 | 改成 **day 5, 8, 9**（见下表） |

---

## 二、整份保留（复制到 skeleton 对应位置）

| 农夫原稿 | 放到终稿（skeleton） |
|----------|---------------------|
| **Title** | 标题（可不改） |
| **Abstract** | Abstract（**改数字**，见第三节） |
| **Introduction** 两段 + aims 段 | **1 Introduction**（合并成 1–3 段即可，不必 1.1/1.2/1.3 太碎） |
| **Methods 2.1** 数据描述 | **2.3 Data preprocessing** |
| **Table 1** | **2.3** 或 **3.3**（保留，数字仍对：10381 cells） |
| **Methods 2.2** BioTIP 流程句 | 拆到 **2.4 / 2.5** |
| **Table 2** 软件 | **2.7 Software** |
| **3.2** 细胞组成生物学叙述 | **3.2 Biological context** + 部分 **3.3** |
| **3.3 Analysis A** 叙述框架 | **3.6–3.7**（**数字必换**） |
| **3.4 Analysis B** 叙述框架 | **3.8**（**数字必换**） |
| **Discussion** 三段 | **4 Discussion** + 拆一点到 **3.9 / 3.10** |
| **References** | References（加 **Pijuan-Sala 2019** case study） |

---

## 三、必须新增（农夫原稿没有或不够）

### 新图（只从 `figures_for_word/` 插，共 11 张）

| Skeleton 节 | 插入文件 | 图号建议 |
|-------------|----------|----------|
| **3.1** BioTIP case study | `fig_case_mci_panels.pdf` | Figure 1 |
| **3.1** | `fig_case_ic_permutation.pdf` | Figure 2 |
| **3.3** Preprocessing | `fig_umap_celltype.pdf` | Figure 3 |
| **3.3** | `fig_umap_timepoint.pdf` | Figure 4（可选） |
| **3.3** | `fig01_cell_composition.pdf` | Figure 5 |
| **3.6** Analysis A | `fig_analysisA_ic_trajectory.pdf` | Figure 6 |
| **3.6** | `fig_analysisA_permutation.pdf` | Figure 7 |
| **3.8** Analysis B | `fig_analysisB_mci.pdf` | Figure 8 |
| **3.8** | `fig_analysisB_ic_trajectory.pdf` | Figure 9 |
| **3.8** | `fig_analysisB_permutation.pdf` | Figure 10 |
| **3.9** Compare A/B | `fig02_biotip_cts_results.pdf` | Figure 11 |

### 新文字（3.1 至少 1 段，可粘贴）

> For validation, BioTIP was applied to mouse gastrulation data (GSE87038 E8.25; MouseGastrulationData). Ten mesodermal clusters were analysed (1,200 genes; ≤35 cells per cluster). A significant CTS was detected in cluster 8 (51 genes; permutation p < 0.01; Figures 1–2), supporting application to lung injury data (GSE141259).

### 新 Methods 句（替换外包「reduced scope」）

> Lung analyses: ≤80 cells per state, 1,500 genes, RTF cutoff 0.1 (bootstrap B = 3), network FDR 0.2, 100 gene permutations. Gastrulation validation: ten largest clusters, 1,200 genes, ≤35 cells per cluster, 50 permutations.

### 新 Table 3（替换农夫原稿整张表）

**Analysis A**

| Candidate | Genes | ΔIc | p | Peak |
|-----------|-------|-----|---|------|
| day 5 | 21 | 1.10 | <0.01 | day 5 |
| day 8 | 23 | 0.44 | 0.01 | day 8 |
| day 9 | 30 | 0.55 | <0.01 | day 9 |

**Analysis B**

| Candidate | Genes | ΔIc | p | Peak |
|-----------|-------|-----|---|------|
| Krt8+ ADI | 11 | 2.61 | 0.01 | Krt8+ ADI |
| AT1 | 11 | 1.18 | 0.04 | Mki67+ Proliferation |
| Mki67+ Proliferation | 22 | 1.00 | <0.01 | Mki67+ Proliferation |

### 新小节标题（guidance 十条，终稿 Results 里要有）

农夫原稿只有 **3.1–3.4**。请改成：

- 3.1 BioTIP successful case study（GSE87038）← **新图 + 新段**
- 3.2 Dataset paper and biological context ← 原 3.2
- 3.3 Data preprocessing ← 原 3.2 后半 + **UMAP + composition 新图**
- 3.4–3.5 Analysis A RTF & networks ← **各 2–3 句**（无图也可）
- 3.6–3.7 Analysis A Ic / interpretation ← 原 3.3 **改数字** + **Fig 6–7**
- 3.8 Analysis B ← 原 3.4 **改数字** + **Fig 8–10**
- 3.9 Compare A and B ← Discussion 前两句 + **Fig 11**
- 3.10 Comment on BioTIP ← Discussion 局限性段

### Supplementary（终稿末尾，不计字数）

> Key scripts: `02_modules/M010.../run_gse87038_case_study_simplified.R`, `M030/preprocess_gse141259.R`, `STEP6_run_lung_AB_figures.R`, `STEP6B_fix_analysisB_figures.R`.

---

## 四、Abstract 必改两句（农夫 day 6/9 → 你 Step6 结果）

**删：** …around days 6 and 9… weaker… day 10… proliferating Mki67+ only…

**改成类似：**

> Time-point analysis identified significant CTS candidates at days 5, 8 and 9. Cell-type analysis highlighted Krt8+ ADI and Mki67+ proliferating states (permutation p ≤ 0.04), consistent with a transitional and proliferative repair programme.

---

## 五、30 分钟操作顺序（照着点）

1. 复制 `农夫原稿.docx` → `GP2_ICA1_终稿.docx`
2. 打开 `GP2_ICA1_report_skeleton.md` **只看左侧标题**
3. 在 Word 里 **删**：Keywords、旧 Fig1/Fig2、3.1 rationale、乱码行、reduced scope 句
4. **改**：Abstract、Table 3、3.3/3.4 里的 day 数字
5. **加标题**：3.1–3.10（没有的空标题也写上，短句填满）
6. **插入 11 张 PDF**（上表顺序）
7. **加** 3.1 段落 + Methods 参数句 + Pijuan-Sala 参考文献
8. 字数统计 → 导出 PDF

---

## 六、和 skeleton 的关系（一句话）

**skeleton = 地图；农夫原稿 = 旧文字；figures_for_word = 新图；终稿 Word = 三样拼起来的成品。**

做完回复：**「Word 拼好了」** 或发字数，我帮你压到 2000 字以内。
