# Step 5 — Word 终稿组装（约 45–60 分钟）

## A. 准备（5 分钟）

1. 访达打开：`figures_for_word/`（6 张 PDF 已集中在此）
2. 复制 `农夫原稿.docx` → 另存为 **`GP2_ICA1_终稿.docx`**
3. 并排打开：`2026 GP2 ICA1 guidance_v3 (2).docx`（对照要求）

---

## B. 改标题与小节结构（10 分钟）

在终稿里用 Word「样式」或加粗标题，按下面顺序排列（可先删原稿 `Keywords` 行）：

1. Title  
2. Abstract  
3. **1 Introduction**（合并原 Introduction）  
4. **2 Methods**（2.1–2.8 见下）  
5. **3 Results**（3.1–3.10，guidance 十条）  
6. **4 Discussion**  
7. **References**  
8. **Supplementary materials**（代码说明，不计入 2000 字）

---

## C. 文字：从农夫原稿搬到哪（复制粘贴）

| 终稿位置 | 农夫原稿来源 |
|----------|----------------|
| Abstract | 全文 Abstract |
| 1 Introduction | Introduction 三段 |
| 2.3 Data preprocessing | Methods §2.1 + Table 1 |
| 2.4–2.5 Analysis A/B | Methods §2.2（拆成两段） |
| 2.7 Software | Table 2 |
| 2.8 Limitations | Methods 降采样那段 |
| **3.1 Case study** | **新建 2–3 句** + 说明用 GSE87038 E8.25、memory-light（10 clusters） |
| 3.2 Biological context | Introduction 第 1 段 + Results §3.2 前半 |
| 3.3 Preprocessing | Results §3.2 |
| 3.6–3.7 Analysis A | Results §3.3 + Table 3（A 部分） |
| 3.8 Analysis B | Results §3.4 + Table 3（B 部分） |
| 3.9 Compare A/B | Discussion 前两段精简 |
| 3.10 BioTIP comment | Discussion 局限性 + 方法评价 |
| 4 Discussion | Discussion 全文（删与 3.9/3.10 重复句） |

**3.1 可粘贴英文模板：**

> To validate BioTIP, I applied it to the published mouse gastrulation scRNA-seq dataset (GSE87038, E8.25; Pijuan-Sala et al., 2019) via MouseGastrulationData. After QC and walktrap clustering, BioTIP was run on the top ten largest mesodermal clusters (memory-reduced settings: 1,200 genes, ≤35 cells per cluster). One CTS module (cluster 8, 51 genes) showed a significant ΔIc.shrink (permutation p < 0.01), supporting that BioTIP detects transition-like signals in a known developmental system before application to lung injury data.

---

## D. 插图：文件名 → 插入位置

| 插入顺序 | 文件（在 `figures_for_word/`） | 放在 | 图注示例 |
|----------|-------------------------------|------|----------|
| Fig 1 | `fig_case_mci_panels.pdf` | **3.1** 末尾 | Figure 1. MCI scores per cluster (GSE87038 case study). |
| Fig 2 | `fig_case_ic_permutation.pdf` | **3.1** Fig1 后 | Figure 2. Ic.shrink and permutation null for top CTS (cluster 8). |
| Fig 3 | `fig_umap_celltype.pdf` | **3.3** | Figure 3. UMAP of retained AT2-lineage cells by cell type. |
| Fig 4 | `fig_umap_timepoint.pdf` | **3.3** 可选 | Figure 4. UMAP coloured by time point. |
| Fig 5 | `fig01_cell_composition.pdf` | **3.3** | Figure 5. Cell composition across sampling states. |
| Fig 6 | `fig02_biotip_cts_results.pdf` | **3.6 或 3.8** | Figure 6. ΔIc.shrink for CTS candidates (Analyses A and B). |

Word：**插入 → 图片 → 选 PDF**（或先导出 PNG 再插）。

Table 3 数字用 `02_modules/.../analysisA_simplified_cts_summary.tsv` 与 B 的 tsv；修正 typo：`1.456` 不是 `l.456`。

---

## E. 字数与导出（10 分钟）

1. Word：**审阅 → 字数统计** → 正文 ≤ **2000**（不含 References、图、表、Supplementary）
2. 超了：删 Introduction 重复句、缩短 3.2
3. **文件 → 导出 → PDF** → `GP2_ICA1_提交.pdf`

---

## F. 完成标准（发给 Cursor 验收）

- [ ] `GP2_ICA1_终稿.docx` 存在  
- [ ] Results 有 **3.1–3.10** 小标题（至少 3.1–3.3、3.6–3.10）  
- [ ] 至少 **6 张图**  
- [ ] 已导出 PDF  

回复：**「第5步完成」** + 字数统计截图或「约 XXXX 字」。
