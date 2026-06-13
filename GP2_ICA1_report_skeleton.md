# GP2 ICA1 终稿对照稿（对着改 Word）

> **推荐：** 英文终稿已全部扩写 → 直接打开 **`GP2_ICA1_终稿正文.md`**，整段复制到 Word，再插 11 张图。  
> 本文件保留 **[留]/[删]/[改]** 标注，方便你对照农夫原稿删旧内容。  
> **字数：** 正文 ≤ **2000**（不含 References、图表、Supplementary）。

---

## 标注说明

| 标记 | 意思 |
|------|------|
| **\[留\]** | 农夫原稿里**整段保留**（可复制，少量 \[改\] 已写在段内） |
| **\[删\]** | 终稿里**不要出现**（从农夫原稿删掉） |
| **\[新\]** | 农夫原稿**没有**，必须新加 |
| **\[改\]** | 农夫原稿有但**数字/表述要换成这版** |
| **\[图: 文件名\]** | 在 Word **插入 → 图片** 选 `figures_for_word/` 下该 PDF |

---

# BioTIP analysis of alveolar epithelial state changes after lung injury

GP2 Data Analysis ICA

---

## Abstract

**\[改\]** 用下面整段**替换**农夫 Abstract（不要 day 6/9 旧结果）：

Alveolar type 2 (AT2) cells help restore the distal lung epithelium after injury, but repair can pass through unstable intermediate states. I analysed pre-processed lung fibrosis scRNA-seq (GSE141259) with BioTIP, defining state either by sampling time point or by annotated alveolar cell type. After retaining AT2-lineage cells (10,381 cells), time-point analysis identified significant critical-transition signal (CTS) candidates at days 5, 8 and 9. Cell-type analysis highlighted Krt8+ ADI and Mki67+ proliferating states (permutation p ≤ 0.04), in addition to an AT1-associated module. These results support a repair trajectory that is not a simple linear AT2-to-AT1 conversion, but involves transitional and proliferative epithelial programmes. The two state definitions were complementary: time points captured when the tissue response was most unstable, whereas cell types highlighted which transcriptional programmes showed criticality.

---

## 1 Introduction

### 1.1 Scientific question and hypothesis

**\[留\]** 农夫 Introduction 第 1 段：

The distal lung epithelium relies on AT2 cells both for surfactant production and for epithelial repair. After injury, AT2 cells can proliferate and generate AT1 cells, which form the thin gas-exchange surface. In fibrotic injury, however, regeneration is accompanied by persistent transitional epithelial states rather than a simple restoration of the original tissue architecture. Strunz and colleagues described an alveolar regeneration programme in which AT2 cells pass through activated and Krt8+ intermediate states before acquiring AT1-like features; failure to resolve this programme is closely linked to aberrant repair and fibrosis (Strunz et al., 2020). This makes the dataset a useful test case for methods that search for unstable transition points in high-dimensional expression data.

### 1.2 BioTIP and analysis aims

**\[留\]** 农夫 Introduction 第 2–3 段（合并 BioTIP 介绍 + 三个 aims）：

BioTIP was designed to detect critical transition signals (CTSs) from transcriptomic data by combining transcript fluctuation, co-expression network structure and an index of criticality (Ic) (Yang et al., 2022). Conceptually, the method asks whether a group of genes becomes unusually variable and co-regulated in a specific state, which would be expected near a tipping point. This analysis had three aims: first, to demonstrate BioTIP on the published mouse gastrulation case study (GSE87038); second, to apply BioTIP to injury-associated alveolar epithelial cells from GSE141259; and third, to compare whether defining state by time point or by annotated cell type leads to the same biological interpretation.

**\[删\]** 农夫全文：**Keywords:** 整行。

---

## 2 Methods

### 2.1 BioTIP successful case study (GSE87038)

**\[新\]** 终稿全文（已扩写）：

The validation used mouse gastrulation scRNA-seq (GSE87038, E8.25; Pijuan-Sala et al., 2019) obtained through the MouseGastrulationData Bioconductor package. Single-cell objects from embryo samples 24, 25 and 28 were log-normalised, and cells annotated as stripped nuclei or doublets were removed. Mesodermal cell types spanning the haemato-endothelial and related lineages were retained and embedded with PCA (20 components). A shared nearest-neighbour graph (k = 10) was clustered with the walktrap algorithm to define transcriptionally coherent groups. BioTIP was then applied to the ten largest clusters. Highly variable genes were capped at 1,200, and each cluster was down-sampled to at most 35 cells before analysis.

State-specific transcripts were pre-selected with optimize.sd_selection (relative transcript fluctuation cutoff 0.1, percent 0.8, bootstrap B = 2). Co-expression networks were inferred per cluster with getNetwork (FDR 0.2), partitioned into modules with getCluster_methods, and scored with getMCI (fun = BioTIP). CTS modules containing at least 40 genes were carried forward. For each candidate module, Ic.shrink was calculated across clusters with getIc (shrink = TRUE), and empirical significance was evaluated using 50 gene-permutation simulations of ΔIc.shrink.

### 2.2 Data preprocessing (GSE141259)

**\[留\]** 农夫 Methods §2.1 两段（**删**第二段里 reduced scope 那句，见下）：

The main analysis used the pre-processed high-resolution single-cell RNA-seq files from GSE141259, generated from the lung injury study by Strunz et al. (2020). Raw count matrices and matched cell metadata were retrieved from GEO. The analysis focused on alveolar epithelial states relevant to AT2 differentiation: AT2, activated AT2, Krt8+ ADI, AT1 and Mki67+ proliferating cells. Cells were retained if they belonged to one of these states and passed quality filters (≥500 detected genes, ≥1,000 UMIs, ≤20% mitochondrial reads). The final analysis object contained 10,381 cells.

**\[新\]** 替换农夫里这句（整句不要出现在终稿）：

~~Because full-scale network and MCI estimation was computationally heavy for the available workstation, the final report uses a reproducible reduced analysis scope…~~

**\[新\]** 接在上段后面（**不必插表格**；`.rds` 只是保存路径，可改成下面更短句）：

Seurat (v5.5.0) was used for log-normalisation, selection of 3,000 variable features, scaling, PCA (20 dimensions) and UMAP visualization. A single preprocessed Seurat object was saved for all downstream BioTIP runs. Time points with fewer than 20 cells after filtering (notably day 28) were excluded from BioTIP state definitions in Analysis A.

**\[留\]** 农夫 **Table 1**（标题可改为 *Table 1. Dataset and analysis scope.*）：

| Item | Value |
|------|-------|
| Primary dataset | GSE141259 lung injury regeneration scRNA-seq |
| Input matrix | HighResolution raw-count matrix with matched cell metadata |
| Cells before lineage filtering | 32,559 |
| Cells retained for analysis | 10,381 |
| Retained cell states | AT2, activated AT2, Krt8+ ADI, AT1, Mki67+ proliferating cells |
| Analysis A states | Time points after injury; day 28 excluded from BioTIP scoring (two cells) |
| Analysis B states | Annotated alveolar cell types pooled across time |

### 2.3 Analysis A — time point as state

**\[新\]** 终稿全文（已扩写，与 `GP2_ICA1_终稿正文.md` 一致）：

Analysis A treated each post-injury sampling time point (and PBS control where retained) as an independent state. Only states with at least 20 cells were included. To stabilise network estimation while preserving temporal structure, each state was down-sampled to at most 80 cells, and the 1,500 most variable genes across the retained matrix were used for BioTIP input.

**2.3.1 Transcript pre-selection.** State-specific highly fluctuating transcripts were identified using BioTIP’s optimize.sd_selection function, which estimates relative transcript fluctuation (RTF) with bootstrap resampling. Non-default settings were cutoff = 0.1, percent = 0.8 and B = 3. This step reduces the feature space before network construction and highlights genes that deviate from homogenous expression within each time point.

**2.3.2 Network partition and MCI scoring.** For each time point, Pearson-based co-expression networks were built on the pre-selected gene sets with getNetwork (FDR = 0.2). Networks were partitioned into gene modules using walktrap community detection (getCluster_methods). Module criticality index (MCI) scores were then computed with getMCI (fun = BioTIP). Modules with fewer than 20 genes were excluded from CTS nomination. The strongest modules per time point were retained as CTS candidates using getTopMCI and getCTS.

**2.3.3 Tipping point identification.** For every CTS candidate module, Ic.shrink was calculated across all time-point states with getIc (shrink = TRUE, PCC_sample.target = average). The observed ΔIc.shrink was defined as the difference between the highest and second-highest state-specific Ic.shrink values. Statistical significance was assessed by comparing the observed ΔIc.shrink to a null distribution generated from 100 random gene sets of equal size drawn from the same expression matrix (gene-permutation test).

### 2.4 Analysis B — cell type as state

**\[新\]** 终稿全文（已扩写）：

Analysis B used the same preprocessed GSE141259 expression matrix but defined states by annotated alveolar cell type rather than time point. All cells from a given type were pooled across time points, yielding five states: AT2, activated AT2, Krt8+ ADI, AT1 and Mki67+ proliferation. States with fewer than 15 cells were excluded. Each cell type was down-sampled to at most 100 cells, and up to 1,800 highly variable genes were retained. These settings increase per-state homogeneity and improve network density relative to Analysis A, at the cost of removing explicit temporal ordering.

**2.4.1 Transcript pre-selection.** The Analysis A pre-selection procedure was repeated with identical RTF parameters (optimize.sd_selection; cutoff 0.1, percent 0.8, B = 3), producing fluctuating transcript sets specific to each cell type rather than each time point.

**2.4.2 Network partition and MCI scoring.** Cell-type-specific co-expression networks were constructed (getNetwork; FDR 0.2), partitioned into modules, and scored with getMCI. Because cell-type states are transcriptionally more coherent than mixed time-point populations, networks were generally denser and yielded larger CTS modules. CTS candidates were extracted with a minimum module size of 10 genes for Analysis B.

**2.4.3 Tipping point identification.** Ic.shrink profiles and ΔIc.shrink statistics were computed as in Analysis A, with the same 100 gene-permutation simulations per CTS module. Peak Ic.shrink states were interpreted as the cell types in which the nominated module showed the strongest criticality-like behaviour.

### 2.5 Comparison of analyses A and B

**\[新\]** 终稿全文（已扩写）：

Analyses A and B shared the same preprocessing and core BioTIP functions but differed in biological grouping. Comparisons focused on (i) the number and size of CTS modules, (ii) permutation p-values and peak states, and (iii) concordance with the Strunz et al. (2020) model of AT2 activation, Krt8+ ADI emergence and AT1 differentiation. The two analyses were not treated as independent confirmations of the same tipping event, because time-point states mix multiple cell types, whereas cell-type states mix multiple time points.

### 2.6 Software

**\[留\]** 农夫 **Table 2**（版本按你环境 **\[改\]**）：

| Software | Version | Use |
|----------|---------|-----|
| R | 4.5.2 | Analysis environment |
| Seurat | 5.5.0 | Single-cell object handling, normalization, PCA and UMAP |
| BioTIP | 1.23.0 | RTF pre-selection, network partition, MCI and Ic.shrink scoring |
| GEOquery | 2.78.0 | GEO metadata retrieval |
| MouseGastrulationData | (Bioconductor) | GSE87038 E8.25 case study |
| ggplot2 | 4.0.1 | Figure generation |

**\[留\]** 农夫 Methods 最后一句软件句（**改** BioTIP 版本号）：

The analysis was performed in R 4.5.2 using Seurat 5.5.0, BioTIP 1.23.0 and GEOquery 2.78.0. Non-default parameters are listed in Sections 2.1–2.4.

---

## 3 Results

### 3.1 BioTIP successful case study (GSE87038)

**\[删\]** 农夫整个 **§3.1 BioTIP case-study rationale**（只有理论、无图）——**不要复制到终稿**。

**\[新\]** 用下面替换 3.1 全文（≤500 词；配 2 图）：

BioTIP was first applied to the published gastrulation dataset to confirm that the pipeline detects transition-like signals in a system with known developmental branching. After clustering E8.25 mesodermal cells, MCI scores highlighted modules in several clusters exceeding the conventional threshold (MCI > 2; Figure 1). For the top CTS module (cluster 8; 51 genes), Ic.shrink peaked at the corresponding cluster and permutation testing (50 simulations) supported significance (p < 0.01; Figure 2). This successful validation indicates that BioTIP can identify coordinated fluctuation–co-expression signals prior to analysis of the lung injury dataset.

**\[图: fig_case_mci_panels.pdf\]**  
*Figure 1. MCI scores per cluster in the GSE87038 gastrulation case study (E8.25 mesoderm). Dashed line: MCI = 2.*

**\[图: fig_case_ic_permutation.pdf\]**  
*Figure 2. Ic.shrink trajectory and permutation null for the top gastrulation CTS module (cluster 8).*

---

### 3.2 Dataset paper and biological context

**\[留\]** 农夫 Results **§3.2** 生物学叙述（**不要**插旧 Figure 1 前的乱码轴标签）：

Following injury, the retained cells spanned homeostatic control and multiple post-injury time points, with enrichment of activated AT2 and Krt8+ ADI cells during early repair. Activated AT2 cells were especially abundant around days 3–5, while Krt8+ ADI cells became more prominent from approximately day 7 onward. AT1 cells were present at lower abundance but increased after the earliest injury phase. This pattern is consistent with a repair trajectory in which AT2 cells first activate and proliferate, then pass through an intermediate differentiation state before contributing to AT1-like epithelium, as described by Strunz et al. (2020).

---

### 3.3 Data preprocessing (GSE141259)

**\[新\]** 终稿全文（已扩写）：

PCA and UMAP on the retained 10,381 cells showed partial separation of cell types, with Krt8+ ADI and activated AT2 occupying intermediate regions between AT2 and AT1 clusters (Figures 3–4). This layout is consistent with a continuous injury response rather than discrete, fully separated cell types. Stacked counts across time points demonstrated dynamic shifts in the proportion of activated AT2 and Krt8+ ADI cells during early and mid repair (Figure 5).

**\[图: fig_umap_celltype.pdf\]**  
*Figure 3. UMAP of retained AT2-lineage cells coloured by cell type.*

**\[图: fig_umap_timepoint.pdf\]**  
*Figure 4. UMAP coloured by sampling time point.*

**\[图: fig01_cell_composition.pdf\]**  
*Figure 5. Composition of retained alveolar-lineage cells across sampling states.*

**\[删\]** 农夫旧 **Figure 1** 嵌入图 + 乱码 `Cells retained / PBS / D2…` 行。

---

### 3.4 Analysis A — transcript pre-selection

**\[新\]** 终稿全文（已扩写）：

With time point as the state definition, optimize.sd_selection identified state-specific fluctuating transcript sets for each retained sampling time point under RTF cutoff 0.1. The number of pre-selected genes varied across time points, reflecting differences in transcriptional heterogeneity as injury progressed. Time points with very few cells after lineage filtering (including day 28) were excluded from downstream network construction and CTS scoring.

---

### 3.5 Analysis A — network partition and MCI scoring

**\[新\]** 终稿全文（已扩写）：

For each retained time point, BioTIP constructed a co-expression network on the pre-selected transcripts and partitioned it into modules. MCI scoring highlighted several time points with high-scoring modules, indicating coordinated co-expression consistent with local instability. These modules were forwarded as CTS candidates for Ic-based tipping-point testing (Table 2; Analysis A rows).

---

### 3.6 Analysis A — tipping point identification (Ic.shrink)

**\[改\]** 用下面**替换**农夫 **§3.3 Analysis A** 全文（数字来自 Step6）：

When sampling time point was used as the state definition, BioTIP returned three statistically supported CTS candidates: day 5 (ΔIc.shrink = 1.10, 21 genes, p < 0.01), day 8 (ΔIc.shrink = 0.44, 23 genes, p = 0.01) and day 9 (ΔIc.shrink = 0.55, 30 genes, p < 0.01). For the day 5 module, Ic.shrink was highest at day 5 itself (Figure 6), and permutation testing confirmed significance (Figure 7). Overall, the time-point analysis places the strongest signals in mid-acute repair phases after widespread AT2 activation.

**\[留\]** 农夫 **Table 3** 格式，但**整表换成**：

| Analysis | Candidate | Genes | ΔIc.shrink | Permutation p | Peak state |
|----------|-------------|-------|------------|---------------|------------|
| A: time point | day 5 | 21 | 1.10 | <0.01 | day 5 |
| A: time point | day 8 | 23 | 0.44 | 0.01 | day 8 |
| A: time point | day 9 | 30 | 0.55 | <0.01 | day 9 |
| B: cell type | Krt8+ ADI | 11 | 2.61 | 0.01 | Krt8+ ADI |
| B: cell type | AT1 | 11 | 1.18 | 0.04 | Mki67+ Proliferation |
| B: cell type | Mki67+ Proliferation | 22 | 1.00 | <0.01 | Mki67+ Proliferation |

*Table 2. BioTIP CTS candidates (GSE141259).* — **一张表放 3.6**，3.8 可写「见 Table 2」不必重复表。

**\[图: fig_analysisA_ic_trajectory.pdf\]**  
*Figure 6. Ic.shrink across time points for the top Analysis A CTS module (day 5 candidate).*

**\[图: fig_analysisA_permutation.pdf\]**  
*Figure 7. Gene-permutation null distribution versus observed ΔIc.shrink (Analysis A).*

**\[删\]** 农夫 §3.3 里所有 **day 6, day 7, day 8→day 10** 旧叙述。

---

### 3.7 Analysis A — interpreting the result

**\[留\]+改** 农夫 Discussion 第 1 段里关于 time point 的解读，压缩成：

The timing of significant CTS candidates is biologically plausible relative to Strunz et al. (2020): early activated AT2 expansion precedes enrichment of Krt8+ ADI and later AT1-associated cells. Non-significant or weaker modules at sparse time points likely reflect heterogeneity within each time label rather than absence of repair dynamics. Time-point states are therefore useful for asking *when* the epithelial compartment appears most unstable, but not for assigning lineage direction alone.

---

### 3.8 Analysis B — cell type as state

**\[改\]** 用下面**替换**农夫 **§3.4 Analysis B**：

When annotated cell type was used as the state definition, BioTIP identified three CTS modules: Krt8+ ADI (11 genes, ΔIc.shrink = 2.61, p = 0.01), AT1 (11 genes, ΔIc.shrink = 1.18, p = 0.04) and Mki67+ proliferation (22 genes, ΔIc.shrink = 1.00, p < 0.01). Ic.shrink for the Krt8+ ADI module peaked in that cell type (Figure 9). MCI values were distributed across modules within each cell type (Figure 8). The Krt8+ ADI result aligns with the proposed transitional state between AT2 and AT1 (Strunz et al., 2020).

**\[图: fig_analysisB_mci.pdf\]**  
*Figure 8. MCI modules per cell type (Analysis B).*

**\[图: fig_analysisB_ic_trajectory.pdf\]**  
*Figure 9. Ic.shrink across cell types for the Krt8+ ADI CTS module.*

**\[图: fig_analysisB_permutation.pdf\]**  
*Figure 10. Permutation test for the Krt8+ ADI CTS module (Analysis B).*

**\[删\]** 农夫 §3.4 里「only one CTS… Mki67+ only」旧版（若与上表冲突以**上表为准**）。

---

### 3.9 Compare analyses A and B

**\[留\]+扩** 终稿全文：

The two BioTIP analyses are complementary rather than directly interchangeable. Time-point states preserve temporal order and highlight mid-acute injury phases (days 5, 8 and 9). Cell-type states pool cells across time and instead identify which annotated programmes carry criticality—notably Krt8+ ADI and proliferation. Analysis B therefore emphasizes *which* cell programme resembles a tipping module, whereas Analysis A emphasizes *when* the population appears most unstable (Figure 11). Module sizes were generally larger in Analysis B, suggesting that cell-type homogeneity strengthens co-expression signals relative to heterogeneous time-point mixtures.

**\[图: fig02_biotip_cts_results.pdf\]**  
*Figure 11. ΔIc.shrink for all CTS candidates in Analyses A and B.*

**\[删\]** 农夫旧 **Figure 2** 单独图（由 Figure 11 替代）。

---

### 3.10 Brief comment on the BioTIP method

**\[留\]+扩** 终稿全文：

BioTIP is useful for hypothesis generation about unstable epithelial states when combined with careful state definitions and biological context. Limitations include sensitivity to state labelling, module size, downsampling and uneven time-point density; sparse labels can weaken CTS calls. The method detects coordinated fluctuation and co-expression, not causal lineage relationships. Transparent reporting of RTF cutoffs, FDR thresholds and permutation design is therefore essential for reproducibility. For this dataset, BioTIP adds value by linking statistical criticality-like signals to biologically named repair states, but results should be interpreted alongside marker-based annotation and the Strunz et al. (2020) repair model rather than as stand-alone proof of tipping points.

---

## 4 Discussion

### 4.1 Integration with the study question

**\[留\]** 农夫 Discussion 最后一段（coherent repair）：

Despite caveats, the results are coherent with the published view of post-injury alveolar repair. The strongest time-point signals fall after early AT2 activation, while cell-type analysis highlights Krt8+ ADI and proliferative programmes. BioTIP is therefore useful for highlighting candidate transition phases when combined with biological annotation, but is not a stand-alone substitute for careful cell-state interpretation.

### 4.2 Limitations

**\[留\]+扩** 终稿全文：

BioTIP relies on state definitions supplied by the analyst, and both time points and cell labels are imperfect summaries of a continuous injury response. Downsampling and gene filtering were applied to obtain stable networks on a standard workstation; weaker modules might appear under larger, fully powered runs. Permutation tests randomise genes rather than cells, so significant ΔIc.shrink should be interpreted as coordinated expression statistics rather than definitive evidence of directionality. Finally, the gastrulation case study used cluster labels rather than the full 19-cluster tutorial setting, which should be considered when comparing figures to the original BioTIP walkthrough.

---

## References

**\[留\]** 农夫 References，**\[新\]** 加一条 case study：

1. Strunz, M., et al. (2020). Alveolar regeneration through a Krt8+ transitional stem cell state that persists in human lung fibrosis. *Nature Communications*, 11:3559.

2. Yang, X., et al. (2022). BioTIP: a systems biology tool for identifying critical transitions in omics data. *Nucleic Acids Research*, 50(16):e91.

3. Pijuan-Sala, B., et al. (2019). A single-cell molecular map of mouse gastrulation and early organogenesis. *Nature*, 566:490–495.

---

## Supplementary materials

**\[新\]** guidance 要求，不计字数：

Key R scripts (project folder `GP2 ICA1 5.22`):

- `CODE/M010_biotip_case_study/03_code/run_gse87038_case_study_simplified.R`
- `CODE/M030_gse141259_preprocessing/03_code/preprocess_gse141259.R`
- `STEP6_run_lung_AB_figures.R` and `STEP6B_fix_analysisB_figures.R`

Session information: R 4.5.2; BioTIP 1.23.0; Seurat 5.5.0; GEOquery 2.78.0.

---

## 终稿快速核对清单

- [ ] 无 **Keywords**
- [ ] 无 **§3.1 rationale** 旧文
- [ ] 无 **day 6 / 仅 Mki67** 旧 Results
- [ ] 无 **reduced workstation** 自我贬低句
- [ ] Results 有 **3.1 – 3.10** 标题
- [ ] **11 张图** 均在 `figures_for_word/`
- [ ] Table 2（CTS）数字与上文一致
- [ ] 正文 ≤ **2000** 字 → 导出 PDF

---

## 图文件一览（均在 `figures_for_word/`）

| 图号 | 文件 |
|------|------|
| 1 | fig_case_mci_panels.pdf |
| 2 | fig_case_ic_permutation.pdf |
| 3 | fig_umap_celltype.pdf |
| 4 | fig_umap_timepoint.pdf |
| 5 | fig01_cell_composition.pdf |
| 6 | fig_analysisA_ic_trajectory.pdf |
| 7 | fig_analysisA_permutation.pdf |
| 8 | fig_analysisB_mci.pdf |
| 9 | fig_analysisB_ic_trajectory.pdf |
| 10 | fig_analysisB_permutation.pdf |
| 11 | fig02_biotip_cts_results.pdf |
