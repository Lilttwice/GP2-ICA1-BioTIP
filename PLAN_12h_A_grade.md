# 12 小时冲 A 时间表（图 + Word 一起）

## Mac 不是不行

- 32 GB 内存对 scRNA + BioTIP **够用**；你之前爆的是 **23 个 cluster × 大矩阵 × bootstrap 同时占满**，不是 Mac 输给 Windows。
- R 有 **32Gb 向量上限**；Windows 上 R 同样会撞，除非机器 RAM 更大或参数更小。
- **策略**：肺数据 A、B **分开跑**；case study 保持已成功的 10-cluster 版即可。

---

| 时段 | 任务 | 你要做 |
|------|------|--------|
| **现在 ~+1.5h** | **Step 6** 肺数据 A/B 重跑 + 新图 | `source("STEP6_run_lung_AB_figures.R")` |
| **+1.5h ~+4h** | **Step 5** Word 终稿 | `STEP5_word_assembly.md` + 插入 `figures_for_word/` 全部 PDF |
| **+4h ~+6h** | **写作冲 A** | 3.1–3.10 标题齐；Discussion 深化（见下） |
| **+6h ~+7h** | 字数 ≤2000 + Supplementary 代码路径 | 导出 `GP2_ICA1_提交.pdf` |
| **+7h ~+12h** | 缓冲 / 可选润色 | 通读、改 typo、睡一会再交 |

---

## Step 6 成功后会多这些图（`figures_for_word/`）

- `fig_analysisA_ic_trajectory.pdf`
- `fig_analysisA_permutation.pdf`
- `fig_analysisB_ic_trajectory.pdf`
- `fig_analysisB_permutation.pdf`
- `fig_analysisB_mci.pdf`
- 更新版 `fig02_biotip_cts_results.pdf`

加上已有的 case study + UMAP + composition ≈ **10 张图**（A 档厚度）。

---

## Methods 写法（专业、不「自曝偷懒」）

**Case study (3.1):** ten largest mesodermal clusters; 1,200 genes; ≤35 cells/cluster; B = 2; 50 permutations.

**Lung (2.4–2.5):** time_point / cell_type as states; ≤80 cells/state; 1,500 genes; B = 3; 100 permutations; FDR 0.2; min module 20.

---

## Discussion 冲 A 必写 3 点（每点 3–4 句）

1. **Strunz 轨迹**：day 6/9 与 AT2 激活→ADI 时间线是否一致。  
2. **A vs B**：时间点抓「何时」，细胞类型抓「Mki67+ 程序」——互补不可直接比 p 值。  
3. **局限**：state 人为、稀疏时间点、模块基因数——BioTIP 是统计信号非因果。

---

## 若 Step 6 再内存报错

```r
gc()
source("STEP6_run_lung_AB_figures.R")
```

或把脚本里 `max_cells_per_state` 改成 `60`，`max_genes` 改成 `1200` 后重跑。
