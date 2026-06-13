suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

fig_dir <- "03_reports/manuscript/figures"
tab_dir <- "03_reports/manuscript/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

composition <- fread("CODE/M030_gse141259_preprocessing/04_output/timepoint_celltype_counts_final.tsv")
composition[, time_clean := gsub("day ", "D", time_point, fixed = TRUE)]
composition[time_point == "d14_PBS", time_clean := "PBS"]
composition[, day_num := fifelse(time_clean == "PBS", -1L, as.integer(gsub("D", "", time_clean)))]
composition <- composition[order(day_num)]
composition[, time_clean := factor(time_clean, levels = unique(time_clean))]

palette_cells <- c(
  "AT2" = "#4C78A8",
  "AT2 activated" = "#F58518",
  "Krt8+ ADI" = "#54A24B",
  "AT1" = "#E45756",
  "Mki67+ Proliferation" = "#72B7B2"
)

p1 <- ggplot(composition, aes(x = time_clean, y = n_cells, fill = cell_type)) +
  geom_col(width = 0.74, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = palette_cells) +
  labs(x = "Sampling state", y = "Cells retained", fill = "Cell state") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "#1F2933"),
    axis.text.y = element_text(color = "#1F2933"),
    axis.title = element_text(color = "#1F2933"),
    legend.position = "bottom",
    legend.title = element_text(color = "#1F2933"),
    plot.margin = margin(7, 8, 7, 8)
  )

grDevices::cairo_pdf(file.path(fig_dir, "fig01_cell_composition.pdf"), width = 7.1, height = 4.2)
print(p1)
dev.off()

cts <- fread("CODE/M070_compare_and_comment/04_output/simplified_AB_cts_summary.tsv")
cts[, label := fifelse(
  analysis == "A_time_point",
  paste0("A: ", cts_id),
  paste0("B: ", cts_id)
)]
cts[, significant := permutation_p <= 0.05]
cts[, label := factor(label, levels = label[order(analysis, -observed_delta_ic)])]

p2 <- ggplot(cts, aes(x = label, y = observed_delta_ic, fill = significant)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.15) +
  geom_text(aes(label = paste0("p=", format(permutation_p, digits = 2))), hjust = -0.08, size = 3.1, color = "#1F2933") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("TRUE" = "#2F6F73", "FALSE" = "#B9C0C7"), labels = c("FALSE" = "not significant", "TRUE" = "significant")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = expression(Delta~Ic.shrink), fill = "Permutation result") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "#1F2933"),
    axis.title = element_text(color = "#1F2933"),
    legend.position = "bottom",
    plot.margin = margin(7, 35, 7, 8)
  )

grDevices::cairo_pdf(file.path(fig_dir, "fig02_biotip_cts_results.pdf"), width = 7.1, height = 3.8)
print(p2)
dev.off()

dataset_summary <- data.table(
  Item = c(
    "Primary dataset",
    "Input matrix",
    "Cells before lineage filtering",
    "Cells retained for analysis",
    "Retained cell states",
    "Analysis A states",
    "Analysis B states"
  ),
  Value = c(
    "GSE141259 lung injury regeneration scRNA-seq",
    "HighResolution raw-count matrix with matched cell metadata",
    "32,559",
    "10,381",
    "AT2, activated AT2, Krt8+ ADI, AT1, Mki67+ proliferating cells",
    "Time points after injury; day 28 excluded from BioTIP scoring because only two cells remained",
    "Annotated alveolar cell states pooled across time"
  )
)

write_booktabs <- function(dt, path, caption, label, align) {
  lines <- c(
    paste0("\\begin{table}[H]"),
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(names(dt), collapse = " & "),
    "\\\\",
    "\\midrule"
  )
  body <- apply(dt, 1, function(row) paste(row, collapse = " & "))
  lines <- c(lines, paste0(body, "\\\\"), "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
}

write_booktabs(
  dataset_summary,
  file.path(tab_dir, "tab01_dataset_summary.tex"),
  "Dataset and analysis scope.",
  "tab:dataset-summary",
  "p{0.31\\linewidth}p{0.60\\linewidth}"
)

cts_table <- copy(cts)
cts_table[, Analysis := fifelse(analysis == "A_time_point", "A: time point", "B: cell type")]
cts_table[, Candidate := cts_id]
cts_table[, `Genes` := n_genes]
cts_table[, `Delta Ic.shrink` := sprintf("%.3f", observed_delta_ic)]
cts_table[, `Permutation p` := fifelse(permutation_p == 0, "<0.01", sprintf("%.2f", permutation_p))]
cts_table[, `Peak state` := best_state]
cts_table <- cts_table[, .(Analysis, Candidate, Genes, `Delta Ic.shrink`, `Permutation p`, `Peak state`)]

write_booktabs(
  cts_table,
  file.path(tab_dir, "tab02_biotip_results.tex"),
  "BioTIP critical-transition signal candidates from the simplified analysis.",
  "tab:biotip-results",
  "p{0.18\\linewidth}p{0.21\\linewidth}r r r p{0.16\\linewidth}"
)

software <- data.table(
  Software = c("R", "Seurat", "BioTIP", "GEOquery"),
  Version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("Seurat")),
    as.character(packageVersion("BioTIP")),
    as.character(packageVersion("GEOquery"))
  ),
  Use = c(
    "Analysis environment",
    "Single-cell object handling, normalization, PCA and UMAP",
    "RTF pre-selection, network partition, MCI and Ic.shrink scoring",
    "GEO metadata and supplementary-file retrieval"
  )
)

write_booktabs(
  software,
  file.path(tab_dir, "tab03_software.tex"),
  "Software used for the analysis.",
  "tab:software",
  "p{0.18\\linewidth}p{0.16\\linewidth}p{0.56\\linewidth}"
)

message("Report assets generated.")
