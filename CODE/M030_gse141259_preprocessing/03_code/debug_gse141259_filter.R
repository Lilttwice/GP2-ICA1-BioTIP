suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
})

counts_raw <- readMM("01_data/raw/GSE141259/GSE141259_HighResolution_rawcounts.mtx.gz")
barcodes <- readLines(gzfile("01_data/raw/GSE141259/GSE141259_HighResolution_barcodes.txt.gz", "rt"))
genes <- readLines(gzfile("01_data/raw/GSE141259/GSE141259_HighResolution_genes.txt.gz", "rt"))
meta <- as.data.frame(fread("01_data/raw/GSE141259/GSE141259_HighResolution_cellinfo.csv.gz"))
rownames(meta) <- meta$cell_barcode

counts <- t(counts_raw)
rownames(counts) <- genes
colnames(counts) <- barcodes
meta <- meta[colnames(counts), , drop = FALSE]

seu <- CreateSeuratObject(counts = counts, meta.data = meta)
target_cell_types <- c("AT2", "AT2 activated", "Krt8+ ADI", "AT1", "Mki67+ Proliferation")

out <- c(
  paste("n_cells_total:", ncol(seu)),
  paste("meta_columns:", paste(colnames(seu[[]]), collapse = ", ")),
  paste("n_target_lineage:", sum(seu$cell_type %in% target_cell_types)),
  paste("summary_meta_n_genes:", paste(capture.output(summary(seu$n_genes)), collapse = " | ")),
  paste("summary_meta_n_counts:", paste(capture.output(summary(seu$n_counts)), collapse = " | ")),
  paste("summary_meta_percent_mito:", paste(capture.output(summary(seu$percent_mito)), collapse = " | ")),
  paste("summary_seurat_nFeature_RNA:", paste(capture.output(summary(seu$nFeature_RNA)), collapse = " | ")),
  paste("summary_seurat_nCount_RNA:", paste(capture.output(summary(seu$nCount_RNA)), collapse = " | ")),
  paste("pass_lineage:", sum(seu$cell_type %in% target_cell_types)),
  paste("pass_meta_filters:", sum(seu$n_genes >= 500 & seu$n_counts >= 1000 & seu$percent_mito <= 20, na.rm = TRUE)),
  paste("pass_seurat_filters:", sum(seu$nFeature_RNA >= 500 & seu$nCount_RNA >= 1000 & seu$percent_mito <= 20, na.rm = TRUE)),
  paste("pass_both_meta:", sum(seu$cell_type %in% target_cell_types & seu$n_genes >= 500 & seu$n_counts >= 1000 & seu$percent_mito <= 20, na.rm = TRUE)),
  paste("pass_both_seurat:", sum(seu$cell_type %in% target_cell_types & seu$nFeature_RNA >= 500 & seu$nCount_RNA >= 1000 & seu$percent_mito <= 20, na.rm = TRUE))
)

writeLines(out, "CODE/M030_gse141259_preprocessing/04_output/debug_filter_summary.txt")
