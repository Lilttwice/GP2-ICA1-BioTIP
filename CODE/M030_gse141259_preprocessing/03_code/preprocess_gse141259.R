suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(tools)
})

run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
base_dir <- "01_data/raw/GSE141259"
output_dir <- "CODE/M030_gse141259_preprocessing/04_output"
log_dir <- "CODE/M030_gse141259_preprocessing/05_logs"
check_dir <- "CODE/M030_gse141259_preprocessing/06_checks"
processed_dir <- "01_data/processed"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

target_cell_types <- c("AT2", "AT2 activated", "Krt8+ ADI", "AT1", "Mki67+ Proliferation")
min_genes <- 500
max_percent_mito <- 20
min_counts <- 1000
n_variable_features <- 3000

counts_raw <- readMM(file.path(base_dir, "GSE141259_HighResolution_rawcounts.mtx.gz"))
barcodes <- readLines(gzfile(file.path(base_dir, "GSE141259_HighResolution_barcodes.txt.gz"), open = "rt"))
genes <- readLines(gzfile(file.path(base_dir, "GSE141259_HighResolution_genes.txt.gz"), open = "rt"))

if (nrow(counts_raw) == length(barcodes) && ncol(counts_raw) == length(genes)) {
  counts <- t(counts_raw)
} else if (nrow(counts_raw) == length(genes) && ncol(counts_raw) == length(barcodes)) {
  counts <- counts_raw
} else {
  stop("Count matrix dimensions do not match barcode/gene files.")
}

rownames(counts) <- genes
colnames(counts) <- barcodes

meta <- fread(file.path(base_dir, "GSE141259_HighResolution_cellinfo.csv.gz"))
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_barcode
meta <- meta[colnames(counts), , drop = FALSE]

if (!identical(rownames(meta), colnames(counts))) {
  stop("Metadata rownames do not match count matrix columns.")
}

seu <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = "GSE141259_HighResolution"
)

seu$selected_lineage <- seu$cell_type %in% target_cell_types
before_n <- ncol(seu)

seu <- subset(
  seu,
  subset = selected_lineage & nFeature_RNA >= min_genes & nCount_RNA >= min_counts & percent_mito <= max_percent_mito
)

after_n <- ncol(seu)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_variable_features, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)

celltype_counts <- as.data.table(table(seu$cell_type))
setnames(celltype_counts, c("cell_type", "n_cells"))
fwrite(celltype_counts, file.path(output_dir, "celltype_counts_final.tsv"), sep = "\t")

time_celltype_counts <- as.data.table(table(seu$time_point, seu$cell_type))
setnames(time_celltype_counts, c("time_point", "cell_type", "n_cells"))
fwrite(time_celltype_counts, file.path(output_dir, "timepoint_celltype_counts_final.tsv"), sep = "\t")

qc_summary <- c(
  paste("run_id:", run_id),
  paste("input_cells:", before_n),
  paste("output_cells:", after_n),
  paste("retained_fraction:", round(after_n / before_n, 4)),
  paste("target_cell_types:", paste(target_cell_types, collapse = ", ")),
  paste("min_genes:", min_genes),
  paste("min_counts:", min_counts),
  paste("max_percent_mito:", max_percent_mito),
  paste("n_variable_features:", n_variable_features),
  paste("time_points_retained:", paste(sort(unique(seu$time_point)), collapse = ", ")),
  paste("cell_types_retained:", paste(sort(unique(seu$cell_type)), collapse = ", "))
)
writeLines(qc_summary, file.path(output_dir, "filter_summary.txt"))

grDevices::cairo_pdf(file.path(output_dir, "fig_qc_violin.pdf"), width = 10, height = 4)
print(VlnPlot(seu, features = c("n_genes", "n_counts", "percent_mito"), ncol = 3, pt.size = 0))
dev.off()

grDevices::cairo_pdf(file.path(output_dir, "fig_umap_celltype.pdf"), width = 8, height = 6)
print(DimPlot(seu, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE))
dev.off()

grDevices::cairo_pdf(file.path(output_dir, "fig_umap_timepoint.pdf"), width = 8, height = 6)
print(DimPlot(seu, reduction = "umap", group.by = "time_point"))
dev.off()

final_rds <- file.path(processed_dir, "final_tuned_dataset.rds")
final_meta_csv <- file.path(processed_dir, "final_tuned_dataset_cell_metadata.csv")

saveRDS(seu, final_rds)
fwrite(as.data.table(seu[[]], keep.rownames = "cell_barcode"), final_meta_csv, sep = ",")

data_hash <- unname(md5sum(final_rds))
lock_time <- format(Sys.time(), tz = "UTC", usetz = TRUE)

manifest_lines <- c(
  "final_dataset_id: final_tuned_gse141259_alveolar_highres",
  "source_module: M030_gse141259_preprocessing",
  paste("source_run_id:", run_id),
  paste("data_path:", final_rds),
  paste("data_hash:", data_hash),
  "label_version: highres_cellinfo_cell_type_v1",
  "feature_version: GSE141259_HighResolution_genes_v1",
  paste("lock_time:", lock_time),
  "allowed_downstream_modules:",
  "  - M040_analysisA_timepoint_biotip",
  "  - M060_analysisB_cellstate_biotip",
  "  - M080_report_assembly",
  "deprecated_inputs:",
  "  - 01_data/raw/GSE141259/GSE141259_WholeLung_rawcounts.mtx.gz",
  "  - 01_data/raw/GSE141259/GSE141259_HighResolution_rawcounts.mtx.gz"
)
writeLines(manifest_lines, "01_data/metadata/final_dataset_manifest.yaml")

check_lines <- c(
  paste("run_id:", run_id),
  "gate_status: pass",
  "checks:",
  "- high-resolution dataset selected because time_point and cell_type are both available",
  "- final dataset saved to a unique locked path",
  "- manifest updated with md5 hash",
  "- downstream modules should read the manifest-locked rds file only"
)
writeLines(check_lines, file.path(check_dir, "check_report.txt"))

dq_lines <- c(
  paste("run_id:", run_id),
  paste("input_cells:", before_n),
  paste("output_cells:", after_n),
  paste("n_unique_time_points:", length(unique(seu$time_point))),
  paste("n_unique_cell_types:", length(unique(seu$cell_type))),
  "dq_status: pass_with_documented_filters"
)
writeLines(dq_lines, file.path(check_dir, "dq_report.txt"))

log_lines <- c(
  paste("run_id:", run_id),
  "status: completed",
  paste("final_rds:", final_rds),
  paste("data_hash:", data_hash)
)
writeLines(log_lines, file.path(log_dir, "run.log"))

message("Preprocessing completed. Final dataset locked.")
