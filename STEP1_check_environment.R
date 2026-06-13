# GP2 ICA1 — Step 1 only: environment and data check
# Run in RStudio: Source this file, or: Rscript STEP1_check_environment.R

cat("=== GP2 ICA1 Step 1 check ===\n\n")

project_root <- getwd()
cat("Working directory:\n  ", project_root, "\n\n", sep = "")

check_path <- function(label, path) {
  ok <- file.exists(path)
  cat(sprintf("[%s] %s\n  %s\n", if (ok) "OK" else "MISSING", label, path))
  invisible(ok)
}

cat("--- Folder names (scripts expect these) ---\n")
check_path("CODE folder", "CODE")
cat("\n")

cat("--- Data files (required for figures later) ---\n")
check_path("GSE141259 raw folder", "01_data/raw/GSE141259")
check_path("GSE87038 raw folder", "01_data/raw/GSE87038")
check_path("Processed Seurat object", "01_data/processed/final_tuned_dataset.rds")
cat("\n")

cat("--- Existing analysis outputs (from earlier runs) ---\n")
check_path("Analysis A CTS summary", "CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_simplified_cts_summary.tsv")
check_path("Analysis B CTS summary", "CODE/M060_analysisB_cellstate_biotip/04_output/analysisB_simplified_cts_summary.tsv")
cat("\n")

cat("--- R packages ---\n")
pkgs <- c("BioTIP", "Seurat", "GEOquery", "data.table", "ggplot2", "Matrix", "igraph")
for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    ver <- as.character(packageVersion(p))
    cat(sprintf("[OK] %s %s\n", p, ver))
  } else {
    cat(sprintf("[MISSING] %s — install before Step 2\n", p))
  }
}

cat("\nR version: ", R.version.string, "\n", sep = "")
cat("\n=== End Step 1 check ===\n")
cat("Copy ALL text above and send to Cursor for verification.\n")
