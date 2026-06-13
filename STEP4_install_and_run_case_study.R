# Step 4 — install dependencies, then run GSE87038 BioTIP case study
setwd("/Users/liutongtong/Desktop/GP2 ICA1 5.22")

# Avoid interactive "create directory? (yes/no):" from ExperimentHub / MouseGastrulationData
eh_cache <- path.expand("~/Library/Caches/org.R-project.R/R/ExperimentHub")
if (!dir.exists(eh_cache)) {
  dir.create(eh_cache, recursive = TRUE, showWarnings = FALSE)
  message("Created ExperimentHub cache: ", eh_cache)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

need_bioc <- c(
  "MouseGastrulationData", "scran", "scater", "SingleCellExperiment",
  "BiocSingular", "BioTIP", "GEOquery"
)
for (p in need_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " ...")
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

source("02_modules/M010_biotip_case_study/03_code/run_gse87038_case_study_simplified.R")

cat("\n=== Step 4 done if you see fig_case_*.pdf below ===\n")
list.files(
  "02_modules/M010_biotip_case_study/04_output",
  pattern = "fig_case|case_study",
  full.names = TRUE
)
