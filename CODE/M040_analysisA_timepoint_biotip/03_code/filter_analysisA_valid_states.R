suppressPackageStartupMessages({
  library(Seurat)
  library(BioTIP)
})

output_dir <- "CODE/M040_analysisA_timepoint_biotip/04_output"
run_dir <- "CODE/M040_analysisA_timepoint_biotip/05_logs"

seu <- readRDS("01_data/processed/final_tuned_dataset.rds")
counts <- GetAssayData(seu, assay = "RNA", layer = "data")
counts <- as.matrix(counts)

all_states <- as.character(seu$time_point)
samplesL <- split(colnames(seu), all_states)
samplesL <- samplesL[order(names(samplesL))]

valid_states <- names(samplesL)[lengths(samplesL) >= 20]
samplesL_valid <- samplesL[valid_states]
counts_valid <- counts[, unlist(samplesL_valid), drop = FALSE]

writeLines(
  c(
    paste("valid_states:", paste(valid_states, collapse = ", ")),
    paste("excluded_states:", paste(setdiff(names(samplesL), valid_states), collapse = ", ")),
    paste("valid_cells:", ncol(counts_valid))
  ),
  file.path(output_dir, "analysisA_valid_states.txt")
)

subcounts <- optimize.sd_selection(
  counts_valid,
  samplesL_valid,
  cutoff = 0.2,
  percent = 0.8,
  B = 10
)
saveRDS(subcounts, file.path(output_dir, "analysisA_subcounts_validstates.rds"))

igraphL <- getNetwork(subcounts, fdr = 0.1)
saveRDS(igraphL, file.path(output_dir, "analysisA_networks_validstates.rds"))

cluster <- getCluster_methods(igraphL)
saveRDS(cluster, file.path(output_dir, "analysisA_clusters_validstates.rds"))

writeLines(
  c(
    paste("run_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "status: valid-state objects created"
  ),
  file.path(run_dir, "analysisA_validstates.log")
)

message("Analysis A valid-state objects completed.")
