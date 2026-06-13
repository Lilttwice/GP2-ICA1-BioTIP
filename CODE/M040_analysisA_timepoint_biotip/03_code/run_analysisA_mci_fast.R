suppressPackageStartupMessages(library(BioTIP))

output_dir <- "CODE/M040_analysisA_timepoint_biotip/04_output"
subcounts <- readRDS(file.path(output_dir, "analysisA_subcounts.rds"))
cluster <- readRDS(file.path(output_dir, "analysisA_clusters.rds"))

small_cluster <- lapply(cluster, function(x) {
  if (is.null(x) || is.null(x$membership)) {
    return(x)
  }
  x
})

membersL <- getMCI(small_cluster, subcounts, adjust.size = TRUE, fun = "BioTIP")
saveRDS(membersL, file.path(output_dir, "analysisA_membersL_fast.rds"))
message("Analysis A fast MCI step completed.")
