suppressPackageStartupMessages(library(BioTIP))

output_dir <- "CODE/M040_analysisA_timepoint_biotip/04_output"
subcounts <- readRDS(file.path(output_dir, "analysisA_subcounts.rds"))
cluster <- readRDS(file.path(output_dir, "analysisA_clusters.rds"))
membersL <- getMCI(cluster, subcounts, adjust.size = FALSE, fun = "BioTIP")
saveRDS(membersL, file.path(output_dir, "analysisA_membersL.rds"))
message("Analysis A MCI step completed.")
