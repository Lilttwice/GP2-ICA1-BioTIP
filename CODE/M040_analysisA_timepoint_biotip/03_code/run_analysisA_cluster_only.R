suppressPackageStartupMessages(library(BioTIP))

output_dir <- "CODE/M040_analysisA_timepoint_biotip/04_output"
igraphL <- readRDS(file.path(output_dir, "analysisA_networks.rds"))
cluster <- getCluster_methods(igraphL)
saveRDS(cluster, file.path(output_dir, "analysisA_clusters.rds"))
message("Analysis A clustering step completed.")
