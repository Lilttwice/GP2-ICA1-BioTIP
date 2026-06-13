suppressPackageStartupMessages(library(BioTIP))

output_dir <- "CODE/M040_analysisA_timepoint_biotip/04_output"
subcounts <- readRDS(file.path(output_dir, "analysisA_subcounts.rds"))
igraphL <- getNetwork(subcounts, fdr = 0.1)
saveRDS(igraphL, file.path(output_dir, "analysisA_networks.rds"))
message("Analysis A network step completed.")
