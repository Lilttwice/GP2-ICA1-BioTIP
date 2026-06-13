suppressPackageStartupMessages({
  library(BioTIP)
  library(Seurat)
})

run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
output_dir <- "CODE/M040_analysisA_timepoint_biotip/04_output"
log_dir <- "CODE/M040_analysisA_timepoint_biotip/05_logs"
check_dir <- "CODE/M040_analysisA_timepoint_biotip/06_checks"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS("01_data/processed/final_tuned_dataset.rds")
counts <- GetAssayData(seu, assay = "RNA", layer = "data")
counts <- as.matrix(counts)

states <- as.character(seu$time_point)
samplesL <- split(colnames(seu), states)
samplesL <- samplesL[order(names(samplesL))]

params <- list(
  cutoff = 0.2,
  percent = 0.8,
  B = 10,
  network_fdr = 0.1,
  min_size = 50,
  n_state_candidate = 4,
  fun = "BioTIP",
  ic_simulations = 100
)

writeLines(
  c(
    paste("run_id:", run_id),
    paste("n_cells:", ncol(counts)),
    paste("n_genes:", nrow(counts)),
    paste("states:", paste(names(samplesL), collapse = ", ")),
    paste("state_sizes:", paste(paste(names(samplesL), lengths(samplesL), sep = "="), collapse = ", ")),
    paste("cutoff:", params$cutoff),
    paste("percent:", params$percent),
    paste("B:", params$B),
    paste("network_fdr:", params$network_fdr),
    paste("min_size:", params$min_size),
    paste("n_state_candidate:", params$n_state_candidate),
    paste("ic_simulations:", params$ic_simulations)
  ),
  file.path(output_dir, "analysisA_run_summary.txt")
)

subcounts_path <- file.path(output_dir, "analysisA_subcounts.rds")
networks_path <- file.path(output_dir, "analysisA_networks.rds")
clusters_path <- file.path(output_dir, "analysisA_clusters.rds")
members_path <- file.path(output_dir, "analysisA_membersL.rds")
cts_path <- file.path(output_dir, "analysisA_cts_and_scores.rds")

set.seed(20260511)
if (file.exists(subcounts_path)) {
  subcounts <- readRDS(subcounts_path)
} else {
  subcounts <- optimize.sd_selection(
    counts,
    samplesL,
    cutoff = params$cutoff,
    percent = params$percent,
    B = params$B
  )
  saveRDS(subcounts, subcounts_path)
}

if (file.exists(networks_path)) {
  igraphL <- readRDS(networks_path)
} else {
  igraphL <- getNetwork(subcounts, fdr = params$network_fdr)
  saveRDS(igraphL, networks_path)
}

if (file.exists(clusters_path)) {
  cluster <- readRDS(clusters_path)
} else {
  cluster <- getCluster_methods(igraphL)
  saveRDS(cluster, clusters_path)
}

if (file.exists(members_path)) {
  membersL <- readRDS(members_path)
} else {
  membersL <- getMCI(cluster, subcounts, adjust.size = FALSE, fun = params$fun)
  saveRDS(membersL, members_path)
}

topMCI <- getTopMCI(
  membersL[["members"]],
  membersL[["MCI"]],
  membersL[["MCI"]],
  min = params$min_size,
  n = params$n_state_candidate
)

maxMCIms <- getMaxMCImember(
  membersL[["members"]],
  membersL[["MCI"]],
  min = params$min_size,
  n = 3
)

maxMCI <- getMaxStats(membersL[["MCI"]], maxMCIms[["idx"]])
CTS <- getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])

tmp <- unlist(lapply(maxMCIms[["idx"]][names(topMCI)], length))
whoistop2nd <- names(tmp[tmp == 2])
whoistop3rd <- names(tmp[tmp == 3])

if (length(whoistop2nd) > 0) {
  CTS <- append(CTS, maxMCIms[["2topest.members"]][whoistop2nd])
}
if (length(whoistop3rd) > 0) {
  CTS <- append(CTS, maxMCIms[["2topest.members"]][whoistop3rd])
  CTS <- append(CTS, maxMCIms[["3topest.members"]][whoistop3rd])
}

maxMCI_used <- maxMCI[names(topMCI)]
if (length(whoistop2nd) > 0) {
  maxMCI_used <- c(maxMCI_used, getNextMaxStats(membersL[["MCI"]], idL = maxMCIms[["idx"]], whoistop2nd))
}
if (length(whoistop3rd) > 0) {
  maxMCI_used <- c(maxMCI_used, getNextMaxStats(membersL[["MCI"]], idL = maxMCIms[["idx"]], whoistop3rd))
  maxMCI_used <- c(maxMCI_used, getNextMaxStats(membersL[["MCI"]], idL = maxMCIms[["idx"]], whoistop3rd, which.next = 3))
}

if (file.exists(cts_path)) {
  cts_result <- readRDS(cts_path)
  CTS <- cts_result$CTS
  cts_scores <- cts_result$Ic
} else {
  cts_scores <- lapply(CTS, function(x) getIc(counts, samplesL, x, fun = params$fun, PCC_sample.target = "average"))
  saveRDS(
    list(
      topMCI = topMCI,
      maxMCIms = maxMCIms,
      maxMCI = maxMCI,
      maxMCI_used = maxMCI_used,
      CTS = CTS,
      Ic = cts_scores
    ),
    cts_path
  )
}

cts_summary <- data.frame(
  cts_id = names(CTS),
  n_genes = vapply(CTS, length, integer(1)),
  stringsAsFactors = FALSE
)
write.table(
  cts_summary,
  file = file.path(output_dir, "analysisA_cts_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  c(
    paste("run_id:", run_id),
    "status: completed",
    paste("n_cts_candidates:", length(CTS)),
    paste("top_states:", paste(names(topMCI), collapse = ", "))
  ),
  file.path(log_dir, "run.log")
)

writeLines(
  c(
    paste("run_id:", run_id),
    "gate_status: provisional_pass",
    "checks:",
    "- Analysis A completed through CTS and Ic score calculation",
    "- 100-gene permutation significance is pending implementation in a separate verification step",
    "- downstream interpretation should use only finalized outputs from this run"
  ),
  file.path(check_dir, "check_report.txt")
)

message("Analysis A core pipeline completed.")
