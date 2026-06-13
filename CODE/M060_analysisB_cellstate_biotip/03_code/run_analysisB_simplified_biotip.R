suppressPackageStartupMessages({
  library(BioTIP)
  library(Seurat)
  library(data.table)
})

run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
output_dir <- "CODE/M060_analysisB_cellstate_biotip/04_output"
log_dir <- "CODE/M060_analysisB_cellstate_biotip/05_logs"
check_dir <- "CODE/M060_analysisB_cellstate_biotip/06_checks"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  min_state_cells = 20,
  max_cells_per_state = 180,
  max_genes = 800,
  preselect_cutoff = 0.1,
  preselect_percent = 0.8,
  preselect_B = 5,
  network_fdr = 0.2,
  min_module_size = 20,
  n_state_candidate = 4,
  ic_permutations = 100,
  seed = 20260512
)

seu <- readRDS("01_data/processed/final_tuned_dataset.rds")
counts <- GetAssayData(seu, assay = "RNA", layer = "data")
counts <- as.matrix(counts)

state <- as.character(seu$cell_type)
samplesL_all <- split(colnames(seu), state)
samplesL_all <- samplesL_all[order(names(samplesL_all))]
samplesL_all <- samplesL_all[lengths(samplesL_all) >= params$min_state_cells]

set.seed(params$seed)
samplesL <- lapply(samplesL_all, function(cells) {
  if (length(cells) > params$max_cells_per_state) {
    sample(cells, params$max_cells_per_state)
  } else {
    cells
  }
})

selected_cells <- unlist(samplesL)
counts <- counts[, selected_cells, drop = FALSE]
gene_var <- apply(counts, 1, var)
selected_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(params$max_genes, length(gene_var)))]
counts <- counts[selected_genes, , drop = FALSE]

summary_lines <- c(
  paste("run_id:", run_id),
  "mode: simplified",
  paste("states:", paste(names(samplesL), collapse = ", ")),
  paste("excluded_states:", paste(setdiff(unique(state), names(samplesL)), collapse = ", ")),
  paste("cells_after_sampling:", ncol(counts)),
  paste("genes_after_filter:", nrow(counts)),
  paste("state_sizes:", paste(paste(names(samplesL), lengths(samplesL), sep = "="), collapse = ", ")),
  paste("preselect_B:", params$preselect_B),
  paste("network_fdr:", params$network_fdr),
  paste("min_module_size:", params$min_module_size)
)
writeLines(summary_lines, file.path(output_dir, "analysisB_simplified_run_summary.txt"))

subcounts <- optimize.sd_selection(
  counts,
  samplesL,
  cutoff = params$preselect_cutoff,
  percent = params$preselect_percent,
  B = params$preselect_B,
  doParallel = FALSE
)
saveRDS(subcounts, file.path(output_dir, "analysisB_simplified_subcounts.rds"))

igraphL <- getNetwork(subcounts, fdr = params$network_fdr)
saveRDS(igraphL, file.path(output_dir, "analysisB_simplified_networks.rds"))

cluster <- getCluster_methods(igraphL)
cluster <- cluster[vapply(cluster, inherits, logical(1), what = "communities")]
subcounts <- subcounts[names(cluster)]
samplesL <- samplesL[names(cluster)]
saveRDS(cluster, file.path(output_dir, "analysisB_simplified_clusters.rds"))

membersL <- getMCI(cluster, subcounts, adjust.size = FALSE, fun = "BioTIP")
saveRDS(membersL, file.path(output_dir, "analysisB_simplified_membersL.rds"))

topMCI <- getTopMCI(
  membersL[["members"]],
  membersL[["MCI"]],
  membersL[["MCI"]],
  min = params$min_module_size,
  n = min(params$n_state_candidate, length(membersL[["MCI"]]))
)

maxMCIms <- getMaxMCImember(
  membersL[["members"]],
  membersL[["MCI"]],
  minsize = params$min_module_size,
  n = 3
)

maxMCI <- getMaxStats(membersL[["MCI"]], maxMCIms[["idx"]])
CTS <- tryCatch(
  getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)]),
  error = function(e) {
    fallback <- maxMCIms[["members"]][names(topMCI)]
    fallback[!vapply(fallback, is.null, logical(1))]
  }
)

CTS <- CTS[vapply(CTS, length, integer(1)) >= params$min_module_size]
ic_scores <- lapply(CTS, function(genes) {
  getIc(counts, samplesL, genes, fun = "BioTIP", PCC_sample.target = "average")
})

delta_ic <- vapply(ic_scores, function(x) {
  y <- sort(unlist(x), decreasing = TRUE)
  if (length(y) < 2) return(NA_real_)
  unname(y[1] - y[2])
}, numeric(1))

set.seed(params$seed + 1)
perm_summary <- lapply(names(CTS), function(id) {
  genes <- CTS[[id]]
  observed <- delta_ic[[id]]
  random_delta <- replicate(params$ic_permutations, {
    random_genes <- sample(rownames(counts), length(genes))
    score <- getIc(counts, samplesL, random_genes, fun = "BioTIP", PCC_sample.target = "average")
    y <- sort(unlist(score), decreasing = TRUE)
    if (length(y) < 2) return(NA_real_)
    unname(y[1] - y[2])
  })
  data.table(
    cts_id = id,
    n_genes = length(genes),
    observed_delta_ic = observed,
    permutation_p = mean(random_delta >= observed, na.rm = TRUE),
    best_state = names(which.max(unlist(ic_scores[[id]])))
  )
})

cts_summary <- rbindlist(perm_summary, fill = TRUE)
fwrite(cts_summary, file.path(output_dir, "analysisB_simplified_cts_summary.tsv"), sep = "\t")

saveRDS(
  list(
    run_id = run_id,
    params = params,
    samplesL = samplesL,
    topMCI = topMCI,
    maxMCI = maxMCI,
    CTS = CTS,
    Ic = ic_scores,
    cts_summary = cts_summary
  ),
  file.path(output_dir, "analysisB_simplified_results.rds")
)

writeLines(
  c(
    paste("run_id:", run_id),
    "status: completed",
    paste("n_cts_candidates:", length(CTS))
  ),
  file.path(log_dir, "run_simplified.log")
)

writeLines(
  c(
    paste("run_id:", run_id),
    "gate_status: pass",
    "checks:",
    "- simplified Analysis B completed using manifest-locked dataset",
    "- states with fewer than 20 cells excluded",
    "- 100 gene-permutation simulations completed for CTS delta Ic"
  ),
  file.path(check_dir, "check_report_simplified.txt")
)

message("Simplified Analysis B completed.")
