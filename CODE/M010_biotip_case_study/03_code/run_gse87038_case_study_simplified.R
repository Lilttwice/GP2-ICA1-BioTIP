# M010 — BioTIP successful case study (GSE87038 E8.25, memory-light)
# Produces fig_case_mci_panels.pdf and fig_case_ic_permutation.pdf

suppressPackageStartupMessages({
  library(BioTIP)
  library(MouseGastrulationData)
  library(scran)
  library(scater)
  library(SingleCellExperiment)
  library(igraph)
  library(data.table)
})

try(mem.maxVSize(48 * 1024^3), silent = TRUE)
gc()

output_dir <- "CODE/M010_biotip_case_study/04_output"
log_dir <- "CODE/M010_biotip_case_study/05_logs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  mesoderm_types = c(
    "Somitic mesoderm", "Intermediate mesoderm", "ExE mesoderm",
    "Paraxial mesoderm", "Allantois", "Pharyngeal mesoderm", "Cardiomyocytes",
    "Mesenchyme", "Haematoendothelial progenitors", "Endothelium",
    "Blood progenitors 1", "Blood progenitors 2"
  ),
  snn_k = 10,
  max_states = 10,
  min_cells_per_state = 40,
  preselect_cutoff = 0.1,
  preselect_percent = 0.8,
  preselect_B = 2,
  network_fdr = 0.2,
  min_module_size = 40,
  n_cts_candidates = 4,
  ic_permutations = 50,
  max_cells_per_state = 35,
  max_genes = 1200,
  mci_threshold = 2,
  seed = 20260525
)

uniquifyFeatureNames <- function(ensembl, symbol) {
  symbol <- make.unique(as.character(symbol))
  names(symbol) <- ensembl
  symbol
}

message("Loading E8.25 samples 24, 25, 28 (cached if already downloaded)...")
sce.24 <- EmbryoAtlasData(samples = 24)
rownames(sce.24) <- uniquifyFeatureNames(rowData(sce.24)$ENSEMBL, rowData(sce.24)$SYMBOL)
sce.25 <- EmbryoAtlasData(samples = 25)
rownames(sce.25) <- uniquifyFeatureNames(rowData(sce.25)$ENSEMBL, rowData(sce.25)$SYMBOL)
sce.28 <- EmbryoAtlasData(samples = 28)
rownames(sce.28) <- uniquifyFeatureNames(rowData(sce.28)$ENSEMBL, rowData(sce.28)$SYMBOL)

drop.24 <- sce.24$doublet | sce.24$stripped
drop.25 <- sce.25$doublet | sce.25$stripped
drop.28 <- sce.28$doublet | sce.28$stripped
sce.24 <- sce.24[, !drop.24]
sce.25 <- sce.25[, !drop.25]
sce.28 <- sce.28[, !drop.28]

sce.24 <- logNormCounts(sce.24)
sce.25 <- logNormCounts(sce.25)
sce.28 <- logNormCounts(sce.28)

corrected <- cbind(sce.24, sce.25, sce.28)
rm(sce.24, sce.25, sce.28)
gc()

corrected <- corrected[, corrected$celltype %in% params$mesoderm_types]
dec <- modelGeneVar(corrected)
hvgs <- rownames(dec)[dec$bio > 0]
corrected <- corrected[hvgs, ]
rm(dec)
gc()

corrected <- runPCA(corrected, ncomponents = 20)
g <- buildSNNGraph(corrected, k = params$snn_k, use.dimred = "PCA")
clust <- cluster_walktrap(g)$membership
colLabels(corrected) <- factor(clust)
rm(g)
gc()

cell_labels <- colLabels(corrected)
logmat <- as.matrix(assay(corrected, "logcounts"))
rm(corrected)
gc()

samplesL_all <- split(colnames(logmat), cell_labels)
samplesL_all <- samplesL_all[lengths(samplesL_all) >= params$min_cells_per_state]
ord <- order(lengths(samplesL_all), decreasing = TRUE)
samplesL_all <- samplesL_all[ord[seq_len(min(params$max_states, length(ord)))]]

set.seed(params$seed)
samplesL <- lapply(samplesL_all, function(cells) {
  if (length(cells) > params$max_cells_per_state) {
    sample(cells, params$max_cells_per_state)
  } else {
    cells
  }
})

selected_cells <- unlist(samplesL, use.names = FALSE)
logmat <- logmat[, selected_cells, drop = FALSE]
gene_var <- apply(logmat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(params$max_genes, length(gene_var)))]
logmat <- logmat[top_genes, , drop = FALSE]
rm(gene_var)
gc()

message(
  "Running BioTIP on ", length(samplesL), " clusters (",
  nrow(logmat), " genes x ", ncol(logmat), " cells)..."
)

subcounts <- optimize.sd_selection(
  logmat,
  samplesL,
  cutoff = params$preselect_cutoff,
  percent = params$preselect_percent,
  B = params$preselect_B,
  doParallel = FALSE,
  n_cores = 1
)
gc()

igraphL <- getNetwork(subcounts, fdr = params$network_fdr)
cluster <- getCluster_methods(igraphL)
cluster <- cluster[vapply(cluster, inherits, logical(1), what = "communities")]
subcounts <- subcounts[names(cluster)]
samplesL <- samplesL[names(cluster)]
rm(igraphL)
gc()

membersL <- getMCI(cluster, subcounts, adjust.size = FALSE, fun = "BioTIP")
topMCI <- getTopMCI(
  membersL[["members"]],
  membersL[["MCI"]],
  membersL[["MCI"]],
  min = params$min_module_size,
  n = min(params$n_cts_candidates, length(membersL[["MCI"]]))
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
    message("getCTS fallback: ", conditionMessage(e))
    maxMCIms[["members"]][names(topMCI)]
  }
)
CTS <- CTS[vapply(CTS, length, integer(1)) >= params$min_module_size]

state_names <- names(membersL[["MCI"]])
n_states <- length(state_names)
ncol <- min(5, n_states)
nrow <- ceiling(n_states / ncol)
grDevices::cairo_pdf(file.path(output_dir, "fig_case_mci_panels.pdf"), width = 2.4 * ncol, height = 2.2 * nrow)
par(mfrow = c(nrow, ncol), mar = c(2, 2.5, 2, 0.5), oma = c(0, 0, 2, 0))
for (st in state_names) {
  mci_vec <- membersL[["MCI"]][[st]]
  if (length(mci_vec) == 0) {
    plot.new()
    title(main = paste0("Cluster ", st, "\n0 modules"))
    next
  }
  cols <- ifelse(mci_vec >= params$mci_threshold, "#E45756", "#9ECAE1")
  barplot(
    mci_vec,
    col = cols,
    las = 2,
    cex.names = 0.45,
    main = paste0("Cluster ", st, "\n", length(mci_vec), " modules"),
    ylab = "MCI",
    ylim = c(0, max(5, max(mci_vec, na.rm = TRUE) * 1.1))
  )
  abline(h = params$mci_threshold, lty = 2, col = "grey40")
}
mtext(
  "BioTIP gastrulation case study (reduced scope: top 10 clusters)",
  side = 3, outer = TRUE, line = 0, cex = 0.9
)
dev.off()

ic_scores <- lapply(CTS, function(genes) {
  getIc(logmat, samplesL, genes, fun = "BioTIP", PCC_sample.target = "average")
})

delta_ic <- vapply(ic_scores, function(x) {
  y <- sort(unlist(x), decreasing = TRUE)
  if (length(y) < 2) return(NA_real_)
  unname(y[1] - y[2])
}, numeric(1))

set.seed(params$seed + 1)
perm_dt <- rbindlist(lapply(names(CTS), function(id) {
  genes <- CTS[[id]]
  observed <- delta_ic[[id]]
  random_delta <- replicate(params$ic_permutations, {
    random_genes <- sample(rownames(logmat), length(genes))
    score <- getIc(logmat, samplesL, random_genes, fun = "BioTIP", PCC_sample.target = "average")
    y <- sort(unlist(score), decreasing = TRUE)
    if (length(y) < 2) return(NA_real_)
    unname(y[1] - y[2])
  })
  data.table(
    cts_id = id,
    n_genes = length(genes),
    observed_delta_ic = observed,
    permutation_p = mean(random_delta >= observed, na.rm = TRUE)
  )
}))
fwrite(perm_dt, file.path(output_dir, "case_study_cts_summary.tsv"), sep = "\t")

n_cts <- max(1L, length(CTS))
cts_names <- if (length(CTS) > 0) names(CTS) else "none"
grDevices::cairo_pdf(
  file.path(output_dir, "fig_case_ic_permutation.pdf"),
  width = 10,
  height = 2.8 * n_cts
)
if (length(CTS) == 0) {
  plot.new()
  title("No CTS modules passed filters in reduced run")
} else {
  par(mfrow = c(n_cts, 2), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))
  for (id in names(CTS)) {
    ic_vec <- unlist(ic_scores[[id]])
    states <- names(ic_vec)
    plot(
      ic_vec, type = "b", pch = 16, col = "#E45756",
      xlab = "Cluster/state", ylab = "Ic.shrink",
      main = paste0("CTS ", id, " (", length(CTS[[id]]), " genes)"),
      xaxt = "n"
    )
    axis(1, at = seq_along(states), labels = states, las = 2, cex.axis = 0.6)
    observed <- delta_ic[[id]]
    random_delta <- replicate(params$ic_permutations, {
      random_genes <- sample(rownames(logmat), length(CTS[[id]]))
      score <- getIc(logmat, samplesL, random_genes, fun = "BioTIP", PCC_sample.target = "average")
      y <- sort(unlist(score), decreasing = TRUE)
      if (length(y) < 2) return(NA_real_)
      unname(y[1] - y[2])
    })
    hist(
      random_delta, col = "grey85", border = "white",
      main = paste0("Permutation (p=", format(mean(random_delta >= observed, na.rm = TRUE), digits = 3), ")"),
      xlab = expression(Delta * Ic), breaks = 20
    )
    abline(v = observed, col = "#E45756", lwd = 2, lty = 2)
  }
}
dev.off()

writeLines(
  c(
    paste("run_time:", Sys.time()),
    "mode: memory_light",
    paste("n_clusters:", length(samplesL)),
    paste("n_cts:", length(CTS)),
    paste("genes:", nrow(logmat)),
    paste("cells:", ncol(logmat)),
    "outputs:",
    "  fig_case_mci_panels.pdf",
    "  fig_case_ic_permutation.pdf",
    "  case_study_cts_summary.tsv"
  ),
  file.path(log_dir, "case_study_run.log")
)

message("Case study figures saved to ", output_dir)
