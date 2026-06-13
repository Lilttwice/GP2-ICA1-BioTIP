# Fix Analysis B only (Step 6 produced empty CTS) — looser filters + replot
setwd("/Users/liutongtong/Desktop/GP2 ICA1 5.22")
suppressPackageStartupMessages({
  library(BioTIP)
  library(Seurat)
  library(data.table)
  library(ggplot2)
})
try(mem.maxVSize(48 * 1024^3), silent = TRUE)

fig_dir <- "figures_for_word"
params_b <- list(
  min_state_cells = 15,
  max_cells_per_state = 100,
  max_genes = 1800,
  preselect_cutoff = 0.1,
  preselect_percent = 0.8,
  preselect_B = 3,
  network_fdr = 0.2,
  min_module_size = 10,
  n_state_candidate = 5,
  ic_permutations = 100,
  seed = 20260526
)

output_dir <- "02_modules/M060_analysisB_cellstate_biotip/04_output"
seu <- readRDS("01_data/processed/final_tuned_dataset.rds")
counts <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "data"))
state <- as.character(seu$cell_type)

samplesL_all <- split(colnames(seu), state)
samplesL_all <- samplesL_all[lengths(samplesL_all) >= params_b$min_state_cells]

set.seed(params_b$seed)
samplesL <- lapply(samplesL_all, function(cells) {
  if (length(cells) > params_b$max_cells_per_state) sample(cells, params_b$max_cells_per_state) else cells
})

selected <- unlist(samplesL, use.names = FALSE)
counts <- counts[, selected, drop = FALSE]
gv <- apply(counts, 1, var)
genes <- names(sort(gv, decreasing = TRUE))[seq_len(min(params_b$max_genes, length(gv)))]
counts <- counts[genes, , drop = FALSE]
message("Analysis B matrix: ", nrow(counts), " x ", ncol(counts), " (", length(samplesL), " cell types)")
gc()

subcounts <- optimize.sd_selection(
  counts, samplesL,
  cutoff = params_b$preselect_cutoff,
  percent = params_b$preselect_percent,
  B = params_b$preselect_B,
  doParallel = FALSE,
  n_cores = 1
)
igraphL <- getNetwork(subcounts, fdr = params_b$network_fdr)
cluster <- getCluster_methods(igraphL)
cluster <- cluster[vapply(cluster, inherits, logical(1), "communities")]
subcounts <- subcounts[names(cluster)]
samplesL <- samplesL[names(cluster)]

membersL <- getMCI(cluster, subcounts, adjust.size = FALSE, fun = "BioTIP")
topMCI <- getTopMCI(
  membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]],
  min = params_b$min_module_size,
  n = min(params_b$n_state_candidate, length(membersL[["MCI"]]))
)
maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], minsize = params_b$min_module_size, n = 3)
maxMCI <- getMaxStats(membersL[["MCI"]], maxMCIms[["idx"]])
CTS <- tryCatch(
  getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)]),
  error = function(e) maxMCIms[["members"]][names(topMCI)]
)
CTS <- CTS[vapply(CTS, length, integer(1)) >= params_b$min_module_size]

if (length(CTS) == 0) {
  stop("Analysis B still has 0 CTS. Use analysisB_simplified_results.rds fallback in Word text.")
}

ic_scores <- lapply(CTS, function(genes) {
  getIc(counts, samplesL, genes, fun = "BioTIP", PCC_sample.target = "average")
})
delta_ic <- vapply(ic_scores, function(x) {
  y <- sort(unlist(x), decreasing = TRUE)
  if (length(y) < 2) NA_real_ else unname(y[1] - y[2])
}, numeric(1))

set.seed(params_b$seed + 1)
perm_summary <- rbindlist(lapply(names(CTS), function(id) {
  genes <- CTS[[id]]
  observed <- delta_ic[[id]]
  random_delta <- replicate(params_b$ic_permutations, {
    rg <- sample(rownames(counts), length(genes))
    sc <- getIc(counts, samplesL, rg, fun = "BioTIP", PCC_sample.target = "average")
    y <- sort(unlist(sc), decreasing = TRUE)
    if (length(y) < 2) NA_real_ else unname(y[1] - y[2])
  })
  data.table(
    analysis = "analysisB",
    cts_id = id,
    n_genes = length(genes),
    observed_delta_ic = observed,
    permutation_p = mean(random_delta >= observed, na.rm = TRUE),
    best_state = names(which.max(unlist(ic_scores[[id]])))
  )
}))
fwrite(perm_summary, file.path(output_dir, "analysisB_cts_summary.tsv"), sep = "\t")
saveRDS(
  list(counts = counts, samplesL = samplesL, membersL = membersL, CTS = CTS,
       ic_scores = ic_scores, delta_ic = delta_ic),
  file.path(output_dir, "analysisB_results.rds")
)

ctsB <- names(CTS)[which.max(delta_ic)]
ic_vec <- unlist(ic_scores[[ctsB]])
dt <- data.table(state = names(ic_vec), ic = ic_vec)[order(-ic)]
p1 <- ggplot(dt, aes(reorder(state, -ic), ic)) +
  geom_col(fill = "#2F6F73", width = 0.7) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Cell type", y = "Ic.shrink", title = paste0("Analysis B — ", ctsB))
ggsave(file.path(fig_dir, "fig_analysisB_ic_trajectory.pdf"), p1, width = 7, height = 4.5, device = cairo_pdf)

genesB <- CTS[[ctsB]]
set.seed(params_b$seed + 2)
null_vals <- replicate(params_b$ic_permutations, {
  rg <- sample(rownames(counts), length(genesB))
  sc <- getIc(counts, samplesL, rg, fun = "BioTIP", PCC_sample.target = "average")
  y <- sort(unlist(sc), decreasing = TRUE)
  if (length(y) < 2) NA_real_ else unname(y[1] - y[2])
})
p2 <- ggplot(data.table(delta = null_vals), aes(delta)) +
  geom_histogram(fill = "grey85", color = "white", bins = 25) +
  geom_vline(xintercept = delta_ic[[ctsB]], color = "#E45756", linewidth = 1, lty = "dashed") +
  theme_minimal() +
  labs(title = paste0("Analysis B permutation (", ctsB, ")"), x = expression(Delta * Ic), y = "Count")
ggsave(file.path(fig_dir, "fig_analysisB_permutation.pdf"), p2, width = 6, height = 4, device = cairo_pdf)

dt_mci <- rbindlist(lapply(names(membersL[["MCI"]]), function(st) {
  mci <- membersL[["MCI"]][[st]]
  if (length(mci) == 0) return(NULL)
  data.table(state = st, module = names(mci), mci = as.numeric(mci))
}))
p3 <- ggplot(dt_mci, aes(module, mci, fill = state)) +
  geom_col(position = "dodge") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Analysis B — MCI modules per cell type", x = "Module", y = "MCI")
ggsave(file.path(fig_dir, "fig_analysisB_mci.pdf"), p3, width = 9, height = 4.5, device = cairo_pdf)

ctsA <- fread("02_modules/M040_analysisA_timepoint_biotip/04_output/analysisA_cts_summary.tsv")
combined <- rbind(ctsA, perm_summary, fill = TRUE)
fwrite(combined, file.path(fig_dir, "combined_AB_cts_summary.tsv"), sep = "\t")
combined[, label := fifelse(analysis == "analysisA", paste0("A: ", cts_id), paste0("B: ", cts_id))]
combined[, sig := permutation_p <= 0.05]
p4 <- ggplot(combined, aes(label, observed_delta_ic, fill = sig)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("p=", format(permutation_p, digits = 2))), hjust = -0.05, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2F6F73", "FALSE" = "#B9C0C7")) +
  theme_minimal() +
  labs(x = NULL, y = expression(Delta * Ic[shrink]), title = "CTS candidates — Analyses A and B")
ggsave(file.path(fig_dir, "fig02_biotip_cts_results.pdf"), p4, width = 7.5, height = 4.5, device = cairo_pdf)

message("Analysis B fixed. CTS: ", paste(names(CTS), collapse = ", "))
print(list.files(fig_dir, pattern = "analysisB|fig02", ignore.case = TRUE))
