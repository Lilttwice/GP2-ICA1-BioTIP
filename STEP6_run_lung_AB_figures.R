# Step 6 — Re-run lung Analyses A & B (one at a time) + export A-grade figures
# Safe on 32 GB Mac: run in RStudio; expect ~30–90 min total.

setwd("/Users/liutongtong/Desktop/GP2 ICA1 5.22")
suppressPackageStartupMessages({
  library(BioTIP)
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

try(mem.maxVSize(48 * 1024^3), silent = TRUE)

fig_dir <- "figures_for_word"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  min_state_cells = 20,
  max_cells_per_state = 80,
  max_genes = 1500,
  preselect_cutoff = 0.1,
  preselect_percent = 0.8,
  preselect_B = 3,
  network_fdr = 0.2,
  min_module_size = 20,
  n_state_candidate = 4,
  ic_permutations = 100,
  seed = 20260525
)

run_biotip_AB <- function(state_col, analysis_tag, out_mod) {
  output_dir <- file.path(out_mod, "04_output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  seu <- readRDS("01_data/processed/final_tuned_dataset.rds")
  counts <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "data"))
  state <- as.character(seu[[state_col, drop = TRUE]])

  samplesL_all <- split(colnames(seu), state)
  samplesL_all <- samplesL_all[order(names(samplesL_all))]
  samplesL_all <- samplesL_all[lengths(samplesL_all) >= params$min_state_cells]

  set.seed(params$seed)
  samplesL <- lapply(samplesL_all, function(cells) {
    if (length(cells) > params$max_cells_per_state) sample(cells, params$max_cells_per_state) else cells
  })

  selected <- unlist(samplesL, use.names = FALSE)
  counts <- counts[, selected, drop = FALSE]
  gv <- apply(counts, 1, var)
  genes <- names(sort(gv, decreasing = TRUE))[seq_len(min(params$max_genes, length(gv)))]
  counts <- counts[genes, , drop = FALSE]

  message(analysis_tag, ": ", length(samplesL), " states; matrix ", nrow(counts), " x ", ncol(counts))
  gc()

  subcounts <- optimize.sd_selection(
    counts, samplesL,
    cutoff = params$preselect_cutoff,
    percent = params$preselect_percent,
    B = params$preselect_B,
    doParallel = FALSE,
    n_cores = 1
  )
  igraphL <- getNetwork(subcounts, fdr = params$network_fdr)
  cluster <- getCluster_methods(igraphL)
  cluster <- cluster[vapply(cluster, inherits, logical(1), "communities")]
  subcounts <- subcounts[names(cluster)]
  samplesL <- samplesL[names(cluster)]
  rm(igraphL)
  gc()

  membersL <- getMCI(cluster, subcounts, adjust.size = FALSE, fun = "BioTIP")
  topMCI <- getTopMCI(
    membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]],
    min = params$min_module_size,
    n = min(params$n_state_candidate, length(membersL[["MCI"]]))
  )
  maxMCIms <- getMaxMCImember(membersL[["members"]], membersL[["MCI"]], minsize = params$min_module_size, n = 3)
  maxMCI <- getMaxStats(membersL[["MCI"]], maxMCIms[["idx"]])
  CTS <- tryCatch(
    getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)]),
    error = function(e) maxMCIms[["members"]][names(topMCI)]
  )
  CTS <- CTS[vapply(CTS, length, integer(1)) >= params$min_module_size]

  ic_scores <- lapply(CTS, function(genes) {
    getIc(counts, samplesL, genes, fun = "BioTIP", PCC_sample.target = "average")
  })
  delta_ic <- vapply(ic_scores, function(x) {
    y <- sort(unlist(x), decreasing = TRUE)
    if (length(y) < 2) NA_real_ else unname(y[1] - y[2])
  }, numeric(1))

  set.seed(params$seed + 1)
  perm_summary <- rbindlist(lapply(names(CTS), function(id) {
    genes <- CTS[[id]]
    observed <- delta_ic[[id]]
    random_delta <- replicate(params$ic_permutations, {
      rg <- sample(rownames(counts), length(genes))
      sc <- getIc(counts, samplesL, rg, fun = "BioTIP", PCC_sample.target = "average")
      y <- sort(unlist(sc), decreasing = TRUE)
      if (length(y) < 2) NA_real_ else unname(y[1] - y[2])
    })
    data.table(
      analysis = analysis_tag,
      cts_id = id,
      n_genes = length(genes),
      observed_delta_ic = observed,
      permutation_p = mean(random_delta >= observed, na.rm = TRUE),
      best_state = names(which.max(unlist(ic_scores[[id]])))
    )
  }))
  fwrite(perm_summary, file.path(output_dir, paste0(analysis_tag, "_cts_summary.tsv")), sep = "\t")

  saveRDS(
    list(counts = counts, samplesL = samplesL, membersL = membersL, CTS = CTS,
         ic_scores = ic_scores, delta_ic = delta_ic, perm = perm_summary),
    file.path(output_dir, paste0(analysis_tag, "_results.rds"))
  )

  invisible(list(membersL = membersL, CTS = CTS, ic_scores = ic_scores,
                 delta_ic = delta_ic, perm_summary = perm_summary, samplesL = samplesL))
}

plot_ic_trajectory <- function(ic_scores, cts_id, outfile, title_prefix) {
  ic_vec <- unlist(ic_scores[[cts_id]])
  ord <- order(ic_vec, decreasing = TRUE)
  dt <- data.table(state = names(ic_vec)[ord], ic = ic_vec[ord])
  dt[, state := factor(state, levels = state)]
  p <- ggplot(dt, aes(state, ic)) +
    geom_col(fill = "#2F6F73", width = 0.7) +
    geom_point(color = "#E45756", size = 2) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "State", y = "Ic.shrink", title = paste0(title_prefix, " — ", cts_id))
  ggsave(outfile, p, width = 8, height = 4.5, device = cairo_pdf)
}

plot_permutation <- function(counts, samplesL, genes, observed_delta, outfile, title) {
  set.seed(params$seed + 2)
  null_vals <- replicate(params$ic_permutations, {
    rg <- sample(rownames(counts), length(genes))
    sc <- getIc(counts, samplesL, rg, fun = "BioTIP", PCC_sample.target = "average")
    y <- sort(unlist(sc), decreasing = TRUE)
    if (length(y) < 2) NA_real_ else unname(y[1] - y[2])
  })
  p <- ggplot(data.table(delta = null_vals), aes(delta)) +
    geom_histogram(fill = "grey85", color = "white", bins = 25) +
    geom_vline(xintercept = observed_delta, color = "#E45756", linewidth = 1, linetype = "dashed") +
    theme_minimal(base_size = 11) +
    labs(
      x = expression(Delta * Ic[shrink]),
      y = "Count",
      title = paste0(title, " (p = ", format(mean(null_vals >= observed_delta, na.rm = TRUE), digits = 3), ")")
    )
  ggsave(outfile, p, width = 6, height = 4, device = cairo_pdf)
}

plot_mci_celltype <- function(membersL, outfile) {
  dt <- rbindlist(lapply(names(membersL[["MCI"]]), function(st) {
    mci <- membersL[["MCI"]][[st]]
    if (length(mci) == 0) return(NULL)
    data.table(state = st, module = names(mci), mci = as.numeric(mci))
  }))
  if (is.null(dt) || nrow(dt) == 0) return(invisible(NULL))
  top <- dt[, .(max_mci = max(mci)), by = state][order(-max_mci)]
  dt <- dt[state %in% top$state]
  p <- ggplot(dt, aes(module, mci, fill = state)) +
    geom_col(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "Module", y = "MCI", title = "Analysis B — MCI by cell type (top modules)")
  ggsave(outfile, p, width = 9, height = 4.5, device = cairo_pdf)
}

message("=== Analysis A (time point) ===")
resA <- run_biotip_AB("time_point", "analysisA", "02_modules/M040_analysisA_timepoint_biotip")
rm(resA)
gc()

message("=== Analysis B (cell type) ===")
resB <- run_biotip_AB("cell_type", "analysisB", "02_modules/M060_analysisB_cellstate_biotip")
rm(resB)
gc()

# --- Export figures ---
resA_obj <- readRDS("02_modules/M040_analysisA_timepoint_biotip/04_output/analysisA_results.rds")
resB_obj <- readRDS("02_modules/M060_analysisB_cellstate_biotip/04_output/analysisB_results.rds")

ctsA <- names(resA_obj$CTS)[1]
if (!is.na(ctsA)) {
  plot_ic_trajectory(
    resA_obj$ic_scores, ctsA,
    file.path(fig_dir, "fig_analysisA_ic_trajectory.pdf"),
    "Analysis A"
  )
  genesA <- resA_obj$CTS[[ctsA]]
  plot_permutation(
    resA_obj$counts, resA_obj$samplesL, genesA,
    resA_obj$delta_ic[[ctsA]],
    file.path(fig_dir, "fig_analysisA_permutation.pdf"),
    paste0("Analysis A permutation (", ctsA, ")")
  )
}

ctsB <- names(resB_obj$CTS)[1]
if (!is.na(ctsB)) {
  plot_ic_trajectory(
    resB_obj$ic_scores, ctsB,
    file.path(fig_dir, "fig_analysisB_ic_trajectory.pdf"),
    "Analysis B"
  )
  genesB <- resB_obj$CTS[[ctsB]]
  plot_permutation(
    resB_obj$counts, resB_obj$samplesL, genesB,
    resB_obj$delta_ic[[ctsB]],
    file.path(fig_dir, "fig_analysisB_permutation.pdf"),
    paste0("Analysis B permutation (", ctsB, ")")
  )
}
plot_mci_celltype(resB_obj$membersL, file.path(fig_dir, "fig_analysisB_mci.pdf"))

# Update combined CTS bar chart
cts <- fread("02_modules/M040_analysisA_timepoint_biotip/04_output/analysisA_cts_summary.tsv")
cts2 <- fread("02_modules/M060_analysisB_cellstate_biotip/04_output/analysisB_cts_summary.tsv")
combined <- rbind(cts, cts2, fill = TRUE)
fwrite(combined, file.path(fig_dir, "combined_AB_cts_summary.tsv"), sep = "\t")

combined[, label := fifelse(analysis == "analysisA", paste0("A: ", cts_id), paste0("B: ", cts_id))]
combined[, sig := permutation_p <= 0.05]
p <- ggplot(combined, aes(label, observed_delta_ic, fill = sig)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("p=", format(permutation_p, digits = 2))), hjust = -0.05, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2F6F73", "FALSE" = "#B9C0C7")) +
  theme_minimal() +
  labs(x = NULL, y = expression(Delta * Ic[shrink]), title = "CTS candidates — Analyses A and B")
ggsave(file.path(fig_dir, "fig02_biotip_cts_results.pdf"), p, width = 7.5, height = 4, device = cairo_pdf)

message("=== Step 6 done. New PDFs in figures_for_word/ ===")
print(list.files(fig_dir, pattern = "\\.pdf$"))
