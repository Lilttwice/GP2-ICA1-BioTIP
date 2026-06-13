suppressPackageStartupMessages({
  library(data.table)
})

seu <- readRDS("01_data/processed/final_tuned_dataset.rds")
subcounts <- readRDS("CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_subcounts.rds")
clusters <- readRDS("CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_clusters.rds")

state_sizes <- as.data.table(table(seu[["time_point"]]))
setnames(state_sizes, c("time_point", "n_cells"))
state_sizes <- state_sizes[order(n_cells)]

subcount_sizes <- data.table(
  state = names(subcounts),
  n_genes = vapply(subcounts, nrow, integer(1)),
  n_cells = vapply(subcounts, ncol, integer(1))
)

cluster_sizes <- rbindlist(lapply(names(clusters), function(nm) {
  obj <- clusters[[nm]]
  if (!inherits(obj, "communities")) {
    return(data.table(
      state = nm,
      n_nodes = NA_integer_,
      n_communities = NA_integer_,
      status = class(obj)[1]
    ))
  }
  data.table(
    state = nm,
    n_nodes = as.integer(obj$vcount),
    n_communities = length(unique(obj$membership)),
    status = "communities"
  )
}))

write.table(
  state_sizes,
  file = "CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_state_sizes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  subcount_sizes,
  file = "CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_subcount_sizes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cluster_sizes,
  file = "CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_cluster_sizes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_lines <- c(
  paste("inspect_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("states_n:", nrow(state_sizes)),
  paste("smallest_state:", state_sizes$time_point[1], "=", state_sizes$n_cells[1]),
  paste("largest_state:", state_sizes$time_point[.N], "=", state_sizes$n_cells[.N])
)

writeLines(summary_lines, "CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_object_summary.txt")
message("Analysis A object inspection completed.")
