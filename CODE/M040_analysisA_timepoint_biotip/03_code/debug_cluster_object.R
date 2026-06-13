clusters <- readRDS("CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_clusters.rds")

lines <- c(
  paste("top_class:", class(clusters), collapse = ", "),
  paste("top_length:", length(clusters)),
  paste("top_names:", paste(names(clusters), collapse = ", "))
)

for (nm in names(clusters)) {
  obj <- clusters[[nm]]
  lines <- c(
    lines,
    "",
    paste("name:", nm),
    paste("class:", paste(class(obj), collapse = ", ")),
    paste("typeof:", typeof(obj)),
    if (is.list(obj)) paste("names:", paste(names(obj), collapse = ", ")) else "names: <not list>",
    if (!is.null(dim(obj))) paste("dim:", paste(dim(obj), collapse = "x")) else "dim: <none>",
    if (!is.null(length(obj))) paste("length:", length(obj)) else "length: <none>"
  )
}

writeLines(lines, "CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_cluster_object_debug.txt")
message("Cluster object debug saved.")
