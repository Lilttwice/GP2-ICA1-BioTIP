suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

base_dir <- "01_data/raw/GSE141259"

read_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con)
}

highres_barcodes <- read_lines(file.path(base_dir, "GSE141259_HighResolution_barcodes.txt.gz"))
highres_genes <- read_lines(file.path(base_dir, "GSE141259_HighResolution_genes.txt.gz"))
whole_barcodes <- read_lines(file.path(base_dir, "GSE141259_WholeLung_barcodes.txt.gz"))
whole_genes <- read_lines(file.path(base_dir, "GSE141259_WholeLung_genes.txt.gz"))

highres_cellinfo <- fread(file.path(base_dir, "GSE141259_HighResolution_cellinfo.csv.gz"))
whole_cellinfo <- fread(file.path(base_dir, "GSE141259_WholeLung_cellinfo.csv.gz"))

summary_lines <- c(
  paste("inspect_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("highres_barcodes_n:", length(highres_barcodes)),
  paste("highres_genes_n:", length(highres_genes)),
  paste("whole_barcodes_n:", length(whole_barcodes)),
  paste("whole_genes_n:", length(whole_genes)),
  paste("highres_cellinfo_rows:", nrow(highres_cellinfo)),
  paste("highres_cellinfo_cols:", ncol(highres_cellinfo)),
  paste("whole_cellinfo_rows:", nrow(whole_cellinfo)),
  paste("whole_cellinfo_cols:", ncol(whole_cellinfo)),
  "",
  "highres_cellinfo_colnames:",
  paste(colnames(highres_cellinfo), collapse = ", "),
  "",
  "whole_cellinfo_colnames:",
  paste(colnames(whole_cellinfo), collapse = ", ")
)

writeLines(summary_lines, "CODE/M030_gse141259_preprocessing/04_output/input_structure_summary.txt")

fwrite(
  data.table(dataset = "highres", column = colnames(highres_cellinfo)),
  "CODE/M030_gse141259_preprocessing/04_output/highres_cellinfo_columns.tsv",
  sep = "\t"
)

fwrite(
  data.table(dataset = "wholelung", column = colnames(whole_cellinfo)),
  "CODE/M030_gse141259_preprocessing/04_output/wholelung_cellinfo_columns.tsv",
  sep = "\t"
)

capture_unique_values <- function(dt, dataset_name) {
  candidate_cols <- grep(
    "time|day|cell|type|cluster|anno|state|sample|injury|treat|condition",
    colnames(dt),
    ignore.case = TRUE,
    value = TRUE
  )
  out <- lapply(candidate_cols, function(col) {
    vals <- unique(dt[[col]])
    vals <- vals[!is.na(vals)]
    vals <- as.character(vals)
    vals <- vals[seq_len(min(length(vals), 30))]
    data.table(dataset = dataset_name, column = col, example_values = paste(vals, collapse = " | "))
  })
  if (length(out) == 0) {
    return(data.table(dataset = character(), column = character(), example_values = character()))
  }
  rbindlist(out, fill = TRUE)
}

examples <- rbindlist(
  list(
    capture_unique_values(highres_cellinfo, "highres"),
    capture_unique_values(whole_cellinfo, "wholelung")
  ),
  fill = TRUE
)

fwrite(
  examples,
  "CODE/M030_gse141259_preprocessing/04_output/cellinfo_candidate_fields.tsv",
  sep = "\t"
)

writeLines("GSE141259 input inspection completed.", "CODE/M030_gse141259_preprocessing/05_logs/run.log")
message("Input inspection completed.")
