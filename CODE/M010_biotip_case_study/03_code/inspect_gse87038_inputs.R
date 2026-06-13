suppressPackageStartupMessages(library(readxl))

cellinfo <- read_excel(
  "01_data/raw/GSE87038/GSE87038_Cell_Information_and_Barcode_Information.xlsx",
  sheet = "All_single_cells_information",
  skip = 2
)

barcode_info <- read_excel(
  "01_data/raw/GSE87038/GSE87038_Cell_Information_and_Barcode_Information.xlsx",
  sheet = "Barcode_information"
)

umi_head <- read_excel(
  "01_data/raw/GSE87038/GSE87038_Mouse_Organogenesis_UMI_counts_matrix.xlsx",
  sheet = 1,
  n_max = 10
)

summary_lines <- c(
  paste("inspect_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("cellinfo_rows:", nrow(cellinfo)),
  paste("cellinfo_cols:", ncol(cellinfo)),
  paste("barcode_rows:", nrow(barcode_info)),
  paste("barcode_cols:", ncol(barcode_info)),
  paste("umi_head_rows:", nrow(umi_head)),
  paste("umi_head_cols:", ncol(umi_head)),
  "",
  "cellinfo_colnames:",
  paste(names(cellinfo), collapse = ", "),
  "",
  "barcode_colnames:",
  paste(names(barcode_info), collapse = ", "),
  "",
  "umi_colnames_first_20:",
  paste(names(umi_head)[seq_len(min(length(names(umi_head)), 20))], collapse = ", ")
)

writeLines(summary_lines, "CODE/M010_biotip_case_study/04_output/input_structure_summary.txt")

write.table(
  head(cellinfo, 20),
  file = "CODE/M010_biotip_case_study/04_output/cellinfo_head.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  head(barcode_info, 20),
  file = "CODE/M010_biotip_case_study/04_output/barcode_info_head.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("GSE87038 input inspection completed.")
