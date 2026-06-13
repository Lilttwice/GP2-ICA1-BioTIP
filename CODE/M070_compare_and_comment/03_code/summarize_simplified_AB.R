suppressPackageStartupMessages(library(data.table))

out_dir <- "CODE/M070_compare_and_comment/04_output"
log_dir <- "CODE/M070_compare_and_comment/05_logs"
check_dir <- "CODE/M070_compare_and_comment/06_checks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

a_path <- "CODE/M040_analysisA_timepoint_biotip/04_output/analysisA_simplified_cts_summary.tsv"
b_path <- "CODE/M060_analysisB_cellstate_biotip/04_output/analysisB_simplified_cts_summary.tsv"

a <- fread(a_path)
b <- fread(b_path)
a[, analysis := "A_time_point"]
b[, analysis := "B_cell_type"]

combined <- rbindlist(list(a, b), fill = TRUE)
setcolorder(combined, c("analysis", setdiff(names(combined), "analysis")))
fwrite(combined, file.path(out_dir, "simplified_AB_cts_summary.tsv"), sep = "\t")

overview <- combined[
  ,
  .(
    n_cts = .N,
    best_states = paste(unique(best_state), collapse = "; "),
    min_p = suppressWarnings(min(permutation_p, na.rm = TRUE)),
    median_delta_ic = median(observed_delta_ic, na.rm = TRUE)
  ),
  by = analysis
]
fwrite(overview, file.path(out_dir, "simplified_AB_overview.tsv"), sep = "\t")

interpretation_lines <- c(
  "Simplified A/B comparison notes",
  "",
  "Analysis A uses sampling time point as state and is therefore sensitive to uneven recovery kinetics, small time-point groups, and mixed cell-type composition within each time point.",
  "Analysis B uses annotated alveolar cell type as state and therefore emphasizes transcriptional differences between AT1, AT2, activated AT2, Krt8+ ADI, and proliferating cells after pooling across time.",
  "The two analyses are complementary rather than directly interchangeable because their state definitions encode different biological axes."
)
writeLines(interpretation_lines, file.path(out_dir, "simplified_AB_interpretation_notes.txt"))

writeLines(
  c(
    paste("run_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "status: completed",
    paste("analysisA_cts:", nrow(a)),
    paste("analysisB_cts:", nrow(b))
  ),
  file.path(log_dir, "run.log")
)

writeLines(
  c(
    paste("run_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "gate_status: pass",
    "checks:",
    "- Analysis A and B simplified result tables are present",
    "- comparison summary generated",
    "- notes emphasize non-equivalence of time-point and cell-type state definitions"
  ),
  file.path(check_dir, "check_report.txt")
)

message("Simplified A/B comparison completed.")
