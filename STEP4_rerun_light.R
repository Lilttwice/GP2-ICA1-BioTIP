# Step 4 retry after memory error — does NOT reinstall packages
setwd("/Users/liutongtong/Desktop/GP2 ICA1 5.22")

eh_cache <- path.expand("~/Library/Caches/org.R-project.R/R/ExperimentHub")
if (!dir.exists(eh_cache)) {
  dir.create(eh_cache, recursive = TRUE, showWarnings = FALSE)
}

message("Restarting case study with reduced memory settings...")
gc()
source("02_modules/M010_biotip_case_study/03_code/run_gse87038_case_study_simplified.R")

cat("\n=== Check outputs ===\n")
list.files(
  "02_modules/M010_biotip_case_study/04_output",
  pattern = "fig_case|case_study",
  full.names = TRUE
)
