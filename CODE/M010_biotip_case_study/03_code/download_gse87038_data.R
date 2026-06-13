dir.create("01_data/raw/GSE87038", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(GEOquery))

retry_download <- function(url, destfile, attempts = 3, sleep_seconds = 2) {
  for (i in seq_len(attempts)) {
    ok <- tryCatch({
      download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
      file.exists(destfile) && file.info(destfile)$size > 0
    }, error = function(e) {
      message("Download failed on attempt ", i, ": ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) {
      return(TRUE)
    }
    Sys.sleep(sleep_seconds)
  }
  FALSE
}

files <- getGEOSuppFiles(
  GEO = "GSE87038",
  baseDir = "01_data/raw",
  fetch_files = FALSE
)

targets <- subset(files, grepl("xlsx$|RAW.tar$", fname))
target_dir <- "01_data/raw/GSE87038"

status <- lapply(seq_len(nrow(targets)), function(i) {
  fname <- targets$fname[i]
  url <- targets$url[i]
  destfile <- file.path(target_dir, fname)
  existing_ok <- file.exists(destfile) && file.info(destfile)$size > 0
  downloaded <- existing_ok || retry_download(url, destfile)
  data.frame(
    fname = fname,
    url = url,
    destfile = destfile,
    status = if (downloaded) "ok" else "failed",
    size = if (file.exists(destfile)) file.info(destfile)$size else NA_real_
  )
})

manifest <- do.call(rbind, status)

write.table(
  manifest,
  file = "CODE/M010_biotip_case_study/05_logs/download_manifest.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

overall <- if (all(manifest$status == "ok")) "completed" else "partial"
writeLines(
  c(
    paste("download_status:", overall),
    paste("download_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "geo_accession: GSE87038"
  ),
  "CODE/M010_biotip_case_study/05_logs/run.log"
)

if (!all(manifest$status == "ok")) {
  stop("Some GSE87038 files failed to download.")
}

message("GSE87038 supplementary files downloaded.")
