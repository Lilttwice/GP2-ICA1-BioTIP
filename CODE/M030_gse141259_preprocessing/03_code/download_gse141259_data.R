dir.create("01_data/raw/GSE141259", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(GEOquery))

# RStudio default timeout (60s) is too short for ~76MB matrix files.
options(timeout = max(3600, getOption("timeout", 60)))
if (capabilities("libcurl")) {
  options(download.file.method = "libcurl")
  options(download.file.extra = "--retry 5 --retry-delay 5 --connect-timeout 60")
}

is_valid_gzip <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 100) {
    return(FALSE)
  }
  out <- tryCatch({
    con <- gzfile(path, open = "rb")
    on.exit(close(con), add = TRUE)
    readBin(con, "raw", n = 16)
    TRUE
  }, error = function(e) FALSE)
  isTRUE(out)
}

remote_content_length <- function(url) {
  hdr <- tryCatch(
    curl::get_headers(url, redirect = TRUE),
    error = function(e) NULL
  )
  if (is.null(hdr)) {
    return(NA_real_)
  }
  line <- hdr[grep("^content-length:", tolower(hdr))]
  if (length(line) == 0) {
    return(NA_real_)
  }
  as.numeric(sub(".*:\\s*", "", line[[1]]))
}

file_is_complete <- function(destfile, url) {
  if (!file.exists(destfile)) {
    return(FALSE)
  }
  local_size <- file.info(destfile)$size
  expected <- remote_content_length(url)
  if (!is.na(expected) && expected > 0 && local_size != expected) {
    return(FALSE)
  }
  if (grepl("\\.gz$", destfile, ignore.case = TRUE)) {
    return(is_valid_gzip(destfile))
  }
  local_size > 0
}

retry_download <- function(url, destfile, attempts = 6, sleep_seconds = 5) {
  expected <- remote_content_length(url)
  if (!is.na(expected) && expected > 0) {
    message("Expected size for ", basename(destfile), ": ", expected, " bytes")
  }

  for (i in seq_len(attempts)) {
    if (file.exists(destfile) && !file_is_complete(destfile, url)) {
      message("Removing incomplete file: ", destfile)
      unlink(destfile)
    }
    if (file_is_complete(destfile, url)) {
      message("Already complete: ", basename(destfile))
      return(TRUE)
    }

    ok <- tryCatch({
      download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
      file_is_complete(destfile, url)
    }, error = function(e) {
      message("Download failed on attempt ", i, ": ", conditionMessage(e))
      if (file.exists(destfile)) {
        unlink(destfile)
      }
      FALSE
    })

    if (isTRUE(ok)) {
      message("Downloaded OK: ", basename(destfile), " (", file.info(destfile)$size, " bytes)")
      return(TRUE)
    }
    Sys.sleep(sleep_seconds * i)
  }
  FALSE
}

if (!requireNamespace("curl", quietly = TRUE)) {
  install.packages("curl", repos = "https://cloud.r-project.org")
}

files <- getGEOSuppFiles(
  GEO = "GSE141259",
  baseDir = "01_data/raw",
  fetch_files = FALSE
)

target_dir <- "01_data/raw/GSE141259"

status <- lapply(seq_len(nrow(files)), function(i) {
  fname <- files$fname[i]
  url <- files$url[i]
  destfile <- file.path(target_dir, fname)
  downloaded <- retry_download(url, destfile)
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
  file = "CODE/M030_gse141259_preprocessing/05_logs/download_manifest.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

overall <- if (all(manifest$status == "ok")) "completed" else "partial"
writeLines(
  c(
    paste("download_status:", overall),
    paste("download_time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "geo_accession: GSE141259"
  ),
  "CODE/M030_gse141259_preprocessing/05_logs/run.log"
)

failed <- manifest$fname[manifest$status != "ok"]
if (length(failed) > 0) {
  stop(
    "Some GSE141259 files failed to download: ",
    paste(failed, collapse = ", "),
    "\nDelete partial files in 01_data/raw/GSE141259/ and re-run this script."
  )
}

message("GSE141259 supplementary files downloaded.")
