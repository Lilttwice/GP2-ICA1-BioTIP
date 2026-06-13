# Step 3 fix: re-download only the large HighResolution count matrix (~76 MB)
# Run after the main download failed with "length != reported length"

options(timeout = 3600)
if (capabilities("libcurl")) {
  options(download.file.method = "libcurl")
  options(download.file.extra = "--retry 8 --retry-delay 10 --connect-timeout 120")
}

setwd("/Users/liutongtong/Desktop/GP2 ICA1 5.22")

destfile <- "01_data/raw/GSE141259/GSE141259_HighResolution_rawcounts.mtx.gz"
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_HighResolution_rawcounts.mtx.gz"
expected <- 76106578

if (file.exists(destfile)) {
  message("Current size: ", file.info(destfile)$size, " (need ", expected, ")")
  unlink(destfile)
}

message("Downloading (~76 MB). This can take 5–20 minutes. Do not stop R.")
download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)

size <- file.info(destfile)$size
message("Downloaded size: ", size)
if (size != expected) {
  stop("File still incomplete. Try campus VPN or a stable network, then run again.")
}

con <- gzfile(destfile, open = "rb")
readBin(con, "raw", n = 1)
close(con)
message("Gzip OK. Now run: source('CODE/M030_gse141259_preprocessing/03_code/preprocess_gse141259.R')")
