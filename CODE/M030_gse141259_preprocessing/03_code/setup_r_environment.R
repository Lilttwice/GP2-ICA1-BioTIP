required_cran <- c("BiocManager", "Seurat", "readxl", "data.table")
required_bioc <- c("BioTIP", "GEOquery", "Matrix")

install_if_missing <- function(pkgs, installer) {
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)
  if (length(missing) == 0) {
    message("All packages already installed for: ", paste(pkgs, collapse = ", "))
    return(invisible(NULL))
  }
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  installer(missing)
}

install_if_missing(
  required_cran,
  function(pkgs) install.packages(pkgs, repos = "https://cloud.r-project.org")
)

install_if_missing(
  required_bioc,
  function(pkgs) BiocManager::install(pkgs, ask = FALSE, update = FALSE)
)

session_info <- capture.output({
  cat("R.version.string:", R.version.string, "\n")
  cat("Bioconductor.version:", as.character(BiocManager::version()), "\n")
  for (pkg in c(required_cran, required_bioc)) {
    cat(pkg, ":", as.character(packageVersion(pkg)), "\n")
  }
})

writeLines(session_info, "CODE/M030_gse141259_preprocessing/05_logs/env.txt")
writeLines("BioTIP environment setup completed.", "CODE/M030_gse141259_preprocessing/05_logs/run.log")
message("Environment setup completed.")
