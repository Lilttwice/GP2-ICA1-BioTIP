suppressPackageStartupMessages(library(BioTIP))

functions_to_check <- c(
  "optimize.sd_selection",
  "getNetwork",
  "getCluster_methods",
  "getMCI",
  "getTopMCI",
  "getMaxMCImember",
  "getMaxStats",
  "getCTS",
  "getIc",
  "simulation_Ic"
)

lines <- unlist(lapply(functions_to_check, function(fn) {
  c(
    paste("FUNCTION", fn),
    capture.output(args(get(fn))),
    ""
  )
}))

writeLines(lines, "CODE/M040_analysisA_timepoint_biotip/04_output/biotip_api_signatures.txt")
message("BioTIP API signatures saved.")
