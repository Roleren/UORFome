#' Full report
#' @param experiment an ORFik experiment data.frame
#' @return NULL (all output stored to disc and .GlobalEnv)
ORFikReport <- function(experiment, outputDir = ".", type = c("transcript", "cds")) {
  validateExperiments(experiment)
  
  txdbB <- loadTxdb("txdbPath")
  
  computeFeatures()
}

