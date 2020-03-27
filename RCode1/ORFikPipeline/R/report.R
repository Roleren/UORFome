#' Full ORFik report
#' @param experiment an ORFik experiment data.frame
#' @return NULL (all output stored to disc and .GlobalEnv)
ORFikReport <- function(experiment, outputDir = ".", subset = NULL) {
  if (is(experiment, "list")) {
    dfl <- experiment
  } else dfl <- list(experiment)
    
    
  for (experiment in dfl){
    
    ORFik:::validateExperiments(experiment)
    if (is.null(experiment@txdb)) stop("experiment must have defined txdb")
    
    if (!is.null(subset)) {
      experiment <- experiment[experiment$libtype == subset,]
    }
    libTypes <- ORFik:::libraryTypes(experiment)
    outputLibs(experiment)
    varNames <- ORFik:::bamVarName(experiment)
    
    if (length(varNames) == 1) {
      if (libTypes %in% c("RFP")) {
        print("Running single Ribo-seq")
      } else if (libTypes %in% c("RNA")) {
        print("Running single RNA-seq")
      } else if (libTypes %in% c("CAGE")) {
        print("Running single CAGE")
      } else if (libTypes %in% c("LSU")) {
        print("Running single LSU")
      } else if (libTypes %in% c("SSU")) {
        print("Running single SSU")
      } 
    } else {
      variants <- libTypes
      
      if (length(variants) == 1) {
        if (libTypes %in% c("RFP")) {
          print("Running mutliple Ribo-seq")
        } else if (libTypes %in% c("RNA")) {
          print("Running mutliple RNA-seq")
        } else if (libTypes %in% c("CAGE")) {
          print("Running mutliple CAGE")
        } else if (libTypes %in% c("LSU")) {
          print("Running mutliple LSU")
        } else if (libTypes %in% c("SSU")) {
          print("Running mutliple SSU")
        } 
      } else { # mutliple variable
        if (variants %in% all("RNA", "RFP")) {
          print("Running mutliple Ribo-seq & RNA-seq")
          
          computeFeatures()
        } 
      }
    }
  }
  return(NULL)
}
