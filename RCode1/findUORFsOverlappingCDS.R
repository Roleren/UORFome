#arcsRFU = commandArgs(trailingOnly = T)

source("./createfasta.R")
library(GenomicFeatures)
library(GenomicAlignments)




# If called without path, need rangesofUORFs in global scope!!
#' Find all uorf into the cds,filter bad ones
filterORFs <- function(rangesOfuORFs, loadPath = NULL, saveToFile = F,outputFastaAndBed = F, nameSave = NULL){
  
  print("starting to filter out ourfs...")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)
  
  #check start upstream, ending downstream
  rangesOfuORFs = removeORFsWithinCDS(rangesOfuORFs)
  
  rangesOfuORFs <- ORFik:::sortPerGroup(rangesOfuORFs)
  print("finished filtering ourfs")
  
  ################SAVING###############
  ####Save overlaps of cds and uorfs for plotting later
  
  if(saveToFile){ #Not working!
    cat("Saving uorfs rdata to: ",loadPath)
    save(rangesOfuORFs, file = loadPath)
  }
  if(outputFastaAndBed){
    nameU = chooseUORFName(loadPath,nameSave)
    createFastaAndBedFile(rangesOfuORFs,loadPath = NULL,nameUsed = nameU)
  }
  return(rangesOfuORFs)
}
### Use findOverlaps to find equal start sites, this work for - strand ?
removeORFsWithinCDS <- function(grl){
  getCDS()

  overlaps <- findOverlaps(query = grl, cds, type = "within")
  grl <- grl[-unique(from(overlaps))]
  
  print("Removed uorfs that were within cds'")
  return(grl)
}
