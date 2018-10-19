
# If called without path, need rangesofUORFs in global scope!!
#' Find all uorf into the cds,filter bad ones
filterORFs <- function(rangesOfuORFs, loadPath = NULL, saveToFile = F,outputFastaAndBed = F, nameSave = NULL){
  
  print("starting to filter out ourfs...")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)
  
  #check start upstream, ending downstream
  rangesOfuORFs = removeORFsWithinCDS(rangesOfuORFs)
  
  rangesOfuORFs <- ORFik:::sortPerGroup(rangesOfuORFs)
  rangesOfuORFs <- removeORFsWithinSameStop(rangesOfuORFs)
  rangesOfuORFs <- removeORFsWithinSameStart(rangesOfuORFs)
  print("finished filtering ourfs")

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


removeORFsWithinSameStop <- function(grl){
  getCDS()
  
  overlaps <- findOverlaps(query =  stopSites(grl, asGR = TRUE),
                           stopSites(cds, asGR = TRUE), type = "within")
  grl <- grl[-unique(from(overlaps))]
  
  print("Removed uorfs that had same stop as cds'")
  return(grl)
}

removeORFsWithinSameStart <- function(grl){
  getCDS()
  
  # filter out uORFs with same start as cds
  starts <- startSites(grl, asGR = T, is.sorted = T)
  cdsstarts <- startSites(cds, asGR = T, is.sorted = T)
  overlaps <- findOverlaps(starts, cdsstarts, type = "within")
  grl <- grl[-unique(from(overlaps))]
  
  print("Removed uorfs that had same start as cds'")
  return(grl)
}
