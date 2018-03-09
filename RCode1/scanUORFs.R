#5. Find uorfs from sequences in fastafile combined with positions from utrs
##this version will lose the unbound open reading frames
library(ORFik)
source("./HelperLibraries.R")
source("./findUORFsOverlappingCDS.R")

#' Find uorfs
#' 
#' @param fiveUTRs must be GRangesList
#' @param saveToFile set to false by default
scanUORFs = function(fiveUTRs, saveToFile = T, outputName = NULL,
                     assignUorf = T, outputFastaAndBed = T, filterORFs = T){
  
  cat("started scanning for uorfs\n")
  
  rangesOfuORFs <- getUnfilteredUORFs(fiveUTRs,assignUorf)
  gr <- unlist(rangesOfuORFs, use.names = F)
  rangesOfuORFs <- ORFik:::groupGRangesBy(gr, gr$names)
  if (filterORFs) {
    if(is.null(outputName)){
      rangesOfuORFs <- filterORFs(rangesOfuORFs, outputFastaAndBed = outputFastaAndBed)
    }else{
      rangesOfuORFs <- filterORFs(rangesOfuORFs, outputFastaAndBed = outputFastaAndBed, nameSave = outputName)
    }
  }
  
  if(saveToFile){
    print("saving rangesOfuORFs")
    if(!is.null(outputName)){
      save(rangesOfuORFs, file = getUORFRDataName(outputName))
    }else if(exists("thisCage") && !is.null(thisCage)){
      save(rangesOfuORFs, file = getUORFRDataName(thisCage))
    }else if(exists("cageName") && !is.null(cageName)){
      save(rangesOfuORFs, file = getUORFRDataName(cageName))
    }else{
      stop("no name to set RangesOfUorfs with!")
    }
  }
  print("rangesOfuORFs finished")
  return(rangesOfuORFs)
}
