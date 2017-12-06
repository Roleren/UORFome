#5. Find uorfs from sequences in fastafile combined with positions from utrs
##this version will lose the unbound open reading frames
library(ORFik)
source("./HelperLibraries.R")
source("./findUORFsOverlappingCDS.R")

#' Get ranges of uorfs, fiveUTRs must be GRangesList
#' Save to file set to false by default
scanUORFs = function(fiveUTRs,saveToFile = T,outputName = NULL, assignUorf = T){
  
  cat("started scanning for uorfs\n")
  
  rangesOfuORFs = getUnfilteredUORFs(fiveUTRs,assignUorf)
  if(is.null(outputName)){
    rangesOfuORFs = filterORFs(rangesOfuORFs, outputFastaAndBed = T)
  }else{
    rangesOfuORFs = filterORFs(rangesOfuORFs, outputFastaAndBed = T, nameSave = outputName)
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
