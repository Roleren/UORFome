#5. Find uorfs from sequences in fastafile combined with positions from utrs
##this version will lose the unbound open reading frames
library(ORFik)
source("./HelperLibraries.R")
source("./findUORFsOverlappingCDS.R")
#source("./GenomicGetters.R")
###Get ranges of uorfs, fiveUTRs must be GRangesList
### Save to file set to false by default
scanUORFs = function(fiveUTRs,saveToFile = T,outputName = NULL, assignUorf = T){
  
  cat("started scanning for uorfs\n")
  
  rangesOfuORFs = getUnfilteredUORFsFast(fiveUTRs,assignUorf)
  if(is.null(outputName)){
    rangesOfuORFs = removeFalseUORFs(rangesOfuORFs, outputFastaAndBed = T)
  }else{
    rangesOfuORFs = removeFalseUORFs(rangesOfuORFs, outputFastaAndBed = T, nameSave = outputName)
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



###For each 5'utr, find the uorfs, NB! Needs global fiveUTRs
findInFrameUORF = function(i){
  grangesObj = fiveUTRs[i]
  
  ORFdef <- find_in_frame_ORFs(as.character(unlist(seqs[i])),
                               longestORF = F,
                               minimumLength = 8) #function for finding orfs
  
  if(length(ORFdef) > 0 ){ #if found any orfs
    return( map_to_GRanges(ORFdef,unlist(grangesObj),names(grangesObj)) ) # map it
  }
  return ( NULL ) #else return nothing
}

