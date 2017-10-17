#5. Find uorfs from sequences in fastafile combined with positions from utrs
##this version will lose the unbound open reading frames
library(ORFik)
source("./HelperLibraries.R")
source("./findUORFsOverlappingCDS.R")
#source("./GenomicGetters.R")
###Get ranges of uorfs, fiveUTRs must be GRangesList
### Save to file set to false by default
scanUORFs = function(fiveUTRs,saveToFile = T,outputName = NULL, assignUorf = T){
  ###FIX LISTING###
  print(thisCage)
  #unlistfiveUTRs = unlist(fiveUTRs) ##Unlist from GRangeslist to Granges
  #names(unlistfiveUTRs) = names(fiveUTRs) What happend to this part ????
  
  getFasta() #get fasta and fai
  
  cat("finished fiveUTR-seqs\n")
  cat("finding ranges of uorfs, this takes around 4 hours.\n")
  getUnfilteredUORFs(fiveUTRs)
  rangesOfuORFs = removeFalseUORFs(outputFastaAndBed = T) #add possibilities here!
  
  
  if(saveToFile){
    print("saving rangesOfuORFs")
    if(!is.null(outputName))
      save(rangesOfuORFs, file = outputName)
    else
      save(rangesOfuORFs, file = paste0(uorfFolder,getRelativePathName(thisCage),".uorf.rdata" ))
  }
  print("rangesOfuORFs finished")
  return(rangesOfuORFs)
}



###For each 5'utr, find the uorfs, NB! Needs global fiveUTRs
findInFrameUORF = function(i){
  grangesObj = fiveUTRs[i]
  
  dna <- as.character(unlist(getSeq(fa, unlist(grangesObj))))
  dna <- paste(dna,collapse = "")
  
  ORFdef <- find_in_frame_ORFs(dna,
                               longestORF = F,
                               minimumLength = 8) #function for finding orfs
  
  if(length(ORFdef) > 0 ){ #if found any orfs
    transcriptName = names(grangesObj) #get name of tx  
    return( map_granges(ORFdef,grangesObj,transcriptName) ) # map it
  }else return ( NULL ) #else return nothing
}

#Use map_to_granges from the ORFik package by Cornel -CBU
map_granges = function(ORFdef,grangesObj,transcriptName){
  names(grangesObj[[1]]) <- rep(transcriptName, length(unlist(grangesObj)))
  
  # ORFranges <- GRanges(transcriptName, ORFdef)
  
  ORF = map_to_GRanges(ORFdef,unlist(grangesObj),transcriptName)
}
