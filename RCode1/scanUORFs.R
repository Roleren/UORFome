#5. Find uorfs from sequences in fastafile combined with positions from utrs
##this version will lose the unbound open reading frames
library(ORFik)
source("/export/valenfs/projects/uORFome/RCode1/HelperLibraries.R")
source("/export/valenfs/projects/uORFome/RCode1/GenomicGetters.R")

###Get ranges of uorfs, fiveUTRs must be GRangesList
### Save to file set to false by default
scanUORFs = function(fiveUTRs,saveToFile = F,outputName = NULL){
  ###FIX LISTING###
  print("loading fasta files and unlisting fiveUTRs")
  unlistfiveUTRs = unlist(fiveUTRs) ##Unlist from GRangeslist to Granges
  names(unlistfiveUTRs) = names(fiveUTRs)
  
  getFasta() #get fasta and fai
  
  cat("finished fiveUTR-seqs\n")
  cat("finding ranges of uorfs, this takes around 4 hours.\n")
  rangesOfuORFs = lapply(X = 1:length(fiveUTRs), FUN = findInFrameUORF)
  
  
  rangesOfuORFs = GRangesList(unlist(rangesOfuORFs))
  rangesOfuORFs = rangesOfuORFs[width(rangesOfuORFs) > 0]
  
  rangesOfuORFs = removeFalseUORFs(NULL)
  
  assign("rangesOfuORFs",rangesOfuORFs,envir = .GlobalEnv)
  
  if(saveToFile){
    print("saving rangesOfuORFs")
    if(!is.null(outputName))
      save(rangesOfuORFs, file = outputName)
    else
      save(rangesOfuORFs, file = "rangesOfuorfs.rdata")
  }
  print("rangesOfuORFs finished")
  return(rangesOfuORFs)
}
###For each 5'utr, find the uorfs, NB! Needs global fiveUTRs
findInFrameUORF = function(i){
  grangesObj = fiveUTRs[i]
  transcriptName = names(grangesObj)
  dna <- as.character(unlist(getSeq(fa, unlist(grangesObj))))
  dna <- paste(dna,collapse = "")
  
  ORFdef <- find_in_frame_ORFs(dna,
                               longestORF = F,
                               minimumLength = 8)
  
  if(length(ORFdef) > 0 ) map_granges(ORFdef,grangesObj,transcriptName) else NULL
  
}

#Use map_to_granges from the ORFik package by Cornel -CBU
map_granges = function(ORFdef,grangesObj,transcriptName){
  names(grangesObj[[1]]) <- rep(transcriptName, length(unlist(grangesObj)))
  
  ORFranges <- GRanges(transcriptName, ORFdef)
  
  ORF = map_to_GRanges(ORFdef,unlist(grangesObj),transcriptName)
}