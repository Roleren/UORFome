#5. Find uorfs from sequences in fastafile combined with positions from utrs

library(ORFik)
library(GenomicFeatures)
library(GenomicAlignments)
scanUORFs = function(fiveUTRs,saveToFile = F){
  ###FIX LISTING###
  
  unlistfiveUTRs = unlist(fiveUTRs)
  names(unlistfiveUTRs) = names(fiveUTRs)
  faiFileNotexists = !findFF("fai",T)
  if(faiFileNotexists){indexFa(findFF("fa"))}
  
  fa = FaFile(findFF("fa"))
  fiveUTRseqs = getSeq(fa,unlistfiveUTRs)
  
  cat("finished seqs")
  rangesOfuORFs = list()
  for(i in 1:length(fiveUTRs)){
    grangesObj = fiveUTRs[[i]]
    transcriptName = names(fiveUTRs)[i]
    dna <- as.character(unlist(getSeq(fa, grangesObj)))
    
    
    ORFdef <- find_in_frame_ORFs(dna,
                                 longestORF = F,
                                 minimumLength = 8)
    
    ORFcount = length(ORFdef)
    if(ORFcount > 0 ){
      ORF = map_to_GRanges(ORFdef,grangesObj,transcriptName )
    }
    else{
      ORF = NULL
    }
    rangesOfuORFs = c(rangesOfuORFs,ORF)
  }
  rangesOfuORFs = GRangesList(rangesOfuORFs)
  
  if(saveToFile){
    save(rangesOfuORFs, file = "rangesOfuorfs.rdata")
    assign("rangesOfuORFs",rangesOfuORFs,envir = .GlobalEnv)
    }
  print("rangesOfuORFs finished")
  return(rangesOfuORFs)
}
