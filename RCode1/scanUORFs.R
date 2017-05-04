#5. Find uorfs from sequences in fastafile combined with positions from utrs
##this version will lose the unbound open reading frames
library(ORFik)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
###Get ranges of uorfs, fiveUTRs must be GRangesList
### Save to file set to false by default
scanUORFs = function(fiveUTRs,saveToFile = F,outputName = NULL,fastaName = "/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa"){
  ###FIX LISTING###
  print("loading fasta files and unlisting fiveUTRs")
  unlistfiveUTRs = unlist(fiveUTRs) ##Unlist from GRangeslist to Granges
  names(unlistfiveUTRs) = names(fiveUTRs)
  if(exists("fasta") == F){
    fasta =  readDNAStringSet(fastaName) ##Get fasta file
    assign("fasta",fasta,envir = .GlobalEnv)
  }
  faiFileNotexists = !findFF("fai",T)## Get fasta indexed file
  if(faiFileNotexists)indexFa(fastaName)
  #if(faiFileNotexists){indexFa(findFF("fa"))}
  if(exists("fa") == F){
    fa = FaFile(fastaName)
    #fa = FaFile(findFF("fa"))
    assign("fa",fa,envir = .GlobalEnv)
  }
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
  ORFdef
  grangesObj
  ORF = map_to_GRanges(ORFdef,unlist(grangesObj),transcriptName)
}