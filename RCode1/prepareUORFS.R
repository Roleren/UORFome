arcs = commandArgs(trailingOnly = T)
###watch out where you call this from, arcs are scarry
source("/export/valenfs/projects/uORFome/RCode1/HelperFunctions.R")
source("/export/valenfs/projects/uORFome/RCode1/scanUORFs.R")
require(Biostrings)
library(rtracklayer)

createUORFs = function(fiveUTRs = NULL,leaderBED = NULL,uorfName = NULL){
  print("starting on new uorfs")
  
  if(is.null(uorfName)){
    leaderName = gsub(".*/", "", leaderBED)
    leaderName = strsplit(leaderName,".Leader.bed")[[1]][1]#used for uorf too
    uorfName = paste0(leaderName,".RangesUorf.rdata",sep = "")
    cat("new name is: ", uorfName)
  }
  if(!is.null(leaderBED)){
    cat("retrieving 5utrs from bed file: ", uorfName)
    fiveUTRs = import.bed(leaderBED)
    fiveUTRs = as.data.frame(fiveUTRs)
    fiveUTRstest <- split(fiveUTRs, fiveUTRs$name)
    fiveUTRstest = lapply(fiveUTRstest,function(x) as(x,"GRanges"))
    
    fiveUTRstest = GRangesList(fiveUTRstest)
    fiveUTRs = fiveUTRstest
  }
  setwd("/export/valenfs/projects/uORFome/test_results/rangesOfUORFs/")
  assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
  scanUORFs(fiveUTRs,saveToFile = T,outputName = uorfName)
}
if(length(arcs) == 1)
  createUORFs(leaderBED = normalizePath(arcs[1]))

#createUORFs(leaderBED = "/export/valenfs/projects/uORFome/test_results/New_Cage_Leaders/CD34%2b%20stem%20cells%20-%20adult%20bone%20marrow%20derived%2c%20donor1%2c%20tech_rep1.CNhs12588.12225-129F2.Leader.bed")
print("finished uorfs preparation")