arcs = commandArgs(trailingOnly = T)

source("/export/valenfs/projects/uORFome/RCode1/CageDataIntegration.R")
source("/export/valenfs/projects/uORFome/RCode1/HelperFunctions.R")

exportNewCageLeader = function(data = "/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data", cageName = "/export/valenfs/projects/uORFome/DATA/CAGE/human/brain%2c%20adult%2c%20donor1.CNhs11796.10084-102B3.hg38.nobarcode.ctss.bed.gz",outputFolder = "/export/valenfs/projects/uORFome/test_results/New_Cage_Leaders",createUORFS = F){
  #change 5' utr sizes
  print("Using cage")
  print(cageName)
  setwd(data)
  Gtf = makeTxDbFromGFF("/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf")
  assign("Gtf",Gtf,envir = .GlobalEnv)
  fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
  print("finding new fivePRIMe UTRs")
  fiveUTRs = getNewfivePrimeUTRs(fiveUTRs,dataName = cageName)
  print("finished new fivePRIMe UTRs")
  leaderName = gsub(".*/", "", cageName)
  leaderName = strsplit(leaderName,".hg38*.")[[1]][1]#used for uorf too
  leaderNameFull = paste0(leaderName,".Leader.bed",sep = "")#for full name
  newUORFName = paste0(leaderName,".RangesUorf.rdata",sep = "")
  setwd(outputFolder)
  export.bed(unlist(fiveUTRs),leaderNameFull)
  print("leader creation finished")
  if(createUORFS)
    createUORFs(fiveUTRs,uorfName = newUORFName)
}



if(length(arcs) > 0){
  # update uorfomegenerator, to include premade leaders
  if(length(arcs) == 1){
    exportNewCageLeader(cageName = normalizePath(arcs[1]))
  }else if(length(arcs) == 2)
    exportNewCageLeader(cageName = normalizePath(arcs[1]),createUORFS = as.logical(arcs[2]))
}
print("cage leader script finished")

