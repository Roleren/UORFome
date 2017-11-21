source("./uorfomeGeneratorHelperFunctions.R")
library(foreach)
#run these only once
resultsFolder = p(mainFolder,"/zebrafish_results")
helperMainFolder = p(resultsFolder,"/Old_Tests/test_data")
dataFolder = helperMainFolder
fastaName = p(dataFolder,"/danRer10.fa")
faiName = p(dataFolder,"/danRer10.fa")
gtfName = p(dataFolder,"/Danio_rerio.GRCz10.83.gtf")
leadersbedFolder = p(resultsFolder,"/New_Cage_bedLeaders/")
leadersFolder = p(resultsFolder,"/New_Cage_Leaders/")
fastaFolder = p(resultsFolder,"/fasta/")
uorfBedFolder = p(resultsFolder,"/bedUORFS/")
uorfFolder = p(resultsFolder,"/rangesOfUORFs/")

Gtf = makeTxDbFromGFF(gtfName)
assign("Gtf",Gtf,envir = .GlobalEnv)

#change cage name for each and run rest once per cage
cageName = "/export/valenfs/data/processed_data/CAGE/nepal_2013_zebrafish/called_peaks/leaders_zv10_zf_04_cells_512.csv"
fiveUTRs = otherSpeciesCageCSV(cageName)

exportNamebed = paste0(leadersbedFolder,getRelativePathName(p(cageName,".leader.bed")))
exportNamerdata = paste0(leadersFolder,getRelativePathName(p(cageName,".leader.rdata")))
print("exporting new leaders")
export.bed(unlist(fiveUTRs),exportNamebed)
save(fiveUTRs,file = exportNamerdata)

#now cage are finished, run next to find uorfs, once per cage

leadersList = list.files(leadersFolder)
nLeadersList = length(leadersList)

foreach(i=1:nLeadersList) %do% {
  
  usedCage = gsub(pattern = ".leader.rdata",replacement = "",x = leadersList[i]) 
  
  if(UorfRangesNotExists(assignUorf =  F,givenCage = usedCage)){
    load(p(leadersFolder,leadersList[i]))
    scanUORFs(fiveUTRs, outputName = usedCage, assignUorf = F)
  }else{
    i
  }
}