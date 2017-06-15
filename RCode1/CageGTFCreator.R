arcs = commandArgs(trailingOnly = T)
source("./HelperVariables.R")
source("./CageDataIntegration.R")
source("./HelperFunctions.R")
source("./GenomicGetters.R")

exportNewCageLeader = function(cageName = standardCage, outputFolder = leadersFolder,createUORFS = F){
  #change 5' utr sizes
  cat("Using cage:\n",cageName)
  
  getLeaders() #get leaders and make new ones from cage
  print("finished new fivePRIMe UTRs")
  
  if(exists("generalName") == F){
    source("./uorfomeGeneratorHelperFunctions.R")
    makeGeneralName(cageName)
  }
  leaderNameFull = paste0(generalName,".Leader.bed",sep = "")#for full name
  newUORFName = paste0(generalName,".RangesUorf.rdata",sep = "")
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

