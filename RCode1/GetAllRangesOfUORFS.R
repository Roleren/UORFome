rm(list=ls())
#.libPaths(c("/Home/ii/hakontj/R/x86_64-redhat-linux-gnu-library/3.3","/usr/lib64/R/library","/usr/share/R/library" ))
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./uorfomeGeneratorHelperFunctions.R")
library(doParallel)

maxCores = as.integer(detectCores()/2)
cl <- makeCluster(maxCores)
registerDoParallel(cl)
leadersList = list.files(leadersFolder)
nLeadersList = length(leadersList)
#clusterCall(cl, function(x) .libPaths(x), .libPaths())

foreach(i=1:nLeadersList) %dopar% {
  source("./uorfomeGeneratorHelperFunctions.R")
  leadersList = list.files(leadersFolder)
  
  usedCage = gsub(pattern = ".leader.rdata",replacement = "",x = leadersList[i]) 
  
  if(UorfRangesNotExists(assignUorf =  F,givenCage = usedCage)){
    load(p(leadersFolder,leadersList[i]))
    scanUORFs(fiveUTRs, outputName = usedCage, assignUorf = F)
    print("ok")
  }else{
    i
  }
}

stopCluster(cl)

