rm(list=ls())
#.libPaths(c("/Home/ii/hakontj/R/x86_64-redhat-linux-gnu-library/3.3","/usr/lib64/R/library","/usr/share/R/library" ))
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./uorfomeGeneratorHelperFunctions.R")
library(doParallel)

maxCores = as.integer(detectCores()/2)
cl <- makeCluster(maxCores)
registerDoParallel(cl)
uorfsList = list.files(uorfFolder)
nuorfsList = length(uorfsList)
#clusterCall(cl, function(x) .libPaths(x), .libPaths())

foreach(i=1:nuorfsList) %dopar% {
  source("./uorfomeGeneratorHelperFunctions.R")
  load(p(uorfFolder, list.files(uorfFolder)[i]))

  uorfID <- ORFik:::orfID(rangesOfuORFs)
  saveName = paste0(resultsFolder,"/uorfIDs/",gsub("uorf.rdata","",list.files(uorfFolder)[i]),"uorfID.rdata")
  save(uorfID,file = saveName)
}

stopCluster(cl)
