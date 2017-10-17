source("./uorfomeGeneratorHelperFunctions.R")
library(doParallel)

maxCores = as.integer(detectCores()-(detectCores()/2))
cl <- makeCluster(maxCores)
registerDoParallel(cl)
cageList = list.files(cageFolder)
nCageList = length(cageList)

foreach(i=1:nCageList) %dopar% {
  
  source("./uorfomeGeneratorHelperFunctions.R")
  cageList = list.files(cageFolder)
  getLeaders(usingNewCage = T, cageName = p(cageFolder,cageList[i]) , assignLeader = F )

}


stopCluster(cl)

