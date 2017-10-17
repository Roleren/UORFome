source("./uorfomeGeneratorHelperFunctions.R")
library(doParallel)

#maxCores = as.integer(detectCores()-(detectCores()/2))
maxCores = 50
cl <- makeCluster(maxCores)
registerDoParallel(cl)
leadersList = list.files(leadersFolder)
nLeadersList = length(leadersList)

foreach(i=1:nLeadersList) %dopar% {
  source("./uorfomeGeneratorHelperFunctions.R")
  leadersList = list.files(leadersFolder)
  load(p(leadersFolder,leadersList[i]))
  usedCage = gsub(pattern = ".leader.rdata",replacement = "",x = leadersList[i]) 
  decideHowToGetUORFRanges(T,usedCage)
  
}


stopCluster(cl)

