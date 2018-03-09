# Parts:
#1. set up parameters
#2. Find new cage leaders
#3. Find uORFs
#4. Create feature database



##### First set up parameters
# 1. If output dir is not created yet, run:

# orfikDirs(mainPath = , makeDatabase = )

#set working dir correctly to ./RCode1/ location
#setwd()

# set up multithreading options
pipelineCluster()


##### Second Find new cage leaders

source("./uorfomeGeneratorHelperFunctions.R")

cageList = list.files(cageFolder) # make sure this folder is correct
nCageList = length(cageList)

getLeadersFromCage(nCageList)



##### Third Find uORFs

leadersList = list.files(leadersFolder)
nLeadersList = length(leadersList)
#clusterCall(cl, function(x) .libPaths(x), .libPaths())

getUorfsFromLeaders(nLeadersList)


### Fourth Create feature database
#4. 

setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./DataBaseCreator.R")

createCatalogueDB() # fix this to work

# stop cluster
stopCluster(cl)
rm(cl)