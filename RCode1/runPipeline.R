# Parts:
#1. set up parameters
#2. Find new cage leaders
#3. Find uORFs
#4. Create feature database



##### First set up parameters
# 1. If output dir is not created yet, run:

# orfikDirs(mainPath = , makeDatabase = )

#set working dir correctly to ./RCode1/ location
setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path
source("./DataBaseSetup.R")
setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path


# set up multithreading options
pipelineCluster(10) #!! set number of cores, I use 75 usually on furu.


##### Second Find new cage leaders


cageList =  grep(pattern = ".bed", list.files(cageFolder), value = TRUE)# make sure this folder is correct
nCageList = length(cageList)

getLeadersFromCage(nCageList)

##### Third Find uORFs

leadersList = list.files(leadersFolder)
nLeadersList = length(leadersList)
#clusterCall(cl, function(x) .libPaths(x), .libPaths())

getUorfsFromLeaders(nLeadersList)

#### Fourth make uorf IDs

uorfsList = list.files(uorfFolder)
nuorfsList = length(uorfsList)
getIDsFromUorfs(nuorfsList)

### Fifth Create feature database
#5. 
source("./DataBaseSetup.R")
createCatalogueDB()

# stop cluster
stopCluster(cl)
rm(cl)


# foreach(i=1:20) %dopar% {
#   setwd("/export/valenfs/projects/uORFome/RCode1/")
#   source("./uorfomeGeneratorHelperFunctions.R")
#   leadersList = list.files(leadersFolder)
#   load(p(leadersFolder,leadersList[i]))
#   print(mean(ORFik:::widthPerGroup(fiveUTRs)))
# }