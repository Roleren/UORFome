# Parts:
#1. set up parameters
#2. Find new cage leaders
#3. Find uORFs
#4. Create feature database



##### First set up parameters
# 1. If output dir is not created yet, run:

# orfikDirs(mainPath = , makeDatabase = )

#set working dir correctly to ./RCode1/ location
setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path as codeFolder
# source("./HelperFunctions.R")
# updateORFik("update19", "Roleren") # update if needed
source("./pipelineSetup.R")


# set up multithreading options
pipelineCluster(35) #!! set number of cores, I use 69 usually on furu.


##### Second Find new cage leaders

getLeadersFromCage(length(cageFiles))

##### Third Find uORFs

getUorfsFromLeaders()

#### Fourth make uorf IDs

getIDsFromUorfs()


### Fifth Create feature database
#5. 
source(p(codeFolder,"/DataBaseSetup.R"))
createCatalogueDB()

# stop cluster, not needed in prediction, it uses h2o package
stopCluster(cl)
rm(cl)


# Sixth Predict uORFs
# Either split tissues, or use argument "all" for all in one go
#6.
predictUorfs("all")


# r <- as.data.table(RFP, stringsAsFactors = FALSE)
# g <- as.data.table(grl, stringsAsFactors = FALSE)
# 
# setkeyv(r, c("seqnames", "strand", "start", "end"))
# setkeyv(g, c("seqnames", "strand", "start", "end"))
# data.table::foverlaps(r, g, type = "within")
