setwd("/export/valenfs/projects/uORFome/RCode1/")
library(doParallel)

maxCores <- as.integer(detectCores()-(detectCores()/2))
cl <- makeCluster(maxCores)
registerDoParallel(cl)
# load RNA list
load("/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/unfilteredSpeciesGroup.rdata")
rnaList <- SpeciesGroup[SpeciesGroup$Sample_Type == "RNA",]$RnaRfpFolders
rnaList <- grep(pattern = ".bam",x = rnaList, value = T)
nrnaList <- length(rnaList)
rm(SpeciesGroup)
# all rfp features
rnaFPKMs <- foreach(i=1:nrnaList, .combine = 'cbind') %dopar% {
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./DataBaseCreator.R")
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  load("/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/unfilteredSpeciesGroup.rdata")
  rnaList <- SpeciesGroup[SpeciesGroup$Sample_Type == "RNA",]$RnaRfpFolders
  RNAPath <- rnaList[i]
  grlNames <- getORFNamesDB(T, F, F)
  
  # get tx
  tx <- ORFik:::extendLeaders(getTx())
  tx <- tx[grlNames]
  RNA <- readGAlignments(RNAPath)
  rnaFPKM <- ORFik:::fpkm(tx, reads = RNA)
}
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./DataBaseCreator.R")
insertTable(rnaFPKMs, "rnaFPKMs")
insertTable(rnaList, "RNA-seqInfo")

stopCluster(cl)
