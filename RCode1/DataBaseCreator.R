
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./uorfomeGeneratorHelperFunctions.R")
source("./databaseHelpers.R")
source("./DataBaseGetters.R")
source("./DataBaseInfo.R")
source("./DataBaseValidation.R")
source("./TissueTables.R")
source("./CreateCatalogueHelpers.R")


#' Create 1 column of all unique ids from uorfID folder
createUniqueIDs <- function(){
  j = 1
  for(i in idFiles){
    load(p(idFolder, i))
    uorfID <- unique(uorfID)
    if (j == 1) {
      allUniqueIDs <- uorfID
    }else{
      matching <- uorfID %in% allUniqueIDs
      toAdd <- uorfID[which(matching == F)]
      allUniqueIDs <- c(allUniqueIDs,toAdd)
    }
    j <- j+1
  }
  allUniqueIDs <- sort(allUniqueIDs)
  save(allUniqueIDs,file = "allUniqueIDs.rdata")
  insertTable(Matrix = allUniqueIDs,tableName = "uniqueIDs")
}

createUORFAtlas <- function(){
  uorfIDsAllUnique <- readTable("uniqueIDs")
  colnames(uorfIDsAllUnique) = "uorfID"
  uorfAtlas <- as.data.table(matrix(F, nrow = nrow(uorfIDsAllUnique), ncol = length(idFiles)+1))
  uorfAtlas[,1] <- uorfIDsAllUnique
  colnames(uorfAtlas) = c("uorfID", as.character(1:(length(idFiles))))
  j = 1
  for(i in idFiles){
    load(p(idFolder, i))
    
    uorfs <- uorfID[!duplicated(uorfID)]
    
    uorfAtlas[,(j+1)] <- uorfIDsAllUnique$uorfID %in% uorfs
      
    j = j+1
  }
  
  save(uorfAtlas,file = "UORFAtlas.rdata")
}

#' Create tissueTable for cage, 1 row per unique uorf 
getTissueTable <- function(){
  
  cageTable <- getCageInfoTable()
  uniqueTissues <- unique(cageTable$Characteristics.Tissue.)  
  
  cageWeHave <- getAllUsableCage(cageTable)
  
  # load needed tables, and make tissue atlas of cage
  load("UORFAtlas.rdata")
  uorfIDs <- readTable("uniqueIDs")
  colnames(uorfIDs) = "uorfID"
  uorfAtlasRows0 <- which(rowSums(uorfAtlas[,2:ncol(uorfAtlas)]) == 0)
  if(length(uorfAtlasRows0) > 0) stop("uorfAtlas and unique uorf IDs does not match!")
  
  finalMatrix <- as.data.table(matrix(nrow = nrow(uorfIDs), ncol = length(uniqueTissues)+1))
  finalMatrix[,1] <- uorfIDs$uorfID
  colnames(finalMatrix)[1] <- "uorfID"
  uniqueTissues <- as.character(uniqueTissues)
  colnames(finalMatrix)[2:ncol(finalMatrix)] <- uniqueTissues
  
  for(i in 2:(length(uniqueTissues)+1)){
    cageFilestoCheck <- cageWeHave[Characteristics.Tissue. == uniqueTissues[i-1]]$cage_index
    cageFilestoCheck <- cageFilestoCheck + 1
    makeGroupingForColumn <- rowSums(uorfAtlas[,cageFilestoCheck, with = F])

    finalMatrix[,i] <- makeGroupingForColumn > 1
  }
  tissueAtlas <- finalMatrix
  
  save(tissueAtlas,file = "tissueAtlas.rdata")
  insertTable(Matrix = tissueAtlas,tableName = "tissueAtlasByCage")
}



# rfpTables <- function(){
#   setwd("/export/valenfs/projects/uORFome/RCode1/")
#   source("./MatchExperimentsHeader.R")
#   setwd("/export/valenfs/projects/uORFome/dataBase/")
#   
#   grl <- uniqueIdsAsGR()
#   validateExperiments(grl)
#   # now make ribo and rna seq tables
#   SpeciesGroup <- getUnfilteredSpeciesGroups()
#   rpfFilePaths <- getFilteredRFPPaths(SpeciesGroup)
#   insertTable(rpfFilePaths, "RiboSeqInfo")
#   riboTable <- riboAtlasFPKMAll(grl, rpfFilePaths)
#   
#   riboAtlasFPKMTissue(grl,rpfFilesPaths,riboTable,SpeciesGroup)
#   
# }


#createCatalogueDB(databaseName,matrix,tableName)
#0st name of experiments
###1st cage supprt per uorf, bool, 
###2nd rna ribo seq supprt and te per uorf, matrix
###3rd 

###example: find uorf that are only in certain tissue
uorfDB <- createDataBase(databaseName)
#createUniqueIDs()
#createUORFAtlas()
#rfpTables()

# fix the presentation
# Take all uorfs, and cluster the tissue based on uorfs
# can you recover the cell type based on uorf usage
# hclust
clusterByCageTissue <- function(){
  tissueAtlas <- readTable("tissueAtlasByCage")
  binDist <- dist(x = t(tissueAtlas[1:1000, 2:4]), method = "binary")
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  save(binDist,file = "binDist.rdata")
  hclustResult <-  hclust(binDist)
  save(hclustResult,file = "hclustResult.rdata")
}
#clusterByCageTissue()