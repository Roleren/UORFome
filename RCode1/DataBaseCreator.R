
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./CreateCatalogueHelpers.R")
source("./databaseHelpers.R")

createUORFAtlas <- function(){
  
  j = 1
  for(i in idFiles){
    load(p(idFolder, i))
    
    uorfs <- uorfID[!duplicated(uorfID)]
    if(j == 1) {
      uorfsTotal <- data.table(cbind(uorfs,rep(T,length(uorfs))))
      colnames(uorfsTotal) = c("uorfID", "1")
    }else{
      uorfs <- data.table(cbind(uorfs,rep(T,length(uorfs))))
      colnames(uorfs) = c("uorfID", toString(j))
      uorfsTotal = merge(uorfsTotal,uorfs,by = "uorfID", all = T)
      colnames(uorfsTotal) = c("uorfID", as.character(1:(ncol(uorfsTotal)-1)))
    }
    j = j+1
  }
  uorfIDs <- uorfsTotal
  save(uorfIDs,file = "UORFAtlas.rdata")
  insertTable(Matrix = uorfIDs,tableName = "uorfAtlas")
  # for(i in 2:ncol(tesTotal)){
  #   tesTotal[ is.na(tesTotal[,i, with = F]), i] = F 
  # }
}

createUniqueIDs <- function(){
  j = 1
  for(i in idFiles){
    load(p(idFolder, i))
    
    tes = as.data.table(uorfID[!duplicated(uorfID)])
    colnames(tes) = "uorfID"
    if(j == 1) {
      allUniqueIDs = tes
    }else{
      allUniqueIDs = merge(tes,allUniqueIDs, by = "uorfID", all = T)
    }
    j = j+1
  }
  save(allUniqueIDs,file = "allUniqueIDs.rdata")
  insertTable(Matrix = allUniqueIDs,tableName = "uniqueIDs")
}



getTissueTable = function(){
  require(xlsx)
  cageTable = read.xlsx("../HumanSamples2.0.sdrf.xlsx", sheetName = "Sheet1")
  
  cageTable = as.data.table(cageTable)
  cageTable[is.na(Characteristics.Tissue.)] = "unclassifiable"
  
  insertTable(Matrix = cageTable,tableName = "cageInformation")
  
  uniqueTissues = unique(cageTable$Characteristics.Tissue.)  
  
  matchCageIDandCageName = rep("a",length(cageFiles))
  j = 1
  for(i in uorfFiles){
    matchCageIDandCageName[j] = gsub(".*\\.", "", gsub(".hg38.*","",gsub(".*CNhs","",i)))
    j = j + 1
  }
  
  matchCageIDandCageName = as.data.table(matchCageIDandCageName)
  colnames(matchCageIDandCageName) = colnames(cageTable)[1]
  matchCageIDandCageName$cage_index = 1:nrow(matchCageIDandCageName)
  cageWeHave <-  merge(cageTable, matchCageIDandCageName,by = "Source.Name")
  
  if(sum(duplicated(cageWeHave$cage_index)) != 0) stop("duplicated indexes used, check input")
  
  uorfIDs <- readTable("uniqueIDs")
  uorfIDTable <- readTable("uorfAtlas")
  
  finalMatrix <- as.data.table(matrix(nrow = nrow(uorfIDs), ncol = length(uniqueTissues)+1))
  finalMatrix[,1] <- uorfIDs
  colnames(finalMatrix)[1] <- "uorfID"
  uniqueTissues <- as.character(uniqueTissues)
  colnames(finalMatrix)[2:ncol(finalMatrix)] <- uniqueTissues
  
  for(i in 2:(length(uniqueTissues)+1)){
    cageFilestoCheck <- cageWeHave[Characteristics.Tissue. == uniqueTissues[i-1]]$cage_index
    cageFilestoCheck <- cageFilestoCheck + 1
    makeGroupingForColumn <- rowSums(!is.na(tesTotal[,cageFilestoCheck, with = F]))

    finalMatrix[,i] <- makeGroupingForColumn > 1
  }
  #fix duplicates, merge them
  test <- finalMatrix$urethra | finalMatrix$Urethra
  finalMatrix$urethra <- test
  finalMatrix$Urethra <- NULL
  
  save(finalMatrix,file = "tissueAtlas.rdata")
  insertTable(Matrix = finalMatrix,tableName = "tissueAtlasByCage")
}



rfpTables <- function(){
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./MatchExperimentsHeader.R")
  source("./uorfomeGeneratorHelperFunctions.R")
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  
  grl <- uniqueIdsAsGR()
  validateExperiments(grl)
  # now make ribo and rna seq tables
  SpeciesGroup <- getUnfilteredSpeciesGroups()
  rpfFilePaths <- getFilteredRFPPaths(SpeciesGroup)
  
  riboTable <- riboAtlasFPKMAll(grl, rpfFilePaths)
  
  riboAtlasFPKMTissue(grl,rpfFilesPaths,riboTables,SpeciesGroup)
  
  
}


#createCatalogueDB(databaseName,matrix,tableName)
#0st name of experiments
###1st cage supprt per uorf, bool, 
###2nd rna ribo seq supprt and te per uorf, matrix
###3rd 

###example: find uorf that are only in certain tissue
uorfDB = createDataBase(databaseName)
#createUniqueIDs()
#createUORFAtlas()
#rfpAndrnaSupportTable()

createCatalogueDB = function(name, bigMatrix, tableName){
  uorfDB = createDataBase(databaseName)
  createUniqueIDs()
  createUORFAtlas()
  rfpTables()
}