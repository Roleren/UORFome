
createUniqueIDsAtlasNonTranscript <- function(){
  stop("not working, not needed anymore!")
  A = readTable("MergedByExonsuniqueUORFs")
  c = paste(A$seqnames,A$start,A$end,A$strand)
  uniqueIDs = as.data.table(c)
  colnames(uniqueIDs) = "uorfID"
  insertTable(Matrix = uniqueIDs, tableName = "uorfIDNonTranscript")
}

#' get GrangesList from uniqueIDs as strings
uniqueIdsAsGR <- function(){
  if(tableNotExists("SplittedByExonsuniqueUORFs")){
    uniqueIDs <- readTable("uniqueIDs")
    if (sum(duplicated(uniqueIDs)) > 0 ) stop("duplicated uorf names in uniqueIDs")
    grl <- toGRFromUniqueID(uniqueIDs)
    insertTable(Matrix = grl,tableName = "SplittedByExonsuniqueUORFs", rmOld = T)
    #grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
    
    #now make uscs bed 12 format
    bed12(grl, "bedUniqueUorfs.bed", T)
    
  } else {
    grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
  }
  return(grl)
}

#' Riboseq table per experiment(ribo-seq file)
riboAtlasFPKMAll <- function(grl,rpfFilePaths){
  # get all ribo experiment, map reads, get fpkm, put in column
  riboTable <- as.data.table(matrix(nrow = length(grl), ncol = length(rpfFilePaths)+1))
  riboTable[,1] = names(grl)
  j <- 1
  for(i in rpfFilePaths[j:length(rpfFilePaths)]){
    
    riboTable[,(j+1)] <- RiboFPKM(grl, i)
    j = j + 1
  }
  # test <- riboTable[,1:j-1]
  # save(test,file = "riboFPKM.rdata")
  colnames(riboTable)[1] <- "uorfIndex"
  save(riboTable,file = "riboFPKM.rdata")
  insertTable(riboTable,"riboAll", rmOld = T)
  if ((ncol(riboTable) - 1) != length(rpfFilePaths)){
    stop("something wrong in creation of riboTable")
  }
  return(riboTable)
}

allFeaturesAtlas <- function(){
  
}
#' Riboseq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is mean fpkm per tissue
riboAtlasFPKMTissue <- function(grl,rpfFilesPaths,riboTable,SpeciesGroup){
  # now do per tissue true/false
  rpfSamples <- SpeciesGroup[SpeciesGroup$Sample_Type == "RPF",]
  #1. we have tissues in speciesGroup
  tissuesUsed <- rpfSamples[rpfSamples$RnaRfpFolders %in% rpfFilesPaths,]
  uniqueTissues <- as.character(unique(tissuesUsed$Tissue_or_CellLine))
  #2. we have all the tables in riboTables
  #3. So for each tissue, find the group, then for each group ->
  uniqueIDs <- readTable("uniqueIDs")
  riboByTissue <- as.data.table(matrix(nrow = nrow(riboTable), ncol = length(uniqueTissues)+1))
  riboByTissue[,1] <- uniqueIDs
  colnames(riboByTissue)[1] <- "uorfID"
  colnames(riboByTissue)[2:ncol(riboByTissue)] <- uniqueTissues
  
  #4. rowSum(riboColumns > 1) > 1
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf, it will be true
  for(i in uniqueTissues){
    indexes <- which(tissuesUsed$Tissue_or_CellLine == i)
    indexes <- indexes + 1
    riboColumns <- riboTable[,indexes, with=F]
    riboByTissue[,i] <- rowSums(riboColumns > 1) > 1
    
  }
  insertTable(Matrix = riboByTissue,tableName = "RiboByTissueTF", rmOld = T)
  
  #now get mean value instead of true/false
  riboByTissueMean <- riboByTissue
  rm(riboByTissue)
  #4. rowSum(riboColumns > 1) > 1
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf, it will be true
  for(i in uniqueTissues){
    indexes <- which(tissuesUsed$Tissue_or_CellLine == i)
    indexes <- indexes + 1
    riboColumns <- riboTable[,indexes, with=F]
    riboByTissueMean[,i] <- rowMeans(riboColumns)
  }
  insertTable(Matrix = riboByTissueMean,tableName = "RiboByTissueMean", rmOld = T)
}
getCageInfoTable <- function(){
  if (tableNotExists("cageInformation")){
    require(xlsx)
    cageTable <- read.xlsx("../HumanSamples2.0.sdrf.xlsx", sheetName = "Sheet1")
    #filter bed tissues
    cageTable <- as.data.table(cageTable)
    cageTable[is.na(Characteristics.Tissue.)] <- "unclassifiable"
    cageTable[Characteristics.Tissue. == "Urethra"]$Characteristics.Tissue. <- "urethra" 
    cageTable[Characteristics.Tissue. == "adipose tissue"]$Characteristics.Tissue. <-"adipose"
    
    insertTable(Matrix <- cageTable,tableName = "cageInformation")
  } else{
    cageTable <- readTable("cageInformation")
  }
  return(cageTable)
}
#' get all cage files that we have info on
getAllUsableCage <- function(cageFiles){
  matchCageIDandCageName <- rep("a",length(cageFiles))
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
  if(nrow(cageWeHave) < nrow(matchCageIDandCageName)) warning("did not find all cage experiments in info table")
  return(cageWeHave)
}

# get variance between different leader versions
getAllLeaderChanges <- function(){
  if(!file.exists(p(dataFolder,"/leaderOriginalWidths.rdata"))){
    getLeaders()
    widths <- ORFik:::widthPerGroup(fiveUTRs)
    save(widths, file = p(dataFolder,"/leaderOriginalWidths.rdata"))
    rm(fiveUTRs)
  }
  
  library(doParallel)
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  maxCores = as.integer(detectCores()/2)
  cl <- makeCluster(maxCores)
  registerDoParallel(cl)
  leadersList = list.files(leadersFolder)
  nLeadersList = length(leadersList)
  rm(fiveUTRs)
  output <- foreach(i=1:nLeadersList, .combine = 'rbind') %dopar% {
    source("./uorfomeGeneratorHelperFunctions.R")
    leadersList = list.files(leadersFolder)
    
    load(p(dataFolder,"/leaderOriginalWidths.rdata"))
    load(p(leadersFolder,leadersList[i]))
    widthsCage <- ORFik:::widthPerGroup(fiveUTRs)
    
    diffWidths <- widths - widthsCage
    same <- sum(diffWidths == 0)
    bigger <- sum(diffWidths < 0)
    smaller <- sum(diffWidths > 0)
    meanDifBigger <- mean(diffWidths[diffWidths < 0])
    meanDifSmaller <- mean(diffWidths[diffWidths > 0])
    return(c(same,bigger,smaller,meanDifBigger,meanDifSmaller))
  }
  dt <- as.data.table(matrix(output, ncol = 5))
  colnames(dt) <- c("same", "bigger", "smaller", "meanBigger", "meanSmaller")
  stopCluster(cl)
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  save(dt,file = "leaderWidthChanges.rdata")
}


#' Assign transcriptnames to orfs, and find rank for each orf
#' 
#' Given orfs and transcripts, find all transcripts the orfs are within
#' and name them by this. Also the second orf in 
linkORFsToTx <- function(){
  leaders <- leaderAllSpanning()
  grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
  overlaps <- findOverlaps(grl,leaders, type = "within")
  sortedIndeces <- order(to(overlaps))
  from <- from(overlaps)[sortedIndeces]
  to <- to(overlaps)[sortedIndeces]
  txNames <- names(leaders)[to]
  uorfIDs <- ORFik:::orfID(grl)[from]
  dt <- data.table(uorfID = uorfIDs, txNames = txNames)
  insertTable(Matrix = dt, tableName = "linkORFsToTx",rmOld = T)
  
  # now make grl with transcript mapping
  grl <- grl[from]
  names(grl) <- txNames
  grlb <- sortPerGroup(grl)
  grls <- ORFik:::makeORFNames(grlb)
  ranks <- ORFik:::makeExonRanks(grlb, T)
  asGR <- unlist(grlb, use.names = F)
  asGR$names <- paste0(names(asGR), "_", ranks)
  grlf <- groupGRangesBy(asGR, asGR$names)
  insertTable(Matrix = grlf, tableName = "uorfsAsGRWithTx",rmOld = T)
  
  
  # grlf <- sortPerGroup(grlf)
  # insertTable(Matrix = dt, tableName = "uorfsAsGRByTX",rmOld = T)
}

createCatalogueDB <- function(){
  uorfDB <- createDataBase(databaseName)
  createUniqueIDs()
  createUORFAtlas()
  rfpTables()
}
