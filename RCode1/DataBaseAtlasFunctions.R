
createUniqueIDsAtlasNonTranscript <- function(){
  stop("not working, not needed anymore!")
  A = readTable("MergedByExonsuniqueUORFs")
  c = paste(A$seqnames,A$start,A$end,A$strand)
  uniqueIDs = as.data.table(c)
  colnames(uniqueIDs) = "uorfID"
  insertTable(Matrix = uniqueIDs, tableName = "uorfIDNonTranscript")
}

allFeaturesAtlas <- function(){
  # first sequence features
  getSequenceFeatures()
  
  # then computeFeatures
  getAllFeaturesFromUorfs()
  # then rfpTables
  riboAtlasFPKMTissue() # fix this is wrong
  
  # then rna tables
  getRNAFpkms()
  # then do te tables
  getTeFeatures()
  
  # thats it right ?
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
}


createCatalogueDB <- function(){
  uorfDB <- createDataBase(databaseName)
  createUniqueIDs()
  createUORFAtlas()
  allFeaturesAtlas()
  
}
