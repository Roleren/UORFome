

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

#' Make a leader that spans the max leader per transcript by cage files
allLeadersSpanningLeader <- function(){
  
  leadersList = list.files(leadersFolder)
  nLeadersList = length(leadersList)
  
  # Get all width of cage experiment TSS reassignments
  widths <- foreach(i=1:nLeadersList, .combine = 'cbind') %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./uorfomeGeneratorHelperFunctions.R")
    leadersList = list.files(leadersFolder)
    
    usedCage = gsub(pattern = ".leader.rdata",replacement = "",x = leadersList[i]) 
    
    
    load(p(leadersFolder,leadersList[i]))
    
    
    return(ORFik:::widthPerGroup(fiveUTRs, keep.names = F))
  }
  maxWidths <- rowMaxs(widths)
  getCDS()
  if( exists("fiveUTRs", mode = "S4")) {
    rm(fiveUTRs)
  }
  getLeaders()
  
  fiveWithCds <- ORFik:::addFirstCdsOnLeaderEnds(fiveUTRs, cds)
  originalWidths <- ORFik:::widthPerGroup(fiveWithCds, F)
  
  fEWithCds <- ORFik:::firstExonPerGroup(fiveWithCds) 
  fEWithCDSWidths <- ORFik:::widthPerGroup(fEWithCds, F)
  excludedFirstWidth <- originalWidths - fEWithCDSWidths
  
  unl <- unlist(fiveUTRs, use.names = T)
  firstExons <- unl[unl$exon_rank == 1]
  
  changes <- ((maxWidths - excludedFirstWidth) - (originalWidths - excludedFirstWidth))
  posStrand <- strandBool(firstExons)
  start(firstExons[posStrand]) <- start(firstExons[posStrand]) - changes[posStrand]
  end(firstExons[!posStrand]) <- end(firstExons[!posStrand]) + changes[!posStrand]
  
  backListFive <- unl
  backListFive[backListFive$exon_rank == 1] <- firstExons
  
  f <- groupGRangesBy(backListFive)
  fwc <- ORFik:::addFirstCdsOnLeaderEnds(f, cds)
  
  fwcWidths <- ORFik:::widthPerGroup(fwc, F)
  
  if(!all(fwcWidths == maxWidths)) {
    stop("Algorithm is wrong for five extension!")
  }
  
  # all ok, then save
  source("./DataBaseSetup.R")
  CageFiveUTRs <- f
  save(CageFiveUTRs, file = "CageFiveUTRs.rdata")
  CageFiveWithCDS <- fwc
  save(CageFiveWithCDS, file = "CageFiveUTRsWithCDS.rdata")
  
}

#' Assign transcriptnames to orfs, and find rank for each orf
#' 
#' Given orfs and transcripts, find all transcripts the orfs are within
#' and name them by this. Also the second orf in 
linkORFsToTx <- function(){
  leaders <- leaderCage()
  grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
  overlaps <- findOverlaps(grl,leaders, type = "within")
  if( length(unique(from(overlaps))) != length(grl)) {
    stop("leader is not spanning all uORFs, check for uORFs going into cds.")
  }
  sortedIndeces <- order(to(overlaps))
  from <- from(overlaps)[sortedIndeces]
  to <- to(overlaps)[sortedIndeces]
  txNames <- names(leaders)[to]
  uorfIDs <- ORFik:::orfID(grl)[from]
  dt <- data.table(uorfID = uorfIDs, txNames = txNames)
  insertTable(Matrix = dt, tableName = "linkORFsToTx",rmOld = T)
  
  # now make grl with transcript mapping
  grl <- grl[from]
  grlb <- grl
  names(grlb@unlistData) <- NULL
  names(grlb) <- txNames
  asGR <- unlist(grlb, use.names = T)
  names(grlb@unlistData) <- names(asGR)
  
  grls <- ORFik:::makeORFNames(grlb)
  ranks <- ORFik:::makeExonRanks(grlb, T)
  asGR <- unlist(grlb, use.names = F)
  asGR$names <- paste0(names(asGR), "_", ranks)
  grlf <- groupGRangesBy(asGR, asGR$names)
  insertTable(Matrix = grlf, tableName = "uorfsAsGRWithTx",rmOld = T)
}


createCatalogueDB <- function(){
  #uorfDB <- createDataBase(databaseName)
  createUniqueIDs()
  createGRObjects()
  createUORFAtlas()
  allFeaturesAtlas()
  
}
