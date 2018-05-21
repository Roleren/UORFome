
#' Create 1 column of all unique ids from uorfID folder
createUniqueIDs <- function(){
  if (length(idFiles) == 0) {
    stop("idFiles can not have 0 length")
  }
  if (tableNotExists("uniqueIDs")) {
    j = 1
    for(i in idFiles){
      load(p(idFolder, i))
      uorfID <- unique(uorfID)
      if (j == 1) {
        allUniqueIDs <- uorfID
      }else{
        matching <- data.table::`%chin%`(uorfID, allUniqueIDs)
        toAdd <- uorfID[which(matching == F)]
        allUniqueIDs <- c(allUniqueIDs,toAdd)
      }
      j <- j+1
    } 
    allUniqueIDs <- allUniqueIDs[data.table::chorder(allUniqueIDs)]
    
    #save(allUniqueIDs,file = "allUniqueIDs.rdata")
    insertTable(Matrix = allUniqueIDs,tableName = "uniqueIDs")
  } else {
    message("uniqueIDs already exist, skipping remake of them")
  }
}

#' #' Create 1 column of all unique ids from uorfID folder
#' createUniqueIDsFast <- function(){
#'   if (length(idFiles) == 0) {
#'     stop("idFiles can not have 0 length")
#'   }
#'   if (tableNotExists("uniqueIDs")) {
#'     
#'     j = 1
#'     
#'     for(i in idFiles[1:50]){
#'       load(p(idFolder, i))
#'       uorfID <- unique(uorfID)
#'       if (j == 1) {
#'         allUniqueIDs <- uorfID
#'       }else{
#'         
#'         matching <- data.table::`%chin%`(uorfID,allUniqueIDs)
#'         toAdd <- uorfID[which(matching == F)]
#'         allUniqueIDs <- c(allUniqueIDs,toAdd)
#'       }
#'       j <- j+1
#'     }
#'     #data.table::`%chin%`(a,allUniqueIDs[1:100])
#'     c <- foreach(i=seq_along(idFiles), .combine=function(x,y){unique(c(unique(x), unique(y)))}) %dopar% {
#'       setwd("/export/valenfs/projects/uORFome/RCode1/")
#'       source("./uorfomeGeneratorHelperFunctions.R")
#'       load(p(idFolder, idFiles[i]))
#'       return(uorfID)
#'     }
#'     
#'     allUniqueIDs <- sort(allUniqueIDs)
#'     #save(allUniqueIDs,file = "allUniqueIDs.rdata")
#'     insertTable(Matrix = allUniqueIDs,tableName = "uniqueIDs")
#'   } else {
#'     message("uniqueIDs already exist, skipping remake of them")
#'   }
#' }

#' convert to gr from string and filter NB!!! put this  in pipeline!!
createGRObjects <- function(makeBed = T){
  if (!file.exists(paste0(getwd(),"/uniqueUorfsAsGR.rdata"))) {
    uniqueIDs <- readTable("uniqueIDs")
    grl <- toGRFromUniqueID(uniqueIDs$Matrix)
    getCDS()
    
    # filter out uORFs with same start as cds
    starts <- startSites(grl, asGR = T, is.sorted = T)
    cdsstarts <- startSites(cds, asGR = T, is.sorted = T)
    overlaps <- findOverlaps(starts, cdsstarts, type = "within")
    
    grl <- grl[-unique(from(overlaps))]
    uniq <- uniqueIDs$Matrix
    uniq <- uniq[-unique(from(overlaps))]
    save(grl, file = "./uniqueUorfsAsGR.rdata")
    insertTable(Matrix = uniq, tableName =  "uniqueIDs", rmOld = T)
    insertTable(Matrix = grl,tableName = "SplittedByExonsuniqueUORFs", rmOld = T)
    
    if (makeBed)
      bed12(grl, "bedUniqueUorfs.bed", T)
    
    # make all spanning cage leader from cage
    allLeadersSpanningLeader()
    # find tx matching
    linkORFsToTx()
  }
}

#' Create a data.table of true Fase
createUORFAtlas <- function(){
  if (!file.exists(paste0(getwd(),"/UORFAtlas.rdata"))) {
    uorfIDsAllUnique <- readTable("uniqueIDs")
    colnames(uorfIDsAllUnique) = "uorfID"
    uorfAtlas <- as.data.table(matrix(F, nrow = nrow(uorfIDsAllUnique), ncol = length(idFiles)+1))
    uorfAtlas[,1] <- uorfIDsAllUnique
    colnames(uorfAtlas) = c("uorfID", as.character(1:(length(idFiles))))
    j = 1
    for(i in idFiles){
      load(p(idFolder, i))
      
      uorfs <- uorfID[!duplicated(uorfID)]
      
      uorfAtlas[,(j+1)] <- data.table::`%chin%`(uorfIDsAllUnique$uorfID, uorfs)
        
      j = j+1
    }
    
    save(uorfAtlas,file = "UORFAtlas.rdata")
  } else {
    message("UORFAtlas already exist, skipping remake of them")
  }
}

#' Create tissueTable for cage, 1 row per unique uorf 
getTissueTable <- function(){
  if (!tableNotExists("tissueAtlasByCage")) {
    cageTable <- getCageInfoTable()
    uniqueTissues <- unique(cageTable$Characteristics.Tissue.)  
    
    cageWeHave <- getAllUsableCage(cageTable)
    
    # load needed tables, and make tissue atlas of cage
    load("UORFAtlas.rdata")
    uorfIDs <- readTable("uniqueIDs")
    colnames(uorfIDs) = "uorfID"
    uorfAtlasRows0 <- which(rowSums(uorfAtlas[,2:ncol(uorfAtlas)]) == 0)
    if(length(uorfAtlasRows0) > 0) 
      stop("uorfAtlas and unique uorf IDs does not match!")
    
    finalMatrix <- as.data.table(matrix(nrow =
      nrow(uorfIDs), ncol = length(uniqueTissues)+1))
  
    finalMatrix[, 1 := uorfIDs$uorfID]
    colnames(finalMatrix)[1] <- "uorfID"
    uniqueTissues <- as.character(uniqueTissues)
    colnames(finalMatrix)[2:ncol(finalMatrix)] <- uniqueTissues
    
    for(i in 2:(length(uniqueTissues)+1)){
      cageFilestoCheck <- cageWeHave[Characteristics.Tissue. == uniqueTissues[i-1]]$cage_index
      cageFilestoCheck <- cageFilestoCheck + 1
      makeGroupingForColumn <- rowSums(uorfAtlas[,cageFilestoCheck, with = F])
  
      finalMatrix[, uniqueTissues[i-1] := makeGroupingForColumn > 1]
    }
    tissueAtlas <- finalMatrix
    
    save(tissueAtlas,file = "tissueAtlas.rdata")
    insertTable(Matrix = tissueAtlas,tableName = "tissueAtlasByCage")
    return("ok tissueAtlassCage")
  }
  return("tissueAtlasCage already exists, stop if you want new")
}



# fix the presentation
# Take all uorfs, and cluster the tissue based on uorfs
# can you recover the cell type based on uorf usage
# hclust
