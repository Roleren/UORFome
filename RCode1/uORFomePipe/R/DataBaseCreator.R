
#' Create 1 column of all unique ids from uorfID folder
createUniqueIDs <- function() {
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
    insertTable(Matrix = allUniqueIDs, tableName = "uniqueIDs")
  } else {
    message("uniqueIDs already exist, skipping remake of them")
  }
}

#' convert to gr from string and filter NB!!! put this  in pipeline!!
createGRObjects <- function(makeBed = T){
  if (file.exists(p(dataBaseFolder,"/uniqueUorfsAsGRWithTx.rdata"))) {
    grl <- toGRFromUniqueID(readTable("uniqueIDs")$Matrix)
    uniqueIDs <- ORFik:::orfID(grl)
    save(grl, file = "./uniqueUorfsAsGR.rdata")
    insertTable(Matrix = uniqueIDs, tableName =  "uniqueIDs", rmOld = T)
    insertTable(Matrix = grl,tableName = "SplittedByExonsuniqueUORFs", rmOld = T)

    if (makeBed)
      bed12(grl, "bedUniqueUorfs.bed", T)

    # make all spanning cage leader from cage
    allLeadersSpanningLeader()
    # find tx matching
    linkORFsToTx()
  } else {
    message("GRObjects already exist, skipping remake of them")
  }
  return(NULL)
}

#' Create a data.table of true Fase
createUORFAtlas <- function(){
  if (!file.exists(paste0(getwd(),"/UORFAtlas.rdata"))) {
    uorfIDsAllUnique <- readTable("uniqueIDs")
    colnames(uorfIDsAllUnique) = "uorfID"
    uorfAtlas <- as.data.table(matrix(F, nrow = nrow(uorfIDsAllUnique), ncol = length(idFiles)+1))
    uorfAtlas[,"V1" := uorfIDsAllUnique$uorfID]
    colnames(uorfAtlas) = c("uorfID", as.character(1:(length(idFiles))))
    j = 1
    print(paste("creating uorfAtlas at index of max:", length(idFiles)))
    for(i in idFiles){
      load(p(idFolder, i))

      uorfAtlas[, as.character(j) :=  data.table::`%chin%`(uorfIDsAllUnique$uorfID, uorfID)]

      j = j+1
      print(j)
    }

    save(uorfAtlas,file = "UORFAtlas.rdata")
  } else {
    message("UORFAtlas already exist, skipping remake of them")
  }
}

#' Create tissueTable for cage, 1 row per unique uorf
#' Tissue must have at least 2 CAGE libraries supporting the uORF
#' to declare a hit in that tissue. If only one sample, that sample
#' must include uORF to be a hit.
getTissueTable <- function(){
  if (tableNotExists("tissueAtlasByCage")) {
    cageTable <- getCageInfoTable()
    uniqueTissues <- unique(cageTable$tissue)

    # load needed tables, and make tissue atlas of cage
    load("UORFAtlas.rdata")
    uorfIDs <- readTable("uniqueIDs")
    colnames(uorfIDs) = "uorfID"

    if(any(rowSums(uorfAtlas[,2:ncol(uorfAtlas)]) == 0))
      stop("uorfAtlas and unique uorf IDs does not match!")

    finalMatrix <- as.data.table(matrix(nrow =
      nrow(uorfIDs), ncol = length(uniqueTissues)+1))

    finalMatrix[, 1 := uorfIDs$uorfID]
    colnames(finalMatrix)[1] <- "uorfID"
    uniqueTissues <- as.character(uniqueTissues)
    colnames(finalMatrix)[2:ncol(finalMatrix)] <- uniqueTissues
    onlyOneCAGE <- c()
    for(i in 2:(length(uniqueTissues)+1)){
      cageFilestoCheck <- cageTable[tissue == uniqueTissues[i-1]]$cage_index
      if(length(cageFilestoCheck) == 1) onlyOneCAGE <- c(onlyOneCAGE, i)
    }

    for(i in 2:(length(uniqueTissues)+1)){
      cageFilestoCheck <- cageTable[tissue == uniqueTissues[i-1]]$cage_index
      cageFilestoCheck <- cageFilestoCheck + 1
      makeGroupingForColumn <- rowSums(uorfAtlas[,cageFilestoCheck, with = F])
      if (i %in% onlyOneCAGE) {
        finalMatrix[, uniqueTissues[i-1] := makeGroupingForColumn > 0]
      } else {
        finalMatrix[, uniqueTissues[i-1] := makeGroupingForColumn > 1]
      }
    }
    tissueAtlas <- finalMatrix

    save(tissueAtlas,file = "tissueAtlas.rdata")
    insertTable(Matrix = tissueAtlas,tableName = "tissueAtlasByCage", rmOld = T)
    return("ok tissueAtlassCage")
  }
  return("tissueAtlasCage already exists, stop if you want new")
}
