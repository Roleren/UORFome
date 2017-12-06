
createUniqueIDsAtlasNonTranscript <- function(){
  stop("not working, not needed anymore!")
  A = readTable("MergedByExonsuniqueUORFs")
  c = paste(A$seqnames,A$start,A$end,A$strand)
  uniqueIDs = as.data.table(c)
  colnames(uniqueIDs) = "uorfID"
  insertTable(Matrix = uniqueIDs, tableName = "uorfIDNonTranscript")
}


uniqueIdsAsGR <- function(){
  if(tableNotExists("SplittedByExonsuniqueUORFs")){
    uniqueIDs <- readTable("uniqueIDs")
    if (sum(duplicated(uniqueIDs)) > 0 ) stop("duplicated uorf names in uniqueIDs")
    grl <- toGRFromUniqueID(uniqueIDs)
    insertTable(Matrix = grl,tableName = "SplittedByExonsuniqueUORFs")
    #grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
    
    #now make uscs bed 12 format
    bed12(grl, "bedUniqueUorfs.bed", T)
    
  } else {
    grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
  }
  return(grl)
}


#' This is a check to see that pipeline have done everything correctly
#' if redoing the findOverlaps does not find all orfs within fiveUTRs
#' it means that some orfs are outside the mapping area
#' this should not happen!
validateExperiments <- function(grl){
  
  getGTF()
  getLeaders()
  getCDS()
  
  fiveUTRs <- addFirstCdsOnLeaderEnds(
    makeGrlAndFilter(extendsTSSexons(fiveUTRs), fiveUTRs), cds)
  a <- findOverlaps(query = unlist(grl, use.names = F), fiveUTRs)
  a <- a[!duplicated(from(a))]
  if(length(a) != length(unlist(grl))){ 
    stop("Not all orfs where within the FiveUTRs used
         to make them, something is wrong!")
  } else { print("experiments look valid")}
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
  insertTable(riboTable,"riboAll")
  if ((ncol(riboTable) - 1) != length(rpfFilePaths)){
    stop("something wrong in creation of riboTable")
  }
  return(riboTable)
}
#' Riboseq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is mean fpkm per tissue
riboAtlasFPKMTissue <- function(grl,rpfFilesPaths,riboTables,SpeciesGroup){
  # now do per tissue true/false
  rpfSamples <- SpeciesGroup[SpeciesGroup$Sample_Type == "RPF",]
  #1. we have tissues in speciesGroup
  tissuesUsed <- rpfSamples[rpfSamples$RnaRfpFolders %in% rpfFilePaths,]
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
  insertTable(Matrix = riboByTissue,tableName = "RiboByTissueTF")
  
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
  insertTable(Matrix = riboByTissueMean,tableName = "RiboByTissueMean")
}
