

#' Ribo-seq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is row-mean fpkm per tissue
#' @param riboDbName "Ribofpkm"
#' @param dbOutputNames the 2 output names c("RiboByTissueTF", "RiboByTissueMean")
riboAtlasFPKMTissue <- function(riboDbName = "Ribofpkm",
                                dbOutputNames = c("RiboByTissueTF", "RiboByTissueMean")){
  if(length(dbOutputNames) != 2) stop("dbOutputNames must have 2 character elements")
  
  # now do per tissue true/false
  link <- readTable("linkRnaRfp")
  rpfSamples <- link[link$Sample_Type == "RPF",]
  #1. we have tissues in link
  
  uniqueTissues <- as.character(unique(rpfSamples$Tissue_or_CellLine))
  #2. we have all the tables in riboTables
  #3. So for each tissue, find the group, then for each group ->
  riboTable <- readTable(riboDbName)
  idColumns <- getIDColumns(riboTable)
  riboTable <- removeIDColumns(riboTable)
   # number of id columns used
  riboByTissue <- as.data.table(matrix(nrow = nrow(riboTable),
                                       ncol = length(uniqueTissues)))
  colnames(riboByTissue) <- uniqueTissues
  riboByTissueTemp <- riboByTissue
  #4. rowSum(riboColumns > 1) > 1
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf, it will be true
  for(i in uniqueTissues){
    indices <- rpfSamples[(rpfSamples$Tissue_or_CellLine == i),originalIndex]
    riboColumns <- riboTable[,indices, with=F]
    riboByTissue[,i] <- rowSums(riboColumns > 1) > 1
    
  }
  riboByTissue <- data.table(idColumns, riboByTissue)
  insertTable(Matrix = riboByTissue, tableName = dbOutputNames[1])
  
  #now get mean value instead of true/false
  riboByTissueMean <- riboByTissueTemp
  rm(riboByTissue)
  rm(riboByTissueTemp)
  # rowMeans(riboColumns) mean per orf in tissue
  for(i in uniqueTissues){
    indices <- rpfSamples[(rpfSamples$Tissue_or_CellLine == i),originalIndex]
    riboColumns <- riboTable[,indices, with=F]
    riboByTissueMean[,i] <- rowMeans(riboColumns)
  }
  riboByTissueMean <- data.table(idColumns, riboByTissueMean)
  insertTable(Matrix = riboByTissueMean, tableName = dbOutputNames[2])
}

#' RNA-seq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is row-mean fpkm per tissue
#' @param nIDColumns 1L transcript names
rnaAtlasFPKMTissue <- function(nIDColumns = 1L){
  # now do per tissue true/false
  link <- readTable("linkRnaRfp")
  rnaSamples <- link[link$Sample_Type == "RNA",]
  #1. we have tissues in link
  
  uniqueTissues <- as.character(unique(rnaSamples$Tissue_or_CellLine))
  #2. we have all the tables in riboTables
  #3. So for each tissue, find the group, then for each group ->
  rnaTable <- readTable("RNAfpkm")
  rnaByTissue <- as.data.table(matrix(nrow = nrow(rnaTable),
                                      ncol = (length(uniqueTissues)+nIDColumns)))
  rnaByTissue[,1:nIDColumns] <- rnaTable[,1:nIDColumns]
  colnames(rnaByTissue)[1:nIDColumns] <- colnames(rnaTable)[1:nIDColumns]
  colnames(rnaByTissue)[(nIDColumns+1):ncol(rnaByTissue)] <- uniqueTissues
  
  #4. rowSum(riboColumns > 1) > 1
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf, it will be true
  for(i in uniqueTissues){
    indices <- rnaSamples[(rnaSamples$Tissue_or_CellLine == i),originalIndex]
    indices <- indices + nIDColumns
    rnaColumns <- rnaTable[,indices, with=F]
    rnaByTissue[,i] <- rowSums(rnaColumns > 1) > 1
    
  }
  insertTable(Matrix = rnaByTissue, tableName = "RNAByTissueTF", rmOld = T)
  
  #now get mean value instead of true/false
  rnaByTissueMean <- rnaByTissue
  rm(rnaByTissue)
  # rowMeans(riboColumns) mean per orf in tissue
  for(i in uniqueTissues){
    indices <- rnaSamples[(rnaSamples$Tissue_or_CellLine == i),originalIndex]
    indices <- indices + nIDColumns
    rnaColumns <- rnaTable[,indices, with=F]
    rnaByTissueMean[,i] <- rowMeans(rnaColumns)
  }
  insertTable(Matrix = rnaByTissueMean, tableName = "RNAByTissueMean", rmOld = T)
}

#' In tf and
teAtlasTissue <- function(TFDbName = c("RiboByTissueTF", "RNAByTissueTF"),
                          MeanDbName = c("RiboByTissueMean", "RNAByTissueMean"),
                          dbOutputNames = c("TEByTissueTF", "TEByTissueMeanWithInf", "TEByTissueMeanWithoutInf")){
  
  riboTF <- readTable(TFDbName[1])
  rnaTF <- readTable(TFDbName[1]) # check that amount match
  if(nrow(rnaTF) != nrow(ribTF)) stop("not equal nrow in ribo and rna TF")
  
  teTF <- riboTF
  nIDColumns <- 2L # number of id columns used
  uniqueTissues <- colnames(teTF)[(nIDColumns+1):ncol(teTF)]
  
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf,
  # #for both ribo and rna-seq, it will be true
  for(i in uniqueTissues){
    
    teTF[, i] <- rnaTF[, i, with = F] & riboTF[, i, with = F]
  }
  insertTable(Matrix = teTF, tableName = dbOutputNames[1], rmOld = T)
  
  teMean <- teTF
  rm(teTF)
  rm(rnaTF)
  rm(riboTF)
  
  
  riboMean <- readTable(MeanDbName[1])
  rnaMean <- readTable(MeanDbName[2])
  
  # with non-finite values
  for(i in uniqueTissues){
    
    teMean[, i] <- rnaMean[, i, with = F]/riboMean[, i, with = F]
  }
  
  insertTable(Matrix = teMean, tableName = dbOutputNames[2])
  
  # now do without non-finite
  DTtemp <- teMean
  DT <- DTtemp[, (nIDColumns+1):ncol(DTtemp)]
  DT <- removeNonFinite(DT)
  DTtemp[, (nIDColumns+1):ncol(DTtemp)] <- DT
  insertTable(Matrix = DTtemp, tableName = dbOutputNames[3])
}