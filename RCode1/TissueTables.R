


#' Ribo-seq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is row-mean fpkm per tissue
#' @param riboDbName "Ribofpkm"
#' @param dbOutputNames the 2 output names c("RiboByTissueTF", "RiboByTissueMean")
riboAtlasFPKMTissue <- function(riboDbName = "Ribofpkm",
                                dbOutputNames = c("RiboByTissueTF", "RiboByTissueMean")){
  if(length(dbOutputNames) != 2) stop("dbOutputNames must have 2 character elements")
  
  # now do per tissue true/false
  rpfSamples <- getRiboRNAInfoTable()
  #1. we have tissues in link
  
  uniqueTissues <- as.character(unique(rpfSamples$tissue))
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
    print(i)
    riboColumns <- riboTable[,rpfSamples$tissue == i, with = F]
    riboByTissue[,i] <- rowSums(riboColumns > 1) > 1
  }
  riboByTissue <- data.table(idColumns, riboByTissue)
  insertTable(Matrix = riboByTissue, tableName = dbOutputNames[1],  rmOld = T)
  
  #now get mean value instead of true/false
  riboByTissueMean <- riboByTissueTemp
  rm(riboByTissue)
  rm(riboByTissueTemp)
  # rowMeans(riboColumns) mean per orf in tissue
  for(i in uniqueTissues){
    print(i)
    riboColumns <- riboTable[,rpfSamples$tissue == i, with = F]
    riboByTissueMean[,i] <- rowMeans(riboColumns)
  }
  riboByTissueMean <- data.table(idColumns, riboByTissueMean)
  insertTable(Matrix = riboByTissueMean, tableName = dbOutputNames[2], rmOld = T)
}

#' RNA-seq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is row-mean fpkm per tissue
#' @param nIDColumns 1L transcript names
rnaAtlasFPKMTissue <- function(){
  # now do per tissue true/false
  rnaSamples <- getRiboRNAInfoTable()
  #1. we have tissues in link
  
  uniqueTissues <- as.character(unique(rnaSamples$tissue))
  #2. we have all the tables in riboTables
  #3. So for each tissue, find the group, then for each group ->
  rnaTable <- readTable("RNAfpkm")
  idColumns <- getIDColumns(rnaTable)
  rnaTable <- removeIDColumns(rnaTable)
  
  rnaByTissue <- as.data.table(matrix(nrow = nrow(rnaTable),
                                      ncol = length(uniqueTissues)))
  colnames(rnaByTissue) <- uniqueTissues
  rnaByTissueTemp <- rnaByTissue
  
  #4. rowSum(riboColumns > 1) > 1
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf, it will be true
  for(i in uniqueTissues){
    print(i)
    rnaColumns <- rnaTable[,rnaSamples$tissue == i, with = F]
    rnaByTissue[,i] <- rowSums(rnaColumns > 1) > 1
  }
  rnaByTissue <- data.table(idColumns, rnaByTissue)
  insertTable(Matrix = rnaByTissue, tableName = "RNAByTissueTF", rmOld = T)
  
  #now get mean value instead of true/false
  rnaByTissueMean <- rnaByTissueTemp
  rm(rnaByTissue)
  rm(rnaByTissueTemp)
  # rowMeans(riboColumns) mean per orf in tissue
  for(i in uniqueTissues){
    print(i)
    rnaColumns <- rnaTable[,rnaSamples$tissue == i, with = F]
    rnaByTissueMean[,i] <- rowMeans(rnaColumns)
  }
  rnaByTissueMean <- data.table(idColumns, rnaByTissueMean)
  insertTable(Matrix = rnaByTissueMean, tableName = "RNAByTissueMean", rmOld = T)
  return(NULL)
}

#' In tf and
teAtlasTissue <- function(TFDbName = c("RiboByTissueTF", "RNAByTissueTF"),
                          MeanDbName = c("RiboByTissueMean", "RNAByTissueMean"),
                          dbOutputNames = c("TEByTissueTF", "TEByTissueMeanWithInf", "TEByTissueMeanWithoutInf")){
  
  riboTF <- readTable(TFDbName[1])
  rnaTF <- readTable(TFDbName[2]) # check that amount match
  rnaTF <- matchByTranscript(rnaTF, riboTF)
  if(nrow(rnaTF) != nrow(riboTF)) stop("not equal nrow in ribo and rna TF")

  teTF <- copy(riboTF)
  
  uniqueTissues <- colnames(teTF)
  uniqueTissues <- uniqueTissues[-1] # 1 col only, FIX IF NEEDED!!!
  
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf,
  # #for both ribo and rna-seq, it will be true
  for(i in uniqueTissues){
    print(i)
    teTF[, i] <- rnaTF[, i, with = F] & riboTF[, i, with = F]
  }
  
  insertTable(Matrix = teTF, tableName = dbOutputNames[1], rmOld = T)
  
  teMean <- teTF
  rm(teTF)
  rm(rnaTF)
  rm(riboTF)
  
  riboMean <- readTable(MeanDbName[1])
  rnaMean <- readTable(MeanDbName[2])
  rnaMean <- matchByTranscript(rnaMean, riboMean)
  if(nrow(rnaMean) != nrow(riboMean)) stop("not equal nrow in ribo and rna Mean")
  
  # with non-finite values
  for(i in uniqueTissues){
    teMean[, i] <- rnaMean[, i, with = F] / riboMean[, i, with = F]
  }
  
  insertTable(Matrix = teMean, tableName = dbOutputNames[2])
  
  # now do without non-finite
  DTtemp <- teMean
  DT <- DTtemp[, (2):ncol(DTtemp)] # dangerous id remover here! change if needed!
  DT <- removeNonFinite(DT)
  DTtemp[, (2):ncol(DTtemp)] <- DT
  insertTable(Matrix = DTtemp, tableName = dbOutputNames[3])
  return(NULL)
}

teAtlasTissueNew <- function(inputDT, colExclusion = "fpkmRFP_", dbOutputNames = 
                               c("TEByTissueMean")){
  info <- getRiboRNAInfoTable()
  
  pattern <- colExclusion
  inputDTNon <- removeIDColumns(inputDT)
  indices <- as.integer(gsub(pattern = pattern, replacement = "", x = colnames(inputDTNon)))
  if(length(indices) == 0) stop("could not find te indices from colExclusion")
  tissues <- info$tissue[indices] 
  
  uniques <- unique(tissues)
  
  if(!is.null(inputDT$uorfIDs)){
    dt <- data.table(uorfIDs = inputDT$uorfIDs,
                     txNames = inputDT$txNames)
  } else {
    dt <- data.table(txNames = inputDT$txNames)
  }
  
  for(tissue in uniques) {
    which <- tissues == tissue
    dt <- data.table(dt, rowMeans(inputDTNon[,which, with = F]))
  }
  
  if(!is.null(inputDT$uorfIDs)){
    colnames(dt) <- c("uorfIDs", "txNames", uniques )
  } else {
    colnames(dt) <- c("txNames", uniques )
  }
  
  return(dt)
}