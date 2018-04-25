
getTissue <- function(){
  name = NA
  if (exists("tissue")){
    name = tissue
  } else if (exists("cageName")){
    
  }
  cbind(rep(1:length(transcriptName),name))
}

getPassFilter <- function(FPKMRFP,RPKMRNA, passFilter = 0.1){
  pass_filter <- RPKMRNA > passFilter & FPKMRFP > passFilter
  pass_filter[is.na(pass_filter)] <- F
  pass_filter
}

getORFnames <- function(unfilteredNames){
  gsub(".*\\.","", unfilteredNames)
}

clusterUorfFeature <- function(uorfFeature, indices = NULL, saveLocation, folder = "/export/valenfs/projects/uORFome/dataBase/clustering/"){
  # cluster
  if (is.null(indices)) {
    dists <- dist(t(uorfFeature)) # transpose for tissue
  } else {
    dists <- dist(t(uorfFeature[indices])) # transpose for tissue
  }
  pdfPath <- paste0(folder, saveLocation)
  library(fastcluster)
  clustering <- hclust(dists)
  if (ncol(uorfFeature) > 50) {
    pdf(pdfPath, width = 20, height = 12)
  } else {
    pdf(pdfPath)
  }
  
  plot(clustering)
  dev.off()
}

clusterAllUorfFeatures <- function(){
  # first filter out non-feature tables
  list <- listTables()
  list <- list[-grep(x = list, pattern = "cds", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "uorf", ignore.case = T)]
  list <- list[-grep(x = list, pattern = "information", ignore.case = T)]
  list <- list[-grep(x = list, pattern = "Variance", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "link", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "unique", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "Info", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "Tissue", ignore.case = T)]
  list <- list[-grep(x = list, pattern = "teUnfiltered", ignore.case = T)]
  
  riboInfo <- readTable("RiboSeqInfo")
  rnaInfo <- readTable("RNASeqInfo")
  for (feature in list) {
    tab <- readTable(feature, with.IDs = F)
    if(ncol(tab) == 103) { #rfp seq
      colnames(tab) <- riboInfo$Tissue.Cell_line
    } else if (ncol(tab) == 43) {
      colnames(tab) <- rnaInfo$Tissue.Cell_line
    } else if(ncol(tab) == 35) {
      indices <- as.integer(gsub(pattern = "fpkmRFP_", replacement = "", x = colnames(tab)))
      colnames(tab) <- riboInfo$Tissue.Cell_line[indices]
    } else {
      if (ncol(tab) != 1) {
        stop(paste("could not find ncol for feature: ", feature))
      } 
    }
    
    if (ncol(tab) != 1) {
      clusterUorfFeature(uorfFeature = tab, saveLocation = paste0(feature, ".pdf"))
    } 
    
  }
}

#' How good are the features
#' 
#' Te will be used as translation
#' Use cds as positive set and randomized columns of 
#' cds as negative.
#' 
#' Can a random forrest seperate these two sets ?
correlateFeatures <- function(){
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = F)
  cdsKozak <- readTable("cdsKozak", with.IDs = F)
  cdsORFScore <- readTable("cdsORFScore", with.IDs = F)
  cdsRfpFPKMs <- readTable("cdsRfpFPKMs", with.IDs = F)
  
  dt <- data.table(cdsTEs[,1], cdsKozak[,1], cdsORFScore[,1], cdsRfpFPKMs[,1])
  
  # add these to cds
  # ribosomeReleaseScore
  # insideOutsideORF
  # ribosomeStalingScore
  
  dtRandom <- dt
  for (i in 1:ncol(dt)) {
    ran <- sample(nrow(dtRandom))
    dtRandom[, i] <- dtRandom[ran, i, with = F]
  }
  
  merged <- rbindlist(list(dt,dtRandom))
  l <- nrow(dtRandom)
  y <- c(rep(1,l),rep(0,l))
  
  max <- nrow(merged)
  testMerged <- rbindlist(list(dt[1:max,],dtRandom[1:max,]))
  testY <- as.factor(c(rep(1,max),rep(0,max)))
  
  
  rf <- randomForest(x = merged, y = y)
  
}


uidFromCage <- function(cage = standardCage, asUID = TRUE,
                        with.transcriptNames = TRUE){
  
  rm(cageFiveUTRs)
  rm(fiveUTRs)
  getCDS()
  getThreeUTRs()
  getLeaders()
  cageFiveUTRs <- ORFik:::reassignTSSbyCage(fiveUTRs, standardCage, 1000, 1, cds)
  originalUorfsByTx <- getUnfilteredUORFs(cageFiveUTRs, assignRanges = F)
  gr <- unlist(originalUorfsByTx, use.names = F)
  grl <- groupGRangesBy(gr, gr$names)
  grl <- removeORFsWithinCDS(grl)
  
  if (!asUID) {
    return(grl)
  }
  uids <- toUniqueIDFromGR(grl)
  if (with.transcriptNames) {
    return(paste(uids, ORFik:::OrfToTxNames(grl)))
  }
  return(uids)
}

# get only sequence features from orfik
getSequenceFeatures <- function(){
  grl <- getUorfsInDb()
  getAll(extendTx = T)
  
  # kozak
  kozak <- kozakSequenceScore(grl, fa)
  orfID <- ORFik:::orfID(grl)
  dt <- data.table(uorfID = orfID, kozak = kozak)
  insertTable(dt, "kozak")
  # distORFCDS
  distORFCDS <- ORFik:::distToCds(grl, fiveUTRs, cds, 1000)
  insertTable(data.table(uorfID = orfID, distORFCDS = distORFCDS), "distORFCDS")
  # fractionLengths
  tx_len <- ORFik:::widthPerGroup(tx)
  fractionLengths <- fractionLength(grl, tx_len)
 
  insertTable(data.table(uorfID = orfID, fractionLengths = fractionLengths), "fractionLengths")
  # inFrameCDS
  inFrameCDS <- ORFik:::isInFrame(distORFCDS)

  insertTable(data.table(uorfID = orfID, inFrameCDS = inFrameCDS), "inFrameCDS")
  # isOverlappingCds
  isOverlappingCds <- isOverlapping(distORFCDS)
  insertTable(data.table(uorfID = orfID, isOverlappingCds = isOverlappingCds), "isOverlappingCds")
  # rankInTx
  rankInTx <- rankOrder(grl)

  insertTable(data.table(uorfID = orfID, rankInTx = rankInTx), "rankInTx")
  # number of uorfs per tx
  txNames <- txNames(grl)
  numberOfUorfsPerTx <- S4Vectors::Rle(txNames)

  insertTable(data.table(txNames = runValue(numberOfUorfsPerTx),
                         nUorfs = runLength(numberOfUorfsPerTx)), "numberOfUorfsPerTx")
  # start codon
  starts <- startCodons(grl)
  getSequencesFromFasta(starts, isSorted = T)

  insertTable(data.table(uorfID = orfID, startCodon = as.character(seqs)), "StartCodons")
  # stop codon
  stops <- stopCodons(grl)
  getSequencesFromFasta(stops, isSorted = T)
  insertTable( data.table(uorfID = orfID, stopCodon = as.character(seqs)), "StopCodons")
}

#' Get all features from grl
#' @param grl the ORFs to find features on, if null imports from database
getAllFeatures <- function(grl = NULL, RFPPath, RNAPath = NULL, i){  
  saveName <- paste0("featureTablesTemp/dt_",i,".rdata")
  if (file.exists(saveName))
    return(i)

  getAll()  
  #or with extension
  fiveUTRs <- leaderCage(F)
  tx <- extendLeaders(tx, fiveUTRs)

  if (is.null(RNAPath)) {
    RNA <- NULL
  } else {
    RNA <- readGAlignments(RNAPath)
  }
  
  if (is.null(grl)) {
    dt <- ORFik:::computeFeaturesCage(grl = getUorfsInDb(), RFP = ORFik:::fread.bed(RFPPath), RNA = RNA,
                                      fiveUTRs = fiveUTRs, cds = cds, tx = tx,
                                      threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                                      orfFeatures = T, includeNonVarying = F, extension = 0,
                                      grl.is.sorted = T)
  } else {
    dt <- ORFik:::computeFeaturesCage(grl = grl, RFP = ORFik:::fread.bed(RFPPath), RNA = RNA,
                                      fiveUTRs = fiveUTRs, cds = cds, tx = tx,
                                      threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                                      orfFeatures = T, includeNonVarying = F, extension = 0,
                                      grl.is.sorted = T)
  }
  save(dt,file = saveName)
  return(i)
 }
