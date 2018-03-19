
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

clusterUorfFeature <- function(uorfFeature, indices = NULL, saveLocation){
  # cluster
  if (is.null(indices)) {
    dists <- dist(t(uorfFeature)) # transpose for tissue
  } else {
    dists <- dist(t(uorfFeature[indices])) # transpose for tissue
  }
  library(fastcluster)
  clustering <- hclust(dists)
  pdf(saveLocation)
  plot(clustering)
  dev.off()
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
  getAll()
  # kozak
  kozak <- kozakSequenceScore(grl, fa)
  orfID <- ORFik:::orfID(grl)
  dt <- data.table(uorfID = orfID, kozak = kozak)
  insertTable(dt, "kozak")
  # distORFCDS
  distORFCDS <- ORFik:::distOrfToCds(grl, fiveUTRs, cds, 1000)
  dt <- data.table(uorfID = orfID, distORFCDS = distORFCDS)
  insertTable(dt, "distORFCDS")
  # fractionLengths
  tx_len <- ORFik:::widthPerGroup(tx)
  fractionLengths <- fractionLength(grl, tx_len)
  dt <- data.table(uorfID = orfID, fractionLengths = fractionLengths)
  insertTable(dt, "fractionLengths")
  # inFrameCDS
  inFrameCDS <- ORFik:::inFrameWithCDS(distORFCDS)
  dt <- data.table(uorfID = orfID, inFrameCDS = inFrameCDS)
  insertTable(dt, "inFrameCDS")
  # isOverlappingCds
  isOverlappingCds <- isOverlappingCds(distORFCDS)
  dt <- data.table(uorfID = orfID, isOverlappingCds = isOverlappingCds)
  insertTable(dt, "isOverlappingCds")
  # rankInTx
  rankInTx <- OrfRankOrder(grl)
  dt <- data.table(uorfID = orfID, rankInTx = rankInTx)
  insertTable(dt, "rankInTx")
  
}

# for experiment on brain/ hek293
getAllFeatures <- function(grl, RFPPath, RNAPath = NULL, i){  

  getFasta()
  getCDS()
  getThreeUTRs()
  tx <- getTx()
  #tx <- ORFik:::extendLeaders(tx)
  
  #or with extension
  cageFiveUTRs <- leaderAllSpanning()
  
  # only bed here allowed!
  RFPShifted <- ORFik:::cageFromFile(RFPPath)
  if (is.null(RNAPath)) {
    RNA <- NULL
  } else {
    RNA <- readGAlignments(RNAPath)
  }
  
  dt <- ORFik:::computeFeatures(grl = grl, RFP = RFPShifted, RNA = NULL,
                            fiveUTRs = fiveUTRs, cds = cds, tx = tx,
                            threeUTRs = threeUTRs, faFile = fa, riboStart = 26, riboStop = 34,
                            extension = 1000, orfFeatures = T,
                            cageFiveUTRs = cageFiveUTRs, includeNonVarying = F)
  save(dt,file = paste0("featureTablesTemp/dt_",i,".rdata"))
  return(i)
}

