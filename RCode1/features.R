# get only sequence features from orfik
getSequenceFeatures <- function(){
  grl <- getUorfsInDb(T, T, T)
  getAll()
  
  # kozak
  kozak <- kozakSequenceScore(grl, fa)
  orfID <- ORFik:::orfID(grl)
  dt <- data.table(uorfID = orfID, kozak = kozak)
  insertTable(dt, "kozak")
  # distORFCDS
  distORFCDS <- distToCds(grl, cageFiveUTRs, cds)
  insertTable(data.table(uorfID = orfID, distORFCDS = distORFCDS), "distORFCDS")
  # distORFTSS
  distORFTSS <- distToTSS(grl, cageFiveUTRs)
  insertTable(data.table(uorfID = orfID, distORFTSS = distORFTSS), "distORFTSS")
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
  starts <- startCodons(grl, is.sorted = T)
  getSequencesFromFasta(starts, isSorted = T)

  insertTable(data.table(uorfID = orfID, startCodon = as.character(seqs, use.names = F)), "StartCodons")
  # stop codon
  stops <- stopCodons(grl, is.sorted = T)
  getSequencesFromFasta(stops, isSorted = T)
  insertTable( data.table(uorfID = orfID, stopCodon = as.character(seqs, use.names = F)), "StopCodons")
  
  # exon-exon junctions
  eej <- numExonsPerGroup(fiveUTRs, T)
  link <- readTable("linkORFsToTx")
  eej <- as.integer(eej[link$txNames])
  insertTable(data.table(txNames = link$txNames, eej = eej), "exon-exonJunctionsLeader")
  
}

#' Get all features from grl
#' @param grl the ORFs to find features on, if null imports from database
getAllFeatures <- function(grl = NULL, RFPPath, RNAPath = NULL, i){  
  saveName <- paste0("featureTablesTemp/dt_",i,".rdata")
  if (file.exists(saveName))
    return(i)

  getAll()  

  if (is.null(RNAPath)) {
    RNA <- NULL
  } else {
    RNA <- readGAlignments(RNAPath)
  }
  
  if (is.null(grl)) {
    grl = getUorfsInDb(T, T, T)
  } 
  RFP <- fread.bed(filePath = RFPPath)
 
  dt <- ORFik:::computeFeaturesCage(grl = grl, RFP, RNA = RNA,
                                    fiveUTRs = GRangesList(), cds = cds, tx = tx,
                                    threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                                    orfFeatures = T, includeNonVarying = F, grl.is.sorted = T)
  save(dt,file = saveName)
  return(i)
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
