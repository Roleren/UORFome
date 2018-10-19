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
  eejuORF <- numExonsPerGroup(grl)
  insertTable(data.table(uorfID = orfID, eejuORF = eejuORF), "exon-exonJunctionsuORFs")
  
  # gene transcript connections
  genes <- GenomicFeatures::transcriptsBy(Gtf, by = "gene")
  unlGenes <- unlist(genes, use.names = TRUE)
  dt <- data.table(txNames = unlGenes$tx_name, geneNames = names(unlGenes))
  matches <- data.table::chmatch(ORFik:::txNames(grl), dt$txNames)
  dt <- dt[matches]
  insertTable(dt[matches], "uORFTxToGene")
  # go Terms
  uorfGo <- getORFsGoTerms(dt$geneNames)
  insertTable(data.table(go = uorfGo), "goTerms")
  
}

#' Get all features from grl
#' @param grl the ORFs to find features on, if null imports from database
getAllFeatures <- function(grl = NULL, RFPPath, RNAPath = NULL, i){  
  saveName <- paste0(dataBaseFolder, "/featureTablesTemp/dt_",i,".rdata")
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
  dt[,startCodonCoverage := countOverlaps(startCodons(grl,T), RFP)]
  save(dt,file = saveName)
  return(i)
}
