#' Get sequence features from orfik
#'
getSequenceFeatures <- function(){
  if(!tableNotExists("kozak") & !tableNotExists("gcContent")) {
    print("sequence features exist, stop if wrong!")
    return(NULL)
  }
  grl <- getUorfsInDb()
  getAll()

  # kozak
  kozak <- kozakSequenceScore(grl, tx, fa)
  dt <- data.table(kozak = kozak)
  insertTable(dt, "kozak")
  # distORFCDS
  distORFCDS <- distToCds(grl, cageFiveUTRs, cds)
  insertTable(data.table(distORFCDS = distORFCDS), "distORFCDS")
  # distORFTSS
  distORFTSS <- distToTSS(grl, cageFiveUTRs)
  insertTable(data.table(distORFTSS = distORFTSS), "distORFTSS")
  # fractionLengths
  tx_len <- ORFik:::widthPerGroup(tx)
  fractionLengths <- fractionLength(grl, tx_len)
  insertTable(data.table(fractionLengths = fractionLengths), "fractionLengths")
  # inFrameCDS
  inFrameCDS <- ORFik:::isInFrame(distORFCDS)
  insertTable(data.table(inFrameCDS = inFrameCDS), "inFrameCDS")
  # isOverlappingCds
  isOverlappingCds <- isOverlapping(distORFCDS)
  insertTable(data.table(isOverlappingCds = isOverlappingCds), "isOverlappingCds")
  # rankInTx
  rankInTx <- rankOrder(grl)

  insertTable(data.table(rankInTx = rankInTx), "rankInTx")
  # number of uorfs per tx
  txNames <- txNames(grl)
  numberOfUorfsPerTx <- S4Vectors::Rle(txNames)

  insertTable(data.table(nUorfs = runLength(numberOfUorfsPerTx)), "numberOfUorfsPerTx")
  # start codon
  starts <- startCodons(grl, is.sorted = T)
  getSequencesFromFasta(starts, isSorted = T)

  insertTable(data.table(startCodon = as.character(seqs, use.names = F)), "StartCodons")
  # stop codon
  stops <- stopCodons(grl, is.sorted = T)
  getSequencesFromFasta(stops, isSorted = T)
  insertTable(data.table(stopCodon = as.character(seqs, use.names = F)), "StopCodons")
  # Stop codon grouping
  insertTable(data.table(stopCodonGrouping = uniqueOrder(stops)), "stopCodonGrouping")
  # exon-exon junctions
  eej <- numExonsPerGroup(fiveUTRs, T)
  link <- readTable("linkORFsToTx")
  eej <- as.integer(eej[link$txNames])
  insertTable(data.table(eej = eej), "exon-exonJunctionsLeader")
  eejuORF <- numExonsPerGroup(grl, )
  insertTable(data.table(eejuORF = eejuORF), "exon-exonJunctionsuORFs")
  # gc content
  gc <- gcContent(grl, fa)
  insertTable(data.table(gc = gc), "gcContent")

  # Gene information

  # gene transcript connections
  dt <- data.table(txNames = txNames(grl), geneNames = ORFik:::txNamesToGeneNames(txNames(grl), Gtf))
  insertTable(dt, "uORFTxToGene")
  # Gene to symbol
  insertTable(getAllORFGeneSymbols(dt$geneNames), "geneSymbols")
  # go Terms
  uorfGo <- getORFsGoTerms(dt$geneNames)
  insertTable(data.table(go = uorfGo), "goTerms")
  return(NULL)
}

#' Get Ribo-seq features
#'
#' Excluding ones that uses RNA-seq normalizations
#' @return NULL (features saved to database)
getGeneralRiboFeatures <- function(grl, cds, threeUTRs, tx, df.rfp) {

  startRegion <- startRegion(grl, tx, T, -3, 9)
  rfps <- list.files("/export/valenfs/data/processed_data/Ribo-seq/All_human_bedo/", full.names = TRUE)[1:2]
  allRiboFeatures <- foreach(RFPPath=rfps, .combine = 'cbind',
                             .export = c("grl", "fiveUTRs", "cds", "threeUTRs", "tx", "startRegion"), .packages = "ORFik") %dopar% {
     RFP <- fimport(RFPPath)
     return(ORFik:::allFeaturesHelper(grl, RFP, RNA = NULL, tx, fiveUTRs, cds , threeUTRs,
                                      faFile = NULL, riboStart = 26, riboStop = 34,
                                      sequenceFeatures = FALSE, grl.is.sorted = TRUE,
                                      weight.RFP = "score", weight.RNA = 1L,
                                      st = startRegion))
                             }

  # libs <-bplapply(rfps,
  #                 function(x, grl, fiveUTRs, threeUTRs, cds, startRegion) {
  #                   RFP <- fimport(x)
  #                   return(ORFik:::allFeaturesHelper(grl, RFP, RNA = NULL, tx, fiveUTRs, cds , threeUTRs,
  #                                                    faFile = NULL, riboStart = 26, riboStop = 34,
  #                                                    sequenceFeatures = FALSE, grl.is.sorted = TRUE,
  #                                                    weight.RFP = "score", weight.RNA = 1L,
  #                                                    st = startRegion))
  #                 }, grl = grl, fiveUTRs = fiveUTRs, threeUTRs = threeUTRs, cds = cds, startRegion = startRegion)

  for(f in unique(colnames(allRiboFeatures))) { # Create one table per feature in DB
    featu <-  allRiboFeatures[, which(colnames(allRiboFeatures) == f), with = F]
    colnames(featu) <- paste0(f,"_", 1:ncol(featu))
    insertTable(data.table(featu), f)
  }
}
