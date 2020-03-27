# get only sequence features from orfik
getSequenceFeatures <- function(){
  if(!tableNotExists("kozak") & !tableNotExists("gcContent")) {
    print("sequence features exist, stop if wrong!")
    return(NULL)
  }
  grl <- getUorfsInDb()
  getAll()
  orfID <- ORFik:::orfID(grl)

  # kozak
  kozak <- kozakSequenceScore(grl, tx, fa)
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
  # Stop codon grouping
  insertTable(data.table(stopCodonGrouping = uniqueOrder(stops)), "stopCodonGrouping")
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
  # Gene to symbol
  insertTable(getAllORFGeneSymbols(dt$geneNames), "geneSymbols")
  # go Terms
  uorfGo <- getORFsGoTerms(dt$geneNames)
  insertTable(data.table(go = uorfGo), "goTerms")
  # gc content
  getSequencesFromFasta(grl, T)
  alf <- alphabetFrequency(seqs, as.prob=TRUE)
  gc <- rowSums(alf[,c("G", "C")])
  insertTable(data.table(gc = gc), "gcContent")

  return(NULL)
}

#' Get Ribo-seq features
#'
#' Excluding ones that uses RNA-seq normalizations
getGeneralRiboFeatures <- function(grl, cds, threeUTRs, tx, name = "",
                                   rfpList = grep(pattern = "merged",x = list.files(rfpFolder), value = T)){

  nrfpList <- length(rfpList)
  startRegion <- windowPerGroup(startSites(grl, T, T, T), tx, -3, 9)

  allRiboFeatures <- foreach(i=1:nrfpList, .combine = 'cbind',
                             .export = c("grl", "cds", "threeUTRs", "tx", "startRegion", "rfpList", "name")) %dopar% {
     source(p(codeFolder,"/DataBaseSetup.R"))

     RFPPath <- p(rfpFolder, rfpList[i])
     RFP <- fread.bed(RFPPath)

     Diseng <- disengagementScore(grl, RFP, tx)
     Ios <- insideOutsideORF(grl, RFP, tx, ds = cdsDiseng)
     RRS <- ribosomeReleaseScore(grl, RFP, threeUTRs)
     RSS <- ribosomeStallingScore(grl, RFP)
     Floss <- floss(grl, RFP, cds)
     Entropy <- entropy(grl, RFP)
     ORFScores <- orfScore(grl, RFP, T)$ORFScores
     RfpFPKMs <- fpkm(grl, RFP)

     if(name == "") {
       startCodonCoverage <- ORFik:::startRegionCoverage(grl, RFP, tx)
     } else {
       startCodonCoverage <- countOverlaps(startCodons(grl, T), RFP)
     }
     Coverage <- countOverlaps(grl, RFP)
     fiveRegion <- countOverlaps(startRegion, fread.bed(RFPPath)) / pmax(widthPerGroup(startRegion), 1) / 3
     fiveRegionRelative <- startCodonCoverage/(fiveRegion + 1)
     return(data.table(Diseng, Ios, RRS, RSS, ORFScores, Floss, Entropy, RfpFPKMs, startCodonCoverage,
                       Coverage, fiveRegion, fiveRegionRelative))
   }
  colnames(allRiboFeatures) <- paste0(name, colnames(allRiboFeatures))
  txNames <- ORFik:::txNames(grl)
  for(f in unique(colnames(allRiboFeatures))) {
    featu <-  allRiboFeatures[, which(colnames(allRiboFeatures) == f), with = F]
    colnames(featu) <- paste0(f,"_", 1:ncol(featu))
    featu <- data.table(txNames, featu)
    insertTable(featu, f)
  }
}
