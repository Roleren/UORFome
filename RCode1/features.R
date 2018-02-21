
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
  grl <- readTable("uorfsAsGRWithTx", asGR = T)
  gr <- unlist(grl, use.names = F)
  names(gr) <- gsub("_[0-9]*", "", names(gr))

  grl <- groupGRangesBy(gr, gr$names)
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
  
  dt <- ORFik:::allFeatures(grl = grl, RFP = RFPShifted, RNA = NULL,
                            fiveUTRs = fiveUTRs, cds = cds, tx = tx,
                            threeUTRs = threeUTRs, faFile = fa, riboStart = 26, riboStop = 34,
                            extension = 1000, orfFeatures = T,
                            cageFiveUTRs = cageFiveUTRs, includeNonVarying = F)
  save(dt,file = paste0("featureTablesTemp/dt_",i,".rdata"))
  return(i)
}

testtest <- function(){
  dt <- ORFik:::allFeatures(grl = grl, RFP = RFPShifted, RNA = NULL,
                            fiveUTRs = fiveUTRs, cds = cds, tx = tx,
                            threeUTRs = threeUTRs, faFile = fa, riboStart = 26, riboStop = 34,
                            extension = 1000, orfFeatures = T,
                            cageFiveUTRs = cageFiveUTRs, includeNonVarying = F)
  a <- ORFik:::RibosomeStallingScore(grl, RFPShifted)
  grl_len <- ORFik:::widthPerGroup(grl, FALSE)
  overlapGrl <- countOverlaps(grl, RFPShifted)
  stopCodons <- ORFik:::ORFStopCodons(grl, TRUE)
  overlapStop <- countOverlaps(stopCodons, RFP)
  
  rss <- ((overlapStop + 1) / 3) / ((overlapGrl + 1) / grl_len)
}
#' validate features using brain and hela to see that features
#'  seperate them
validateAllFeatures <- function(){
  
  #rm(list=ls())
  #setwd("/export/valenfs/projects/uORFome/RCode1/")
  #source("./DataBaseCreator.R")
  
  # goal: what is different between groups and within group
  # variance / aov
  # clustering, all features hclust(dists) , kmean, jaccard index
  # scale to normalize
  # which are important, pca, svd etc.
  # get uorf names , do for both, change standardCage
  uids1 <- uidFromCage()
  
  uids2 <- uidFromCage("./../DATA/CAGE/human/kidney%2c%20adult%2c%20pool1.CNhs10622.10017-101C8.hg38.nobarcode.ctss.bed.gz")
  
  
  # now merge, either all, or cross set
  
  dtUid1 <- data.table(uid = uids1, index = 1:length(uids1))
  dtUid2 <- data.table(uid = uids2, index = 1:length(uids2))
  # cross valid set
  dtMerged <- merge(x = dtUid1, y = dtUid2, by = "uid")
  # all
  #dtMerged <- merge(x = dtUid1, y = dtUid2, by = "uid", all = T)
  
  # now get data
  experiments <- grep(x = list.files(getwd()), pattern = "dt_", value = T)
  experiments <- paste0(getwd(),"/", experiments)
  for(i in experiments){
    assign(x = paste0("dt", i), data.table::fread(i), envir = .GlobalEnv)
  }
  
  # fix na, inf, nan values, set to 0
  dtExper <- paste0("dt", experiments)
  for(i in dtExper){
    invisible(lapply(names(get(i)),function(.name) set(get(i), which(!is.finite(get(i)[[.name]])), j = .name,value =0)))
  }
  # now use merge to make nrows equal for both lists
  for(i in dtExper[1:5]){
    assign(i, get(i)[dtMerged$index.x,], .GlobalEnv)
    invisible(lapply(names(get(i)),function(.name) set(get(i), which(!is.finite(get(i)[[.name]])), j = .name,value =0)))
  }
  for(i in dtExper[6:10]){
    assign(i, get(i)[dtMerged$index.y,], .GlobalEnv)
    invisible(lapply(names(get(i)),function(.name) set(get(i), which(!is.finite(get(i)[[.name]])), j = .name,value =0)))
  }
  
  
  # get a data.table for each feature
  ncols <- ncol(get(dtExper[1]))
  nrows <- nrow(get(dtExper[1]))
  names <- colnames(get(dtExper[1]))
  
  for(i in names){
    assign(i, data.table(matrix(nrow = nrows, ncol = 10)), .GlobalEnv)
  }
  colnames <- colnames(get(i))
  # for each column data.table(ig. floss data.table), update each column
  x <- 1L
  for(i in names){
    
    currentCol <- get(i)
    y <- 1L
    for(j in dtExper){
      currentCol[, colnames[y] :=  get(j)[, x, with = F]]
      y <- y + 1L
    }
    assign(names[x], currentCol, .GlobalEnv)
    x <- x + 1L
  }
  # scale them ?
  
  # now do pca:
  for(i in names[2:length(names)]){
    assign(paste0("pca_",i), prcomp(get(i), scale = T), .GlobalEnv)
  }
  
  # make the pca, do the plots, make quantiles,
  # find overlapping uorfs between quantile sets 
  # a way to find interesting uorfs
  pcaNames <- paste0("pca_",names[2:length(names)])
  for(i in pcaNames){
    pca <- get(i)
    barplot(pca$sdev/sum(pca$sdev), main = i) # variance in percentage
    
    # split tissue groups
    plot(pca$rotation, col = c(rep(1, 5), rep(2, 5)), main = i)
    abline(h = mean(pca$rotation[,2]))
    abline(v = mean(pca$rotation[,1]))
    
    seperator = pca$x[,1]
    nfquant <- quantile(seperator, 0.95)
    zfquant <- quantile(seperator, 0.05)
    
    fivepercentSet <- seperator <= zfquant
    ninetyFivePercentSet <- seperator >= nfquant
    
    seperator = pca$x[, 2]
    nfquant <- quantile(seperator, 0.95)
    zfquant <- quantile(seperator, 0.05)
    
    fivepercentSet <- fivepercentSet | (seperator <= zfquant)
    ninetyFivePercentSet <- ninetyFivePercentSet | (seperator >= nfquant)
    
    fivepercentWhich <- which(fivepercentSet)
    ninetyFiveWhich <- which(ninetyFivePercentSet)
    assign(paste0("fq", i), fivepercentWhich, .GlobalEnv)
    assign(paste0("nfq", i), ninetyFiveWhich, .GlobalEnv)
  }
  # find the quntile overlaps of uorfs, all duplicated
  # these are the potentialy interesting uorfs
  # look for overlaps between the quantile sets
  # remove booleans
  # which features are good, which uorfs are regulated
  pcaNamesNonBool <- pcaNames[-c(9, 11, 12, 13, 14)]
  
  nfqNames <- paste0("nfq", pcaNamesNonBool)
  
  fqNames <- paste0("fq", pcaNamesNonBool)
  tempMatches <- NULL
  for(i in 1:length(fqNames)){
    tempMatches <- c(tempMatches, get(fqNames[i]))
  }
  oriMatches <- tempMatches
  oriLengthMatches <- length(oriMatches)
  fquniqueMatches <- unique(tempMatches)
  lengthUniqueMatches <- length(fquniqueMatches)
  fqdiffMatches <-lengthUniqueMatches / oriLengthMatches
  
  tempMatches <- NULL
  for(i in 1:length(nfqNames)){
    tempMatches <- c(tempMatches, get(nfqNames[i]))
  }
  oriMatches <- tempMatches
  oriLengthMatches <- length(oriMatches)
  nfquniqueMatches <- unique(tempMatches)
  lengthUniqueMatches <- length(nfquniqueMatches)
  nfqdiffMatches <-lengthUniqueMatches / oriLengthMatches
  
  inBoth <- fquniqueMatches[fquniqueMatches %in% nfquniqueMatches]
  
  onlyFirst <- fquniqueMatches[!(fquniqueMatches %in% nfquniqueMatches)]
  onlySecond <- nfquniqueMatches[!(nfquniqueMatches %in% fquniqueMatches)]
  
  
  # for FPKM RFP
  View(fpkmRFP[onlyFirst])
  View(fpkmRFP[onlySecond])
  meanFirst <- mean(colMeans(fpkmRFP[onlyFirst]))
  meanSecond <- mean(colMeans(fpkmRFP[onlySecond]))
  # for FPKM RNA
  meanFirstRNA <- mean(colMeans(fpkmRNA[onlyFirst]))
  meanSecondRNA <- mean(colMeans(fpkmRNA[onlySecond]))
  # the ones that overlap, check
  
  #for ioscore
  meanFirstIO <- mean(colMeans(ioScore[onlyFirst]))
  meanSecondIO <- mean(colMeans(ioScore[onlySecond]))
  
  #for orfsScore
  meanFirstORFScore <- mean(colMeans(ORFScores[onlyFirst]))
  meanSecondORFScore <- mean(colMeans(ORFScores[onlySecond]))
  
  #for fractionLengths
  meanFirstfractionLengths <- mean(colMeans(fractionLengths[onlyFirst]))
  meanSecondfractionLengths <- mean(colMeans(fractionLengths[onlySecond]))
  
  #for te
  meanFirstTe <- mean(colMeans(te[onlyFirst]))
  meanSecondTe <- mean(colMeans(te[onlySecond]))
  nInBoth <- sum(fquniqueMatches %in% nfquniqueMatches)
  print(paste(round(nInBoth/max(length(fquniqueMatches),
                                       length(fquniqueMatches)), 2)*100, 
                "% of uorfs are in both quantile sets"))
  
  # cluster
  for(i in names){
    clusterUorfFeature(get(i), inBoth, paste0("./clustering/","cluster_",i,".pdf"))
  }
 
}

#' Test features
#' 
#' Just if something is wrong
test_all_features <- function(){
  #these are slow:
  #1. floss
  #2. entropy
  
  #rm(list=ls())
  #setwd("/export/valenfs/projects/uORFome/RCode1/")
  #source("./DataBaseCreator.R")
  #detach("package:ORFik", unload=TRUE)
  getCDS()
  getThreeUTRs()
  getLeaders()
  tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)
  
  #originalUorfsByTx <- getUnfilteredUORFs(fiveUTRs, assignRanges = F)
  
  
  #or with extension
  cageFiveUTRs <- ORFik:::reassignTSSbyCage(fiveUTRs, standardCage, 1000, 1, cds)
  originalUorfsByTx <- getUnfilteredUORFs(cageFiveUTRs, assignRanges = F)
  tx <- ORFik:::extendLeaders(tx, extension = cageFiveUTRs)
  RFPShifted <- ORFik:::cageFromFile("/export/valenfs/data/processed_data/Ribo-seq/fantom_human_bed/Gonzalez_C_2014.Human.brain.RPF.GRCh38.SRR1562543.reads_merged.bed")
  gr <- unlist(originalUorfsByTx, use.names = F)
  grl <- groupGRangesBy(gr, gr$names)
  grl <- removeORFsWithinCDS(grl)
  #tx_len <- ORFik:::widthPerGroup(tx)
  #RFP <- readGAlignments("/export/valenfs/data/processed_data/Ribo-seq/gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/Gonzalez_C_2014.Human.brain.RPF.GRCh38.SRR1562540.bam")
  RNA <- readGAlignments("/export/valenfs/data/processed_data/RNA-seq/gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/Gonzalez_C_2014.Human.brain.RNA.GRCh38.SRR1562546.bam")
  
  # floss <- ORFik:::floss(grl, RFP, cds)
  # entropy <- ORFik:::entropy(grl, RFP)
  # tedt <- ORFik:::te(grl, RNA, RFP, tx_len, with.fpkm = T)
  # fracLengths <- ORFik:::fractionLength(grl, tx_len)
  # dScore <- ORFik:::disengagementScore(grl, RFP, tx)
  # releaseScore <- ORFik:::RibosomeReleaseScore(grl, RFP, threeUTRs, RNA)
  # orfScore <- ORFik:::ORFScores(grl, RFP)
  # distOrfToCds <- ORFik:::distOrfToCds(grl, fiveUTRs, extension = 0)
  # kozak <- ORFik:::kozakSequenceScore(grl, fa)
  # ioScore <- ORFik:::insideOutsideORF(grl, RFP, tx)
  # frame <- ORFik:::inFrameWithCDS(distOrfToCds)
  # overlapping <- ORFik:::isOverlappingCds(distOrfToCds)
  # ranks <- ORFik:::OrfRankOrder(grl)
  # ioScore <- ORFik:::insideOutsideORF(grl, RFP, tx)
  # dt <- ORFik:::allFeatures(grl = grl, RFP = RFP, RNA = RNA,
  #                            fiveUTRs = fiveUTRs, cds = cds, tx = tx,
  #                           threeUTRs = threeUTRs, faFile = fa, riboStart = 26, riboStop = 34,
  #                           extension = 0, orfFeatures = T)
  #or
  dt <- ORFik:::allFeatures(grl = grl, RFP = RFPShifted, RNA = RNA,
                            fiveUTRs = fiveUTRs, cds = cds, tx = tx,
                            threeUTRs = threeUTRs, faFile = fa, riboStart = 26, riboStop = 34,
                            extension = 1000, orfFeatures = T,
                            cageFiveUTRs = cageFiveUTRs)
}
