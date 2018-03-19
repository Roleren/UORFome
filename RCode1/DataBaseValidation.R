
#' Validate uorfs of data-base
#' 
#' This is a check to see that pipeline have done everything correctly
#' if redoing the findOverlaps does not find all orfs within fiveUTRs
#' it means that some orfs are outside the mapping area
#' this should not happen!
validateExperiments <- function(grl){
  
  fiveUTRs <- leaderAllSpanning()
  a <- findOverlaps(query = unlist(grl, use.names = F), fiveUTRs)
  a <- a[!duplicated(from(a))]
  if(length(a) != length(unlist(grl))){ 
    stop("Not all orfs where within the FiveUTRs used
         to make them, something is wrong!")
  } else { print("experiments look valid")}
}

#' Get variance between different leader versions
#' 
#' For validation
getAllLeaderChanges <- function(){
  if(!file.exists(p(dataFolder,"/leaderOriginalWidths.rdata"))){
    getLeaders()
    widths <- ORFik:::widthPerGroup(fiveUTRs)
    save(widths, file = p(dataFolder,"/leaderOriginalWidths.rdata"))
    rm(fiveUTRs)
  }
  
  library(doParallel)
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  maxCores = as.integer(detectCores()/2)
  cl <- makeCluster(maxCores)
  registerDoParallel(cl)
  leadersList = list.files(leadersFolder)
  nLeadersList = length(leadersList)
  rm(fiveUTRs)
  output <- foreach(i=1:nLeadersList, .combine = 'rbind') %dopar% {
    source("./uorfomeGeneratorHelperFunctions.R")
    leadersList = list.files(leadersFolder)
    
    load(p(dataFolder,"/leaderOriginalWidths.rdata"))
    load(p(leadersFolder,leadersList[i]))
    widthsCage <- ORFik:::widthPerGroup(fiveUTRs)
    
    diffWidths <- widths - widthsCage
    same <- sum(diffWidths == 0)
    bigger <- sum(diffWidths < 0)
    smaller <- sum(diffWidths > 0)
    meanDifBigger <- mean(diffWidths[diffWidths < 0])
    meanDifSmaller <- mean(diffWidths[diffWidths > 0])
    return(c(same,bigger,smaller,meanDifBigger,meanDifSmaller))
  }
  dt <- as.data.table(matrix(output, ncol = 5))
  colnames(dt) <- c("same", "bigger", "smaller", "meanBigger", "meanSmaller")
  stopCluster(cl)
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  save(dt,file = "leaderWidthChanges.rdata")
}

#' validate features using brain and hela to see that features
#'  seperate them
#'  using pca and quantiles
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

# look like a good features test, the cds does not seperate well
kozakVsRiboseq <- function(){
  kozak <- as.numeric(unlist(readTable("kozak", with.IDs = FALSE), use.names = FALSE))
  q99 <- quantile(kozak, 0.99)
  best <- which(kozak > q99)
  ribo <- readTable("Ribofpkm", with.IDs = FALSE)
  bestMean <- mean(colMeans(ribo[best,]))
  worst <- which(kozak < quantile(kozak, 0.01))
  worstMean <- mean(colMeans(ribo[worst,]))
  
  getCDS()
  
  cdsFPKM <- readTable("cdsRfpFPKMs", with.IDs = F)
  kozakCDS <- kozakSequenceScore(cds, fa)
  bestCDS <-  which(kozakCDS > quantile(kozakCDS, 0.90))
  bestMeanCDS <- mean(colMeans(cdsFPKM[bestCDS,]))
  worstCDS <-  which(kozakCDS < quantile(kozakCDS, 0.10))
  worstMeanCDS <- mean(colMeans(cdsFPKM[worstCDS,]))
}

# Looks like the distance increases with better kozak sequence
kozakVsdistanceToCDS <- function(){
  kozak <- as.numeric(unlist(readTable("kozak", with.IDs = FALSE), use.names = FALSE))
  q99 <- quantile(kozak, 0.99)
  best <- which(kozak > q99)
  distCDS <- readTable("distORFCDS", with.IDs = FALSE)
  bestMean <- mean(colMeans(distCDS[best,]))
  worst <- which(kozak < quantile(kozak, 0.01))
  worstMean <- mean(colMeans(distCDS[worst,]))
}

# Looks like a good feature test, must check cds
kozakVsORFScores <- function(){
  kozak <- as.numeric(unlist(readTable("kozak", with.IDs = FALSE), use.names = FALSE))
  q99 <- quantile(kozak, 0.99)
  best <- which(kozak > q99)
  ORFScores <- readTable("ORFScores", with.IDs = FALSE)
  bestMean <- mean(colMeans(ORFScores[best,]))
  worst <- which(kozak < quantile(kozak, 0.01))
  worstMean <- mean(colMeans(ORFScores[worst,]))
}

#' How much does the TE go down for CDS with uorfs in tx
uorfTeVsCDSTe <- function(){
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  grl <- getUorfsInDb()
  uorfNames <- unique(OrfToTxNames(grl))
  uorfTXCDS <- cdsTEs$txNames %in% uorfNames
  cdsTEUORFs <- cdsTEs[uorfTXCDS, 2:ncol(cdsTEs)]
  withUorfCDS <- sum(colMeans(cdsTEUORFs))
  withoutUorfCDS <- sum(colMeans(cdsTEs[!uorfTXCDS, 2:ncol(cdsTEs)]))
  
  # quantile
  cdsTxUorfs <- OrfToTxNames(grl) %in% cdsTEs$txNames
  uorfTEs <- readTable("teFiltered", with.IDs = F)
  
  
  # uorfTEs <- data.table(riboFPKM[,1:2], uorfTEs)
  uorfTEs <- uorfTEs[, lapply(.SD, mean), by = OrfToTxNames(grl)]
  colnames(uorfTEs)[1] <- "txNames"
  # best
  rowMeansUorfs <- rowMeans(uorfTEs[,2:ncol(uorfTEs)])
  q90 <- quantile(rowMeansUorfs, 0.90)
  best <- which(rowMeansUorfs > q90)
  bestCDSTEs <- sum(colMeans(cdsTEUORFs[best, ]))
  #worst
  q10 <- quantile(rowMeansUorfs, 0.10)
  worst <- which(rowMeansUorfs < q10)
  worstCDSTEs <- sum(colMeans(cdsTEUORFs[worst,]))
  # conclusion, no q90/q10 correlation it looks like
  
  rowMeansCDS <- rowMeans(cdsTEs[uorfTXCDS, 2:ncol(cdsTEs)])
  
  corResult <- cor.test(rowMeansUorfs, rowMeansCDS)
  
  
  return((1 - withoutUorfCDS/withUorfCDS)*100)
}

distVSCDSTe <- function(){
  
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  dists <- readTable("distORFCDS")
  riboFPKM <- readTable("Ribofpkm")
  
  q90 <- quantile(dists$distORFCDS, 0.90)
  best <- which(dists$distORFCDS > q90)
  bestDists <- dists[best,]
  
  longestTXNames <- cdsTEs$txNames %in% riboFPKM$txNames[riboFPKM$uorfID %in% bestDists$uorfID]
  longestCDSTe <- sum(colMeans(cdsTEs[longestTXNames, 2:ncol(cdsTEs)]))
  
  q10 <- quantile(dists$distORFCDS, 0.10)
  worst <- which(dists$distORFCDS < q10)
  worstDists <- dists[worst,]
  
  shortestTXNames <- cdsTEs$txNames %in% riboFPKM$txNames[riboFPKM$uorfID %in% worstDists$uorfID]
  shortestCDSTe <- sum(colMeans(cdsTEs[shortestTXNames, 2:ncol(cdsTEs)]))
}
