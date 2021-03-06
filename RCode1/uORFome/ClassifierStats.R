checkTopPred <- function(){
  # load(paste0("forests/prediction_", "all", ".rdata"))
  load(paste0("forests/finalPrediction_filtered",tissue, ".rdata"))
  uorfData <- getAllSequenceFeaturesTable()
  uorfTable <- makeUORFPredicateTable()
  grl <- getUorfsInDb()
  index <- which(uorfTable$startCodonPerGroupBest & uorfData$StartCodons == "AGG" & prediction$p1 >= 0.50)
  t <- 18
  for(i in index[t:(t+2)]){
    print(uorfData[i,])
    print(uorfTable[i,])
    print(grl[i])
  }
}

checkBadPred <- function(){
  # region check
  uorfData <- getAllSequenceFeaturesTable()
  agg <- uorfData$StartCodons == "AGG" & prediction$predict == 1
  
  starts <- startSites(grl[agg], T, T, T)
  startRegion <- windowPerGroup(starts, tx, 6, 9)
  
  seqs <- getSequencesFromFasta(startRegion, T)
  allOver <- c("CTG|ATG|ACG|GTG|TTG")
  notBestStart <- (uorfTable$startCodonPerGroupBest == F)[agg]
  hits <- grep(x = seqs, pattern = allOver)
  valid <- (!(seq.int(1,sum(agg)) %in% hits))
  index <- which((!notBestStart & valid))
  t <- 81
  for(i in index[t:(t+2)]){
    print(i)
    print(prediction$predict[agg][i])
    print(uorfData[agg,][i])
    print(uorfTable[agg,][i])
    print(grl[agg][i])
  }
  toKeep <- which(agg)[index]
}

findClassificationBoundary <- function(tissues){
  x <- seq(0, 1, 0.1)[2:10]
  # for riboseq prediction
  for(tissue in tissues) {
    load(paste0("forests/prediction_",tissue, ".rdata"))
    hits <- unlist(lapply(x, function(y) sum(as.logical(prediction$p1 > y))),
                   use.names = F)
    plot(x, hits, main = tissue)
  }
  # we pick 0.5 from here
  # for sequence prediction
  for(tissue in tissues) {
    load(paste0("forests/finalPrediction_",tissue, ".rdata"))
    hits <- unlist(lapply(x, function(y) sum(as.logical(uorfPrediction$p1 > y))),
                   use.names = F)
    plot(x, hits, main = tissue,
         xlab = "prediction cutoff", ylab = "# of predicted uORFs")
  }
  # we pick 0.75 from here
  
  # validate boundaries
  for(j in seq(0.5, 1, 0.05)[2:10]) {
    print(paste("pred:", j))
    for(i in 1:5) {
      d <- getBestIsoformStartCodonCoverage()[readHits >= i & prediction$predict == 1 & prediction$p1 > j,]
      d <- d[, .SD[which.max(readHits)], by = group]
      print(i)
      print(round((table(uorfData$StartCodons[d$index])[6]/table(uorfData$StartCodons[d$index])[1]), 2))
    }
  }
  # validate read count
  
  start <- 0.5
  stop <- 1.5
  
  for(i in start:stop) {
    d <- getBestIsoformStartCodonCoverage()[readHits >= i & prediction$p1 >= 0.65,]
    d <- d[, .SD[which.max(readHits)], by = group]
    print(i)
    len[(i-start+1)] <- length(d$index)
  }
  plot(start:stop, len, main = paste("read count # change in", tissue),
       xlab = "# of reads as cutoff", ylab = "# of predicted uORF groups")
  
  x <- seq(0, 1, 0.05)[2:20]
  for(i in x){
    print(paste("x:   ",i))
    hits <- which(as.logical(uorfPrediction$p1 > i))
    # good: ATG, CTG, procaryote: GTG, TTG
    # bad: AAG AND AGG
    ySeq <- rep(0, nrow(uorfPrediction))
    ySeq[hits] <- 1
    StartResultsSequences <- chisq.test(table(data.frame(uorfData$StartCodons, prediction = as.factor(ySeq))))
    print(paste("number of uORFs predicted translated:", length(hits)))
    print(round(StartResultsSequences$residuals,1))
    print(round(table(uorfData$StartCodons[hits])/length(hits), 2))
    print("\n\n")
  }
  # We choose 0.75 here
}

#' Get start codon bias usage
#' good: ATG, CTG, ACG and GTG 
#' bad: AAG AND AGG
startCodonMetrics <- function(hits){
  if(!exists("uorfData")) uorfData <- getAllSequenceFeaturesTable()
  # if(!exists("uorfTable")) uorfTable <- 
  ySeq <- rep(0, length(hits))
  ySeq[hits] <- 1
  StartResultsSequences <- chisq.test(table(data.frame(uorfData$StartCodons, prediction = as.factor(ySeq))))
  # prin(cbind(a[order(a[,2], decreasing = T),], relative = a[order(a[,2], decreasing = T),2]/max(a[,2])))
  res <- round(StartResultsSequences$residuals,1)
  res <- res[order(res[,2],decreasing = T),] 
  count <- table(uorfData$StartCodons[hits])[rownames(res)]
  relativeCount <- round(table(uorfData$StartCodons[hits])/sum(hits), 2)[rownames(res)]
  res <- cbind(res, relative = round(res[,2]/max(res[,2]), 2), count, relativeCount)
  print(paste("number of uORFs predicted translated:", sum(hits)))
  print(res)

}

outputDivisionsOfPred <- function(prediction){
  print("orf score correlation with ribo prediction")
  print(cor.test(uorfTable$ORFScores, as.double(prediction$p1)))
  uorfTable <- makeUORFPredicateTable(tissue)
  # remove cds overlaps
  uorfTable <- makeUORFPredicateTable()
  uts <- uorfTable[!overCDS()]
  uorfData <- getAllSequenceFeaturesTable()
  uds <- uorfData[!overCDS()]
  
  pred <- prediction$predict
  print("Output feature summaries")
  print("For ribo seq prediction")
  print("On Ribo seq feature table")
  for(i in colnames(uorfTable)){
    print(i)
    print(summary(uts[pred == 1,i, with = F]))
    print(summary(uts[pred == 0,i, with = F]))
  }
  print("On sequence feature table")
  for(i in colnames(uorfData)){
    print(i)
    print(summary(uds[pred == 1,i, with = F]))
    print(summary(uds[pred == 0,i, with = F]))
  }
  print("For sequence prediction")
  print("On Ribo seq feature table")
  for(i in colnames(uorfTable)){
    print(i)
    print(summary(uorfTable[uorfPrediction$predict == T,i, with = F]))
    print(summary(uorfTable[uorfPrediction$predict == F,i, with = F]))
  }
  print("On sequence feature table")
  for(i in colnames(uorfData)){
    print(i)
    print(summary(uorfData[uorfPrediction$predict == T,i, with = F]))
    print(summary(uorfData[uorfPrediction$predict == F,i, with = F]))
  }
  
  # check orfscore vs start codon usage
  for(i in unique(uorfData$StartCodons)){
    print(paste(i , round(mean(uorfTable$ORFScores[uorfData$StartCodons == i]), 3)))
  }
  for(i in unique(uorfData$StartCodons)){
    print(paste(i , round(mean(uorfTable$entropyRFP[uorfData$StartCodons == i]), 3)))
  }
  for(i in unique(uorfData$StartCodons)){
    print(paste(i , round(mean(uorfTable$RRS[uorfData$StartCodons == i]), 3)))
  }
  for(i in unique(uorfData$StartCodons)){
    print(paste(i , round(mean(uorfTable$ioScore[uorfData$StartCodons == i]), 3)))
  }
  for(i in unique(uorfData$StartCodons)){
    print(paste(i , round(mean(uorfTable$RFPFpkm[uorfData$StartCodons == i]), 3)))
  }
  
  # more checks on orf score
  x <- seq(0,1, 0.05)[2:20]
  for(i in x){
    print(paste(i , round(mean(uorfTable$ORFScores[prediction$p1 > i]), 3)))
  }
  # cor test ORF score
  for(i in colnames(uorfTable)){
    print(i)
    print(cor.test(uorfTable$ORFScores, as.double(unlist(uorfTable[,i, with = F]))))
  }
}