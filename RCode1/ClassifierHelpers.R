
#' Get cds and 3'UTR table of ribo seq features
#' 
#' Positive set is cds, negative is 3' UTRs
#' @param tissue Tissue to train on, use all if you want all in one
makePredicateTable <- function(tissue) {
  if(file.exists(paste0("forests/predicateTables/table_cds3utr_",tissue,".rdata"))) {
    load(paste0("forests/predicateTables/table_cds3utr_",tissue,".rdata"))
    return(predicate)
  }
  # group all features by tissue
  # pick 1 tissue
  # make the predicate table
  
  posFeatureNames <- grep(pattern = "cds", x = listTables(), value = T)
  # posFeatureNames <- posFeatureNames[-grep(pattern = "Unfiltered", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "Tissue", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "Kozak", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "FractionLengths", x = posFeatureNames, value = F)]
  # posFeatureNames <- posFeatureNames[-grep(pattern = "FPKM", x = posFeatureNames, value = F)]
  
  negFeatureNames <- grep(pattern = "three", x = listTables(), value = T)
  # negFeatureNames <- negFeatureNames[-grep(pattern = "Unfiltered", x = negFeatureNames, value = F)]
  # negFeatureNames <- negFeatureNames[-grep(pattern = "FPKM", x = negFeatureNames, value = F)]
  
  if ( length(posFeatureNames) != length(negFeatureNames)) stop("Not equal length of pos and neg feature names")
  
  pos <- foreach(i = 1:length(posFeatureNames), .combine = 'cbind', .noexport = "uorfDB", .export = c("codeFolder")) %dopar% {
    setwd(codeFolder)
    source("./DataBaseSetup.R")
    
    getTissueFromFeatureTable(tableName = posFeatureNames[i], tissue = tissue)
  }
  validCDSByRiboSeq <- which(pos[,6] > 0.0189589) #FPKM
  pos <- pos[validCDSByRiboSeq, ]
  
  neg <- foreach(i = 1:length(negFeatureNames), .combine = 'cbind', .noexport = "uorfDB", .export = c("codeFolder")) %dopar% {
    setwd(codeFolder)
    source("./DataBaseSetup.R")
    
    getTissueFromFeatureTable(tableName = negFeatureNames[i], tissue = tissue)
  }
  predicate <- data.table(rbind(pos, neg))
  colnames(predicate) <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
                           "ORFScores","RFPFpkm", "RRS", "RSS", "startCodonCoverage")
  
  predicate <- fixNAandINFTable(predicate)
  y <- as.factor(c(rep(1, nrow(pos)), rep(0, nrow(neg))))
  predicate <- data.table(y, predicate)
  
  save(predicate, file = paste0("forests/predicateTables/table_cds3utr_",tissue,".rdata"))
  return(predicate)
}

#' Predict table for uORFs
#' 
#' Get ribo-seq features for uORFs
makeUORFPredicateTable <- function(tissue) {
  if(file.exists(paste0("forests/predicateTables/table_",tissue,".rdata"))) {
    load(paste0("forests/predicateTables/table_",tissue,".rdata"))
    return(predictors)
  }
  
  # group all features by tissue
  # pick 1 tissue
  # make the predicate table
  # featureNames <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
  #                   "ORFScores", "RRS", "RSS", "teFiltered")
  featureNames <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
                    "ORFScores","Ribofpkm", "RRS", "RSS", "startCodonCoverage")
  
  
  pos <- foreach(i = 1:length(featureNames), .combine = 'cbind', .noexport = "uorfDB", .export = "featureNames") %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path
    source("./DataBaseSetup.R")
    
    getTissueFromFeatureTable(tableName = featureNames[i], tissue = tissue)
  }
  
  predicate <- data.table(pos)
  colnames(predicate) <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
                           "ORFScores","RFPFpkm", "RRS", "RSS", "startCodonCoverage")
  predictors <- fixNAandINFTable(predicate)
  save(predictors, file = paste0("forests/predicateTables/table_",tissue,".rdata"))
  
  return(predictors)
}

#' Get all sequence features for ORFs 
#' 
#' Extracted from the data base
getAllSequenceFeaturesTable <- function(){
  if (file.exists("forests/uORFSequenceFeatures.rdata")) {
    load(file = "forests/uORFSequenceFeatures.rdata")
  } else {
    # load features
    rankInTx <- readTable("rankInTx", with.IDs = F)
    
    kozak <- readTable("kozak", with.IDs = F)
    fractionLengths <- readTable("fractionLengths", with.IDs = F)
    StartCodons <- readTable("StartCodons")[,2]
    StopCodons <- readTable("StopCodons")[,2]
    distORFCDS <- readTable("distORFCDS", with.IDs = F)
    distORFTSS <- readTable("distORFTSS", with.IDs = F)
    inFrameCDS <- readTable("inFrameCDS", with.IDs = F)
    isOverlappingCds <- readTable("isOverlappingCds", with.IDs = F)
    grl <- getUorfsInDb(T, T, T)
    numberOfUorfsPerTx <- readTable("numberOfUorfsPerTx", with.IDs = T)
    order <- chmatch(txNames(grl), numberOfUorfsPerTx$txNames)
    numberOfUorfsPerTx <- numberOfUorfsPerTx$txNames[order]
    lengths <- widthPerGroup(grl, F)
    exonExonJunctions  <- readTable("exon-exonJunctionsuORFs", with.IDs = F)
    
    # make Predicates:
    uorfData <- data.table(rankInTx, lengths, kozak, fractionLengths,
                           StartCodons = as.factor(StartCodons$startCodon),
                           StopCodons = as.factor(StopCodons$stopCodon),
                           distORFCDS, distORFTSS, inFrameCDS, isOverlappingCds, 
                           exonExonJunctions)
    save(uorfData, file = "forests/uORFSequenceFeatures.rdata")
  }
  return(uorfData)
}

fixNAandINFTable <- function(predicate){
  if(!all((c("RSS", "entropyRFP") %in% colnames(predicate))))
    stop("predicate did not contain names RSS and entropyRFP")
  
  isNA <- is.na(predicate[,"RSS"])[,1] #RSS
  predicate <- as.matrix(predicate)
  isINF <- is.infinite(predicate[, "entropyRFP"])
  predicate[isNA, "RSS"] <- 0
  predicate[isINF, "entropyRFP"] <- 0
  return(as.data.table(predicate))
}

#' A filter per stop codon group
#' Uses 3 large ribo-seq libraries for ribo-seq validation
getBestIsoformStartCodonCoverage <- function() {
  if(file.exists(paste0("forests/bestStartCodons.rdata"))) {
    load(paste0("forests/bestStartCodons.rdata"))
    return(d)
  }
  # reduce isoform overlaps by highest start codon reads per group
  g <- getUorfsInDb(T, T, T)
  
  sg <- stopCodons(g, is.sorted = T)
  uo <- uniqueOrder(sg) # <- grouping 
  
  # get start for each in group
  # count overlaps
  # return orf with highest per group
  # which start codon does it have ?
  counts <- rowMeans(readTable("startCodonCoverage"))
  d <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))
  save(d, file = paste0("forests/bestStartCodons.rdata"))
  return(d)
}

#' Train h2o rf model.
#' negDT if you want own samples for that
forest <- function(dt, cv = 10, ntrees = 50){
  # new h2o model training
  indices <- sample(1:nrow(dt), 0.6*nrow(dt), replace = F)
  validationIndices <- seq.int(nrow(dt))[-indices]
  trainingFrame <-  as.h2o(dt[indices, ])
  validationFrame <- as.h2o(dt[validationIndices,])
  forestH2o <- h2o.randomForest(y = "y", training_frame = trainingFrame,
                                validation_frame = validationFrame,
                                nfolds = cv, ntrees = ntrees)
  print(data.table(feature = forestH2o@model$variable_importances[,1],
                   round(forestH2o@model$variable_importances[,3:4], 2)))
  print(data.frame(name = rownames(forestH2o@model$cross_validation_metrics_summary)[c(1,17,19)], 
                   score = round(as.double(forestH2o@model$cross_validation_metrics_summary[c(1,17,19),1]),
                                 2)))
  return(forestH2o)
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

#' Combine classifier and CAGE data, for final prediction table
#' 
makeCombinedPrediction <- function(tissue, cutOff = 0.75) {
  # load data
  load(paste0("forests/finalPrediction_",tissue, ".rdata"))
  cageTissues <- readTable("tissueAtlasByCage", with.IDs = F)
  
  some <- uorfPrediction$p1 > cutOff # set value
  
  cageTissuesPrediction <- copy(cageTissues)
  for(i in colnames(cageTissuesPrediction)) {
    cageTissuesPrediction[, paste(i) := (cageTissues[,i, with=F] & some)]
  }
  insertTable(cageTissuesPrediction, "tissueAtlasByCageAndPred", rmOld = T)
  
  # sums <- colSums(cageTissuesPrediction)
  # cageTissuesPrediction[, names(which(sums == 0)) := NULL]
  
  finalCagePred <- rowSums(cageTissuesPrediction) > 0
  insertTable(finalCagePred, "finalCAGEuORFPrediction", rmOld = T)
  
  hits <- finalCagePred
  # good: ATG, CTG, procaryote: GTG, TTG
  # bad: AAG AND AGG
  ySeq <- rep(0, nrow(uorfPrediction))
  ySeq[hits] <- 1
  StartResultsSequences <- chisq.test(table(data.frame(uorfData$StartCodons, prediction = as.factor(ySeq))))
  print(paste("number of uORFs predicted translated:", sum(hits)))
  print(round(StartResultsSequences$residuals,1))
  print(round(table(uorfData$StartCodons[hits])/sum(hits), 2))
  print(table(uorfData$StartCodons[hits]))
}

# Get only specific tissues, or all.
# Grouped by rowMeans
getTissueFromFeatureTable <- function(tableName, tissue) {
  if (tableNotExists(tableName)) stop(paste("table does not exist:",
                                            tableName))
  rpfSamples <- getRiboRNAInfoTable()
  #1. we have tissues in link
  
  uniqueTissues <- as.character(unique(rpfSamples$tissue))
  if (is.null(tissue) | (tissue == "all")){
    print("Grouping all together")
  } else if ((tissue %in% uniqueTissues)){
    indices <- rpfSamples$tissue == tissue
    riboTable <- riboTable[,indices, with = F]
  } else stop("tissue does not exist in db")
  
  riboTable <- readTable(tableName, with.IDs = F)
  
  return(rowMeans(riboTable))
}


