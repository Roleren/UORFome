#' Positive set is cds, negative is 3' UTRs
makePredicateTable <- function(tissue) {
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
  
  pos <- foreach(i = 1:length(posFeatureNames), .combine = 'cbind', .noexport = "uorfDB") %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path
    source("./DataBaseSetup.R")
    
    getTissueFromFeatureTable(tableName = posFeatureNames[i], tissue = tissue)
  }
  validCDSByRiboSeq <- getValidCDS()
  pos <- pos[validCDSByRiboSeq, ]
  
  
  neg <- foreach(i = 1:length(negFeatureNames), .combine = 'cbind', .noexport = "uorfDB") %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path
    source("./DataBaseSetup.R")
    
    getTissueFromFeatureTable(tableName = negFeatureNames[i], tissue = tissue)
  }
  
  
  y <- c(rep(1, nrow(pos)), rep(0, nrow(neg)))
  
  predicate <- data.table(y, rbind(pos, neg))
  # colnames(predicate)[2:ncol(predicate)] <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
  #                                             "ORFScores", "RRS", "RSS", "teFiltered")
  colnames(predicate)[2:ncol(predicate)] <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
                                              "ORFScores","RFPFpkm", "RRS", "RSS")
  return(predicate)
}

#' Predict table for uORFs
makeUORFPredicateTable <- function(tissue) {
  if(file.exists(paste0("forests/predicateTables/table_",tissue,".rdata"))) {
    load(paste0("forests/predicateTables/table_",tissue,".rdata"))
    return(predictors)
  }
  
  # group all features by tissue
  # pick 1 tissue
  # make the predicate table
  featureNames <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
                    "ORFScores", "RRS", "RSS", "teFiltered")
  featureNames <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
                    "ORFScores", "RRS", "RSS", "teFiltered")
  
  
  pos <- foreach(i = 1:length(featureNames), .combine = 'cbind', .noexport = "uorfDB", .export = "featureNames") %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path
    source("./DataBaseSetup.R")
    
    getTissueFromFeatureTable(tableName = featureNames[i], tissue = tissue)
  }
  
  predicate <- data.table(pos)
  colnames(predicate) <- featureNames
  isNA <- is.na(predicate[,6])[,1]
  predictors <- as.matrix(predicate)
  isINF <- is.infinite(predictors[, 2])
  predictors[isNA, 6] <- 0 # RRS
  predictors[isINF, 2] <- 0 # Entropy
  
  predictors <- as.data.table(predictors)
  save(predictors, file = paste0("forests/predicateTables/table_",tissue,".rdata"))
  
  return(predictors)
}

getValidCDS <- function() {
  if(!exists("validCDSByRiboSeq")) {
    if(file.exists(p(dataBaseFolder,"/validCDSByRiboSeq.rdata")) ) {
      load(p(dataBaseFolder,"/validCDSByRiboSeq.rdata"))
      return(validCDSByRiboSeq)
    }
  } else {
    return(validCDSByRiboSeq)
  }
  
  RFP <- readTable("cdsRfpFPKMs")
  validCDSByRiboSeq <- which(rowMeans(RFP[,2:ncol(RFP)]) > 0.0189589 )
  save(validCDSByRiboSeq, file = p(dataBaseFolder,"/validCDSByRiboSeq.rdata"))
  return(validCDSByRiboSeq)
}

getTissueFromFeatureTable <- function(tableName, tissue) {
  if (tableNotExists(tableName)) stop(paste("table does not exist:",
                                            tableName))
  rpfSamples <- getRiboRNAInfoTable()
  #1. we have tissues in link
  
  uniqueTissues <- as.character(unique(rpfSamples$tissue))
  
  if ( !(tissue %in% uniqueTissues)) stop("tissue does not exist in db")
  
  riboTable <- readTable(tableName, with.IDs = F)
  
  indices <- rpfSamples$tissue == tissue
  
  riboColumns <- riboTable[,indices, with = F]
  return(rowMeans(riboColumns))
}


findClassificationBoundary <- function(tissues){
  x <- seq(0, 1, 0.1)[2:11]
  # for riboseq prediction
  for(tissue in tissues) {
    load(paste0("forests/prediction_",tissue, ".rdata"))
    hits <- unlist(lapply(x, function(y) sum(as.logical(uorfPrediction > y))),
                   use.names = F)
    plot(x, hits, main = tissue)
  }
  # we pick 0.5 from here
  # for sequence prediction
  for(tissue in tissues) {
    load(paste0("forests/finalPrediction_",tissue, ".rdata"))
    hits <- unlist(lapply(x, function(y) sum(as.logical(uorfPrediction > y))),
                   use.names = F)
    plot(x, hits, main = tissue,
         xlab = "prediction cutoff", ylab = "# of predicted uORFs")
  }
  # we pick 0.69 from here
  
  # validate boundaries
  for(j in seq(0, 1, 0.1)[2:10]) {
    print(paste("pred:", j))
    for(i in 15:30) {
      d <- getBestIsoformStartCodonCoverage()[readHits >= i & pred$predict >= j,]
      d <- d[, .SD[which.max(readHits)], by = group]
      print(i)
      print(round((table(df$startCodon[d$index])[6]/table(df$startCodon[d$index])[1]), 2))
    }
  }
  # validate read count
  pred <- prediction[uniqueOrder, ]
  start <- 10
  stop <- 30
  len <- vector("numeric",length = (stop-start + 1))
  for(i in start:stop) {
    d <- getBestIsoformStartCodonCoverage()[readHits >= i & pred$predict >= 0.7,]
    d <- d[, .SD[which.max(readHits)], by = group]
    print(i)
    len[(i-start+1)] <- length(d$index)
  }
  plot(start:stop, len, main = paste("read count # change in", tissue),
       xlab = "# of reads as cutoff", ylab = "# of predicted uORF groups")
}


getBestIsoformStartCodonCoverage <- function(grl = NULL) {
  if(file.exists(paste0("forests/bestStartCodons.rdata"))) {
    load(paste0("forests/bestStartCodons.rdata"))
    return(d)
  }
  # reduce isoform overlaps by highest start codon reads per group
  g <- uniqueGroups(grl)
  
  sg <- stopCodons(g, is.sorted = T)
  usg <- uniqueGroups(sg)
  uo <- uniqueOrder(sg) # <- grouping 
  
  # get start for each in group
  # count overlaps
  # return orf with highest per group
  # which start codon does it have ?
  usg <- startCodons(g)
  RFP <- fread.bed(p(rfpFolder,list.files(rfpFolder)[97]))
  RFP2 <- fread.bed(p(rfpFolder,list.files(rfpFolder)[7]))
  RFP3 <- fread.bed(p(rfpFolder,list.files(rfpFolder)[14]))
  counts <- countOverlaps(usg, RFP) + countOverlaps(usg, RFP2) + countOverlaps(usg, RFP3)
  names(counts) <- NULL
  
  d <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))
  save(d, file = paste0("forests/bestStartCodons.rdata"))
  return(d)
}

makeCombinedPrediction <- function(tissues, cutOff = 0.7) {
  
  if (tableNotExists("uorfPredictions")) {
    for(tissue in tissues) {
      load(paste0("forests/finalPrediction_",tissue, ".rdata"))
      
      if(tissue == tissues[1]) {
        uorfPred <- as.logical(uorfPrediction > cutOff)
      } else {
        uorfPred <- cbind(uorfPred, as.logical(uorfPrediction > cutOff))
      }               
    }
    uorfPred <- as.data.table(uorfPred)
    colnames(uorfPred) <- tissues
    insertTable(uniqueUorfPred, "uorfPredictions")
  } else {
    readTable("uorfPredictions")
  }
  
  # tests
  tab <- table(uorfPred$Ovary, uorfPred$brain)
  chi <- chisq.test(tab)
  chi$residuals
  # plan:
  # do pairwise tests, see that things are ok.
  # venn diagram:
  
  library(VennDiagram)
  grid.newpage()
  boOver <- draw.pairwise.venn(sum(uorfPred$Ovary), sum(uorfPred$brain),
                               tab[2,2], category = c("Ovary", "Brain"),
                               lty = rep("blank", 2), fill = c("light blue", "yellow"),
                               alpha = rep(0.5, 2), cat.pos = c(0, 0),
                               cat.dist = rep(0.025, 2), title = "abc")
  # report feature difference on 0-1
  
  # with cage
  uniqueOrder <- readTable("uniqueOrder")$Matrix
  uniqueUorfPred <- uorfPred[uniqueOrder, ]
  #1. 
  cageTissues <- readTable("tissueAtlasByCage", with.IDs = F)
  brainRes <- cageTissues$brain & uniqueUorfPred$brain
  
  tab <- table(cageTissues$brain, uniqueUorfPred$brain)
  chi <- chisq.test(tab)
  chi$residuals
  
  grid.newpage()
  draw.pairwise.venn(sum(cageTissues$brain), sum(uniqueUorfPred$brain),
                     tab[2,2], category = c("Brain", "Brain"),
                     lty = rep("blank", 2), fill = c("light blue", "yellow"),
                     alpha = rep(0.5, 2), cat.pos = c(0, 0),
                     cat.dist = rep(0.025, 2), title = "abc")
  
  # filtering by cage
  tab <- table(uniqueUorfPred$Ovary & cageTissues$ovary,
               uniqueUorfPred$brain & cageTissues$brain)
  chi <- chisq.test(tab)
  chi$residuals
  grid.newpage()
  draw.pairwise.venn(sum(uniqueUorfPred$Ovary & cageTissues$ovary),
                     sum(uniqueUorfPred$brain & cageTissues$brain),
                     tab[2,2], category = c("Ovary", "Brain"),
                     lty = rep("blank", 2), fill = c("light blue", "purple"),
                     alpha = rep(0.5, 2), cat.pos = c(0, 0),
                     cat.dist = rep(0.025, 2), title = "Overlap of predicted uORFs and uORFs defined by CAGE.")
  
  
  
  uniqueUorfPred$all <- rowSums(uniqueUorfPred) == 5
  uniqueUorfPred$some <- rowSums(uniqueUorfPred) >= 2
  
  some <- rowSums(uniqueUorfPred) >= 2
  
  cageTissuesPrediction <- copy(cageTissues)
  for(i in colnames(cageTissuesPrediction)) {
    cageTissuesPrediction[, paste(i) := (cageTissues[,i, with=F] & some)]
  }
  
  
  insertTable(cageTissuesPrediction, "tissueAtlasByCageAndPred", rmOld = T)
  
  sums <- colSums(cageTissuesPrediction)
  cageTissuesPrediction[, names(which(sums == 0)) := NULL]
  sums <- colSums(cageTissuesPrediction)
  
  finalCagePred <- rowSums(cageTissuesPrediction) > 1
  insertTable(finalCagePred, "allUorfsByCageAndPred")
  
  grid.newpage()
  draw.triple.venn(area1 = 22, area2 = 20, area3 = 13, n12 = 11, n23 = 4, n13 = 5, 
                   n123 = 1, category = c("Ovary", "Brain"),
                   lty = rep("blank", 2), fill = c("light blue", "purple"),
                   alpha = rep(0.5, 2), cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2), title = "Overlap of predicted uORFs and uORFs defined by CAGE.")
}



# testForrest <- function(predicate, tissue) {
#   y <- factor(as.matrix(data.table(predicate[,1])))
#   predictors <- predicate[,-1]
#   isNA <- is.na(predictors[,6])[,1]
#   predictors <- as.matrix(predictors)
#   isINF <- is.infinite(predictors[,2 ])
#   predictors[isNA, 6 ] <- 0
#   predictors[isINF, 2] <- 0
#   
#   indices <- sample(1:nrow(predictors), 0.6*nrow(predictors), replace = F)
#   trainingSet <- predictors[indices, ]
#   trainingY <- y[indices]
#   testSet <- predictors[-indices, ]
#   testY <- y[-indices]
#   
#   forest <- randomForest(trainingSet, trainingY)
#   
#   prediction <- predict(forest, testSet)
#   
#   fp <- sum(prediction == 1 & testY == 0)
#   tp <- sum((prediction == 1) & (testY == 1))
#   
#   ppv <- tp / (tp + fp) 
#   
#   result <- sum(prediction == testY)/length(prediction)
#   print(paste("ppv was:", round(ppv, 3), "for tissue: ", tissue))
# }

# if(F) {
#   hits <- which(as.logical(uorfPrediction > 0.7))
#   StartCodons <- readTable("StartCodons")[,2]
#   uorfDataHits <- uorfData[hits, ]
#   
#   library(reshape2)
#   facet.df <- melt(facet.df, id.vars="prediction")
#   # lengths <- widthPerGroup(tx[txNames(grl)])
#   t <- data.frame(starts = StartCodons$startCodon, fractionLengths$fractionLengths)
#   
#   hp <- ggplot(t, aes(y=lengths, x = starts)) +
#     geom_boxplot() + scale_y_log10()
#   facet_wrap(~variable) +
#     theme(axis.text.x=element_text(angle=45))
#   plot(hp)
#   # grl <- getUorfsInDb()
#   
#   # check usage of sequence features
#   starts <- StartCodons$startCodon[uniqueOrder]
#   sMax <- starts[d$index][y == 1]
#   sNotMax <- starts[d$index][y == 0]
#   
#   groupsToCheck <- d$group[which(startsUnique[d$index] == "AAG")]
#   
#   pred <- as.data.table(prediction)
#   pred[pred > 0.5] <- 1
#   pred[pred <= 0.5] <- 0
#   startResults <- chisq.test(table(data.frame(StartCodons, prediction = as.factor(pred$predict))[uniqueOrder,]))
#   
#   ySeq <- rep(0, nrow(uorfPrediction))
#   ySeq[hits] <- 1
#   StartResultsSequences <- chisq.test(table(data.frame(StartCodons[uniqueOrder,], prediction = as.factor(ySeq)[uniqueOrder])))
#   StartResultsSequences$residuals
# }