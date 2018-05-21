#' Positive set is cds, negative is 3' UTRs
makePredicateTable <- function(tissue) {
  # group all features by tissue
  # pick 1 tissue
  # make the predicate table
  
  posFeatureNames <- grep(pattern = "cds", x = listTables(), value = T)
  posFeatureNames <- posFeatureNames[-grep(pattern = "Unfiltered", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "Tissue", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "Kozak", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "FractionLengths", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "FPKM", x = posFeatureNames, value = F)]
  
  negFeatureNames <- grep(pattern = "three", x = listTables(), value = T)
  negFeatureNames <- negFeatureNames[-grep(pattern = "Unfiltered", x = negFeatureNames, value = F)]
  negFeatureNames <- negFeatureNames[-grep(pattern = "FPKM", x = negFeatureNames, value = F)]
  
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
  colnames(predicate)[2:ncol(predicate)] <- c("disengagementScores", "entropyRFP", "floss", "ioScore",
                                              "ORFScores", "RRS", "RSS", "teFiltered")
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
  
  
  pos <- foreach(i = 1:length(featureNames), .combine = 'cbind', .noexport = "uorfDB") %dopar% {
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

trainClassifier <- function(tissue = NULL) {
  
  if(file.exists(paste0("forests/randomForrest_",tissue))) {
    forest <- h2o.loadModel(path = paste0("forests/randomForrest_",tissue,"/",
                                          list.files(paste0("forests/randomForrest_",tissue)[1])))
    return(forest)
  }
  predicate <- makePredicateTable(tissue)
  
  # fix style of matrix
  y <- predicate[,1]
  predictors <- predicate[,-1]
  isNA <- is.na(predictors[,6])[,1]
  predictors <- as.matrix(predictors)
  isINF <- is.infinite(predictors[,2 ])
  predictors[isNA, 6 ] <- 0
  predictors[isINF, 2] <- 0
  predictors <- as.data.table(predictors)
  predictors$y <- y
  # define training control
  indices <- sample(1:nrow(predictors), 0.6*nrow(predictors), replace = F)
  validationIndices <- seq.int(nrow(predictors))[-indices]
  trainingTable <- as.h2o(predictors[indices, ])
  validationTable <- as.h2o(predictors[validationIndices, ])
  
  # train the model
  forest <- h2o.randomForest(y = "y", training_frame = trainingTable,
                             validation_frame = validationTable,
                             nfolds = 10)
  if(!is.null(tissue)) {
    h2o.saveModel(forest, path = paste0("forests/randomForrest_",tissue))
  }
  return(forest)
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
  
  link <- readTable("linkRnaRfp")
  rpfSamples <- link[link$Sample_Type == "RPF",]
  #1. we have tissues in link
  
  uniqueTissues <- as.character(unique(rpfSamples$Tissue_or_CellLine))
  
  if ( !(tissue %in% uniqueTissues)) stop("tissue does not exist in db")
  
  riboTable <- readTable(tableName, with.IDs = F)
  
  if (length(grep(pattern = "tefiltered", tableName, ignore.case = T)) > 0) {
    indices <- rpfSamples[(rpfSamples$Tissue_or_CellLine == tissue),
                          matching]
  } else {
    indices <- rpfSamples[(rpfSamples$Tissue_or_CellLine == tissue),
                          originalIndex]
  }
  
  riboColumns <- riboTable[,indices, with = F]
  return(rowMeans(riboColumns))
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