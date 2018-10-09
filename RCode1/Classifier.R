#' random forrest classification:
#' pos set is cds
#' Neg seg is 3'utrs
predictUorfs <- function() {
  # pick tissue
  tissues <- c("brain","fibroblasts", "kidney", "prostate", "Ovary")
  library(foreach)
  library(h2o)
  h2o.init(nthreads = 40, max_mem_size = "50G")
  
  for( tissue in tissues) {
    if(!file.exists( paste0("forests/finalPrediction_",tissue, ".rdata"))) {
      print(paste("running for tissue:", tissue))
      
      #testForrest(predicate, tissue)
      print("starting riboseq classifier for uORFs")
      
      
      # make uORFTable
      if(file.exists(paste0("forests/prediction_", tissue, ".rdata"))) {
        load(paste0("forests/prediction_", tissue, ".rdata"))
      } else {
        forest <- trainClassifier(tissue)
        uorfTable <- makeUORFPredicateTable(tissue)
        prediction <- h2o.predict(forest,  as.h2o(uorfTable))
        prediction <- as.data.table(prediction)
        save(prediction, file = paste0("forests/prediction_", tissue, ".rdata"))
      }
      
      sequenceClassifier(prediction, tissue)
    }
    
    makeCombinedPrediction(tissues)
  }
  # make some plots here on ribo seq prediction
  # Plot predictied sequences features from the ribo prediction mapped to uORFs
  # check some examples
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


sequenceClassifier <- function(prediction, tissue){
  print("predicting sequence classifier")
  if (file.exists("forests/uORFSequenceFeatures.rdata")) {
    load(file = "forests/uORFSequenceFeatures.rdata")
    uniqueOrder <- readTable("uniqueOrder")$Matrix
  } else {
    # load features
    rankInTx <- readTable("rankInTx", with.IDs = F)
    numberOfUorfsPerTx <- readTable("numberOfUorfsPerTx", with.IDs = F)
    kozak <- readTable("kozak", with.IDs = F)
    fractionLengths <- readTable("fractionLengths", with.IDs = F)
    StartCodons <- readTable("StartCodons")[,2]
    StopCodons <- readTable("StopCodons")[,2]
    distORFCDS <- readTable("distORFCDS", with.IDs = F)
    inFrameCDS <- readTable("inFrameCDS", with.IDs = F)
    isOverlappingCds <- readTable("isOverlappingCds", with.IDs = F)
    grl <- getUorfsInDb()
    lengths <- widthPerGroup(grl, F)
    
    
    
    # make Predicates:
    uorfData <- data.table(rankInTx, lengths, kozak, fractionLengths,
                           StartCodons = as.factor(StartCodons$startCodon),
                           StopCodons = as.factor(StopCodons$stopCodon),
                           distORFCDS, inFrameCDS, isOverlappingCds)
    save(uorfData, file = "forests/uORFSequenceFeatures.rdata")
    
    uniqueOrder <- uniqueOrder(grl)
    uniqueOrder <- seq.int(length(grl))[!duplicated(uniqueOrder)]
    insertTable(uniqueOrder, "uniqueOrder")
  }
  
  dt <- uorfData
  dt <- dt[uniqueOrder, ]
  # filter on isoforms
  pred <- prediction[uniqueOrder,]
  d <- getBestIsoformStartCodonCoverage(grl)
  d <- d[readHits >= 20 & pred$predict >= 0.7,]
  d <- d[, .SD[which.max(readHits)], by = group]
  # combine filter with ribo-seq prediction
  

  negIndices <- sample((1:nrow(dt))[-d$index], 3*nrow(d), replace = F)
  # y <- y[d$index]
  y <- pred
  y[d$index,] <- 1
  y <- rbindlist(list(y[d$index,], y[-d$index,][negIndices,]))
  y[y < 1 ,] <- 0
  dt <- dt[c(d$index, negIndices),]
  dt$y <- y
  
  
  # new h2o model training
  indices <- sample(1:nrow(dt), 0.6*nrow(dt), replace = F)
  validationIndices <- seq.int(nrow(dt))[-indices]
  trainingFrame <-  as.h2o(dt[indices, ])
  validationFrame <- as.h2o(dt[validationIndices,])
  forestH2o <- h2o.randomForest(y = "y", training_frame = trainingFrame,
                                validation_frame = validationFrame,
                                nfolds = 15, balance_classes = F, ntrees = 150)
  print(forestH2o@model$variable_importances)
  h2o.saveModel(forestH2o, path = paste0("forests/finalForest_",tissue))
  # prediction
  uorfPrediction <- h2o.predict(forestH2o, newdata = as.h2o(uorfData))
  uorfPrediction <- as.data.table(uorfPrediction)
  save(uorfPrediction, file = paste0("forests/finalPrediction_",tissue, ".rdata"))
  
  # checking
  
  hits <- which(as.logical(uorfPrediction > 0.7))
  StartCodons <- readTable("StartCodons")[,2]
  ySeq <- rep(0, nrow(uorfPrediction))
  ySeq[hits] <- 1
  StartResultsSequences <- chisq.test(table(data.frame(StartCodons[uniqueOrder,], prediction = as.factor(ySeq)[uniqueOrder])))
  print(StartResultsSequences$residuals)
  return(NULL)
}