#' random forrest classification:
#' Training data
#' pos set is cds
#' Neg seg is 3'utrs
#' @param tissues Tissues to train on, use all if you want all in one
predictUorfs <- function(tissues = c("brain","fibroblasts", "kidney", "prostate", "Ovary"),
                         nthreads = 40,  max_mem_size = "60G") {
  # pick tissue
  
  pipelineCluster(8)
  library(h2o)
  h2o.init(nthreads = 40, max_mem_size = "60G")
  
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

#' The training model with cds and 3' UTRs as random forest
#' @param tissue Tissues to train on, use all if you want all in one
trainClassifier <- function(tissue = NULL) {
  
  if(file.exists(paste0("forests/randomForrest_",tissue))) {
    forest <- h2o.loadModel(path = paste0("forests/randomForrest_",tissue,"/",
                                          list.files(paste0("forests/randomForrest_",tissue)[1])))
    return(forest)
  }
  predictors <- makePredicateTable(tissue)
  
  # define training control
  forest <- forest(predictors)
  if(!is.null(tissue)) {
    h2o.saveModel(forest, path = paste0("forests/randomForrest_",tissue))
  }
  return(forest)
}

#' ORF sequence classifier
#' @param prediction data.table of predictions
#' @param tissue Tissues to train on, use all if you want all in one
sequenceClassifier <- function(prediction, tissue){
  print("predicting sequence classifier")
  uorfData <- getAllSequenceFeaturesTable()
  
  dt <- uorfData
  # filter on isoforms
  d <- getBestIsoformStartCodonCoverage()
  # Here is lines with filter
  pred <- prediction$predict
  
  dPos <- d[readHits >= 0.9 & pred == 1 & prediction$p1 > 0.65,]
  dPos <- dPos[, .SD[which.max(readHits)], by = group]
  
  # combine filter with ribo-seq prediction
  dNeg <- d[readHits <= 0.1 & (pred == 0) & prediction$p0 > 0.65,]
  #dNeg <- dNeg[sample(x = dNeg$index, size = nrow(dPos)*2),]
  dNeg <- dNeg[, .SD[which.min(readHits)], by = group]
  dt[,y := as.factor(pred)]
  dt <- dt[c(dPos$index, dNeg$index),]
  # table(dt$StartCodons)
  # make classification
  forestH2o <- forest(dt, cv = 6, ntrees = 150)
  print(forestH2o@model)
  h2o.saveModel(forestH2o, path = paste0("forests/finalForest_",tissue))
  # prediction
  uorfPrediction <- h2o.predict(forestH2o, newdata = as.h2o(uorfData))
  uorfPrediction <- as.data.table(uorfPrediction)
  save(uorfPrediction, file = paste0("forests/finalPrediction_",tissue, ".rdata"))
  
  # checking
  hits <- which(as.logical(uorfPrediction$p1 > 0.75))
  # good: ATG, CTG, procaryote: GTG, TTG
  # bad: AAG AND AGG
  ySeq <- rep(0, nrow(uorfPrediction))
  ySeq[hits] <- 1
  StartResultsSequences <- chisq.test(table(data.frame(uorfData$StartCodons, prediction = as.factor(ySeq))))
  print(paste("number of uORFs predicted translated:", length(hits)))
  print(round(StartResultsSequences$residuals,1))
  print(round(table(uorfData$StartCodons[hits])/length(hits), 2))
  print(table(uorfData$StartCodons[hits]))
  return(NULL)
}