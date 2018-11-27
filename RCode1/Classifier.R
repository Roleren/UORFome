#' random forrest classification:
#' Training data
#' pos set is cds
#' Neg seg is 3'utrs
#' @param tissues Tissues to train on, use all if you want all in one
predictUorfs <- function(tissues = c("brain","fibroblasts", "kidney", "prostate", "Ovary"),
                         nthreads = 40,  max_mem_size = "60G") {
  # pick tissue
  
  pipelineCluster(8)
  if(!file.exists( paste0("forests/finalPrediction_",tissue, ".rdata"))) {
    print(paste("running for tissue:", tissue))
    
    #testForrest(predicate, tissue)
    print("starting riboseq classifier for uORFs")
    
    
    # make uORFTable
    if(file.exists(paste0("forests/prediction_", tissue, ".rdata"))) {
      load(paste0("forests/prediction_", tissue, ".rdata"))
    } else {
      forestRibo <- trainClassifier(tissue)
      uorfTable <- makeUORFPredicateTable()
      # remove cds overlaps
      uorfTable <- uorfTable[!overCDS()]
      
      prediction <- as.data.table(h2o.predict(forestRibo,  as.h2o(uorfTable)))
      save(prediction, file = paste0("forests/prediction_", tissue, ".rdata"))
    }
      
      sequenceClassifier(prediction, tissue)
    }
    
  makeCombinedPrediction(tissues)
  
  # make some plots here on ribo seq prediction
  # Plot predictied sequences features from the ribo prediction mapped to uORFs
  # check some examples
}

#' The training model with cds and 3' UTRs as random forest
#' @param tissue Tissues to train on, use all if you want all in one
trainClassifier <- function(tissue = NULL) {
  
  if(file.exists(paste0("forests/randomForrest_",tissue))) {
    forestRibo <- h2o.loadModel(path = paste0("forests/randomForrest_",tissue,"/",
                                          list.files(paste0("forests/randomForrest_",tissue)[1])))
    return(forest)
  }
  predictors <- makePredicateTable(tissue)
  
  # define training control
  forestRibo <- forest(predictors, ntrees = 100)
  if(!is.null(tissue)) {
    h2o.saveModel(forestRibo, path = paste0("forests/randomForrest_",tissue))
  }
  return(forestRibo)
}

#' ORF sequence classifier
#' @param prediction data.table of predictions
#' @param tissue Tissues to train on, use all if you want all in one
sequenceClassifier <- function(prediction, tissue){
  print("predicting sequence classifier")
  
  uorfData <- getAllSequenceFeaturesTable()
  
  dt <- uorfData[!overCDS()]
  dt[,y := as.factor(prediction$p1 > 0.75)]
  
  # table(dt$StartCodons)
  # make classification
  forestH2o <- forest(dt, cv = 5, ntrees = 150)
  print(forestH2o@model)
  h2o.saveModel(forestH2o, path = paste0("forests/finalForest_",tissue))
  # prediction
  uorfPrediction <- as.data.table(h2o.predict(forestH2o, newdata = as.h2o(uorfData)))
  
  save(uorfPrediction, file = paste0("forests/finalPrediction_",tissue, ".rdata"))
  
  # checking
  hits <- as.logical(uorfPrediction$p1 > 0.75)
  startCodonMetrics(hits)
  return(NULL)
}