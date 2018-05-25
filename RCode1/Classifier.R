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
  draw.pairwise.venn(sum(uorfPred$Ovary), sum(uorfPred$brain),
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