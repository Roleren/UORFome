
#' Get cds and 3'UTR table of ribo seq features
#'
#' Positive set is cds, negative is 3' UTRs
#' @param tissue Tissue to train on, use all if you want all in one
makePredicateTable <- function(tissue) {
  if (file.exists(paste0("forests/predicateTables/table_cds3utr_",tissue,".rdata"))) {
    load(paste0("forests/predicateTables/table_cds3utr_",tissue,".rdata"))
    return(predicate)
  }
  # group all features by tissue
  # pick 1 tissue
  # make the predicate table

  posFeatureNames <- grep(pattern = "cds", x = listTables(), value = T)
  posFeatureNames <- posFeatureNames[-grep(pattern = "Tissue", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "Kozak", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "FractionLengths", x = posFeatureNames, value = F)]
  posFeatureNames <- posFeatureNames[-grep(pattern = "cdsTE", x = posFeatureNames, value = F)]

  # posFeatureNames <- posFeatureNames[-grep(pattern = "FPKM", x = posFeatureNames, value = F)]

  negFeatureNames <- grep(pattern = "three", x = listTables(), value = T)
  # negFeatureNames <- negFeatureNames[-grep(pattern = "FPKM", x = negFeatureNames, value = F)]

  if ( length(posFeatureNames) != length(negFeatureNames)) stop("Not equal length of pos and neg feature names")

  pos <- foreach(i = 1:length(posFeatureNames), .combine = 'cbind', .noexport = "uorfDB") %dopar% {
    setwd(codeFolder)
    source("./DataBaseSetup.R")

    getTissueFromFeatureTable(tableName = posFeatureNames[i], tissue = tissue)
  }
  neg <- foreach(i = 1:length(negFeatureNames), .combine = 'cbind', .noexport = "uorfDB") %dopar% {
    setwd(codeFolder)
    source("./DataBaseSetup.R")

    getTissueFromFeatureTable(tableName = negFeatureNames[i], tissue = tissue)
  }

  # Filter out translating 3' UTRs
  filterThree <- (((neg[,1] < quantile(neg[,1], 0.95)) |
                     (neg[,12] < quantile(neg[,12], 0.981))) | neg[,5] < 1.1)
  # Filter out not needed columns of CDS
  filterCDS <- validByRibio(pos[,1], pos[,9], pos[,12], pos[,5])

  negR <- data.table(rbind(pos[!filterCDS, ], neg[filterThree, ])) # add bad cds to neg
  pos <- data.table(rbind(neg[!filterThree, ], pos[filterCDS, ]))  #!!! UPDATE if new features
  neg <- negR

  predicate <- data.table(rbind(pos, neg))
  colnames(predicate) <- c("coverage", "disengagementScores", "entropyRFP", "fiveRegion","fiveRegionRelative",
                           "floss", "ioScore", "ORFScores","RFPFpkm", "RRS", "RSS", "startCodonCoverage")

  predicate <- fixNAandINFTable(predicate)
  y <- as.factor(c(rep(1, nrow(pos)), rep(0, nrow(neg))))
  predicate <- data.table(y, predicate)

  dCDSThree <- getBestIsoformStartCodonCoverage(cdsAndThree = T)
  pos <- dCDSThree[[1]]
  neg <- dCDSThree[[2]]
  negR <- data.table(rbind(pos[!filterCDS, ], neg[filterThree, ]))      # add bad cds to neg
  pos <- data.table(rbind(neg[!filterThree, ], pos[filterCDS, ]))  #!!! UPDATE if new features
  neg <- negR

  dCDSThree <- data.table(rbind(pos, neg))
  dCDSThree <- dCDSThree[readHits >= quantile(dCDSThree$readHits, 0.849),]
  dInts <- dCDSThree[dCDSThree[, .I[readHits == max(readHits)], by=group]$V1]
  ints <- c(length(filterCDS) + which(!filterThree),
            which(filterCDS),
            which(!filterCDS),
            length(filterCDS) + which(filterThree))
  if (length(ints) != nrow(predicate)) stop("wrong making of ints!")
  predicate$startCodonPerGroupBest <- ints %in% dInts$index
  #analysis
  print(cor(data.matrix(predicate), use = "complete.obs"))

  save(predicate, file = paste0("forests/predicateTables/table_cds3utr_",tissue,".rdata"))
  return(predicate)
}

#' Predict table for uORFs
#'
#' Get ribo-seq features for uORFs
makeUORFPredicateTable <- function(tissue = "all") {
  if(file.exists(paste0("forests/predicateTables/table_",tissue,".rdata"))) {
    load(paste0("forests/predicateTables/table_",tissue,".rdata"))
    return(predictors)
  }
  # group all features by tissue
  # pick 1 tissue
  # make the predicate table

  if(file.exists(paste0("forests/uorfTablePre_", tissue,".rdata"))){
    load(paste0("forests/uorfTablePre_", tissue,".rdata"))
  } else {
    featureNames <- c("coverage", "disengagementScores", "entropyRFP", "fiveRegion","fiveRegionRelative",
                      "floss", "ioScore", "ORFScores","Ribofpkm", "RRS", "RSS", "startCodonCoverage")
    pos <- foreach(i = 1:length(featureNames), .combine = 'cbind', .noexport = "uorfDB", .export = "featureNames") %dopar% {
      setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path
      source("./DataBaseSetup.R")

      getTissueFromFeatureTable(tableName = featureNames[i], tissue = tissue)
    }
    save(pos, file = paste0("forests/uorfTablePre_", tissue,".rdata"))
  }

  predicate <- data.table(pos)
  colnames(predicate) <- c("coverage", "disengagementScores", "entropyRFP", "fiveRegion","fiveRegionRelative",
                           "floss", "ioScore", "ORFScores","RFPFpkm", "RRS", "RSS", "startCodonCoverage")
  predictors <- fixNAandINFTable(predicate)
  # Here is lines with filter
  # filter on isoforms
  d <- getBestIsoformStartCodonCoverage()
  # combine filter with ribo-seq prediction
  d <- d[readHits > quantile(d$readHits, 0.973),]
  d <- d[d[, .I[readHits == max(readHits)], by=group]$V1]
  predictors$startCodonPerGroupBest <- seq.int(1, nrow(predictors)) %in% d$index
  print(cor(data.matrix(predictors), use = "complete.obs"))
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
    goTerms <- readTable("goTerms")
    gc <- readTable("gcContent")
    # make Predicates:
    uorfData <- data.table(rankInTx, lengths, kozak, fractionLengths,
                           StartCodons = as.factor(StartCodons$startCodon),
                           StopCodons = as.factor(StopCodons$stopCodon),
                           distORFCDS, distORFTSS, inFrameCDS, isOverlappingCds,
                           exonExonJunctions, gc, goTerms = as.factor(goTerms))
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
#' get start for each in group, count overlaps, return orf with
#' highest per group
getBestIsoformStartCodonCoverage <- function(cdsAndThree = F) {
  # reduce isoform overlaps by highest start codon reads per group
  if (cdsAndThree) {
    if(file.exists(paste0("forests/bestStartCodonsCDSTHREE.rdata"))) {
      load(paste0("forests/bestStartCodonsCDSTHREE.rdata"))
      return(dCDSThree)
    }
    if(is.null(filterCDS) | is.null(filterThree)) stop("filterCDS or filterThree is null!")

    g <- getCDS(assignIt = F)
    sg <- stopCodons(g, is.sorted = T)
    uo <- uniqueOrder(sg) # <- grouping
    counts <- rowMeans(readTable("cdsstartCodonCoverage"))
    dCDS <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))

    getThreeUTRs()
    threeUTRs <- threeUTRs[widthPerGroup(threeUTRs) > 5]
    sg <- stopCodons(threeUTRs, is.sorted = T)
    uo <- uniqueOrder(sg) # <- grouping
    counts <- rowMeans(readTable("threestartCodonCoverage")) # already at correct length
    dThree <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))
    dThree$group <- dThree$group + max(dCDS$group)
    dThree$index <- dThree$index + max(dCDS$index)
    dCDSThree <- list(dCDS, dThree)
    save(dCDSThree, file = paste0("forests/bestStartCodonsCDSTHREE.rdata"))
    return(dCDSThree)
  }
  if(file.exists(paste0("forests/bestStartCodons.rdata"))) {
    load(paste0("forests/bestStartCodons.rdata"))
    return(d)
  }
  uo <- readTable("stopCodonGrouping")$stopCodonGrouping

  counts <- rowMeans(readTable("startCodonCoverage", with.IDs = F))
  d <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))
  save(d, file = paste0("forests/bestStartCodons.rdata"))
  return(d)
}

#' Train h2o rf model.
#' negDT if you want own samples for that
forest <- function(dt, cv = 10, ntrees = 64, nthreads = 40,  max_mem_size = "200G"){
  library(h2o)
  h2o.init(nthreads = nthreads, max_mem_size = max_mem_size, port = 20050)
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



#' Combine classifier and CAGE data, for final prediction table
#'
makeCombinedPrediction <- function(tissue) {
  # load data
  load(paste0("forests/finalPrediction_filtered",tissue, ".rdata"))
  load(paste0(dataBaseFolder,"/tissueAtlas.rdata"))

  some <- prediction$predict == 1 # set value

  cageTissuesPrediction <- copy(tissueAtlas)
  for(i in colnames(cageTissuesPrediction)[-1]) {
    cageTissuesPrediction[, paste(i) := (tissueAtlas[,i, with=F] & some)]
  }
  insertTable(cageTissuesPrediction, "tissueAtlasByCageAndPred", rmOld = T)

  # sums <- colSums(cageTissuesPrediction)
  # cageTissuesPrediction[, names(which(sums == 0)) := NULL]

  finalCagePred <- rowSums(cageTissuesPrediction[,-1]) > 0
  insertTable(finalCagePred, "finalCAGEuORFPrediction", rmOld = T)

  startCodonMetrics(finalCagePred)
}

validByRibio <- function(coverage, fpkm, startCodonCoverage, fiveRegionRelative) {
  filter <- (coverage > min(quantile(coverage, 0.25), 10)) & (fpkm > quantile(fpkm, 0.15)) &
    (startCodonCoverage > quantile(startCodonCoverage, 0.75)) & fiveRegionRelative > 0.95
  return(filter)
}

overCDS <- function(){
  if(file.exists(paste0("forests/uORFsOverlappingCDS.rdata"))) {
    load(paste0("forests/uORFsOverlappingCDS.rdata"))
    return(isOverCDS)
  }
  grl <- getUorfsInDb()
  getCDS()
  isOverCDS <- countOverlaps(grl, cds) > 0
  save(isOverCDS, file = paste0("forests/uORFsOverlappingCDS.rdata"))
  return(isOverCDS)
}

# Get only specific tissues, or all.
# Grouped by rowMeans
getTissueFromFeatureTable <- function(tableName, tissue) {
  if (tableNotExists(tableName)) stop(paste("table does not exist:",
                                            tableName))
  rpfSamples <- getRiboRNAInfoTable()
  #1. we have tissues in link
  uniqueTissues <- as.character(unique(rpfSamples$tissue))
  riboTable <- readTable(tableName, with.IDs = F)

  if (is.null(tissue) | (tissue == "all")){
    print("Grouping all together")
  } else if ((tissue %in% uniqueTissues)){
    indices <- rpfSamples$tissue == tissue
    riboTable <- riboTable[,indices, with = F]
  } else stop("tissue does not exist in db")

  return(rowMeans(riboTable))
}

filterHardOnes <- function(prediction, tissue = "all"){
  prediction$filtered <- rep(F, nrow(prediction))
  uorfTable <- makeUORFPredicateTable()
  uorfData <- getAllSequenceFeaturesTable()
  grl <- getUorfsInDb()
  getCDS()
  getCageTx()
  getFasta()
  table <- startCodonMetrics(as.logical(prediction[,3] >= 0.50))
  badStarts <- table[,1] > 0 & table[,2] < 0
  badStartCodons <- rownames(table)[badStarts]
  goodStartCodons <- paste(rownames(table)[!badStarts], collapse = "|")
  goodUORFs <- grl[uorfData$StartCodons %in%  rownames(table)[!badStarts]]
  res = c()
  for(codon in badStartCodons) {
    # find region
    agg <- uorfData$StartCodons == codon & prediction[,3] >= 0.75
    starts <- startSites(grl[agg], T, T, T)
    # make string
    hits <- grep(x = ORFik:::startRegionString(grl[agg], tx, fa, 6, 9), pattern = goodStartCodons)
    hitsUp <- grep(x = ORFik:::startRegionString(grl[agg], tx, fa, 5, 0), pattern = stopDefinition(1))
    hitsOverlapsBetter <- starts %over% goodUORFs
    hitsCDS <- to(findOverlaps(startSites(cds, T, T, T), starts, maxgap = 3))
    valid <- (!(seq.int(1,sum(agg)) %in% c(hits, hitsUp, hitsOverlapsBetter, hitsCDS)))

    # find indices
    notBestStart <- (uorfTable$startCodonPerGroupBest == F)[agg]
    index <- which((!notBestStart & valid))
    toKeep <- which(agg)[index]
    res <- c(res, toKeep)
  }
  if(length(res) != length(unique(res))) stop("error in res creation!")
  hits <- (grep(x = uorfData$StartCodons, pattern = paste(badStartCodons, collapse = "|")))

  prediction$predict[hits] <- 0
  prediction$p0[hits] <- 1
  prediction$p1[hits] <- 0
  prediction$filtered[hits] <- T

  prediction$predict[res] <- 1
  prediction$p0[res] <- 0
  prediction$p1[res] <- 1
  startCodonMetrics(as.logical(prediction[,3] >= 0.50))
  save(prediction, file = paste0("forests/finalPrediction_filtered",tissue, ".rdata"))
  return(prediction)
}
