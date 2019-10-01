#' Make a summerizedExperiment object from bam files
#' 
#' If txdb or gtf path is added, it is a rangedSummerizedExperiment
#' @param df a data.frame with minimum 3 columns (2 are static).
#' 1. stage, stage of sample
#' 2. rep, replicate of sample
#' 3. Any name, paths to files like RNA = c(path/to/file), that are used are column in assay
#' @inheritParams loadTxdb
#' @param saveName a character (default NULL), if set save experiment to path given. Always saved as .rds., 
#' it is optional to add .rds, it will be added for you if not present.
#' @param longestPerGene a logical (default TRUE), if FALSE all transcript isoforms per gene.
#' @param geneOrTxNames a character (default gene), if tx use with tx names
#' @return a summerizedExperiment object
makeSummarizedExperimentFromBam <- function(df, saveName = NULL, longestPerGene = TRUE, 
                                            geneOrTxNames = "gene", region = "mrna", type = "count") {
  library(SummarizedExperiment)
  if(file.exists(saveName)) {
    return(readRDS(saveName))
  }
  libTypes <- libraryTypes(df)
  validateExperiments(df)
  
  txdb <- loadTxdb(df@txdb)
  tx <- loadRegion(txdb, region)
  if (geneOrTxNames == "gene") names(tx) <- ORFik:::txNamesToGeneNames(names(tx), txdb)
  
  varNames <- bamVarName(df)
  outputLibs(df, tx)
  
  rawCounts <- data.table(matrix(0, ncol = length(varNames), nrow = length(tx)))
  for (i in 1:length(varNames)) { # For each sample
    print(varNames[i])
    co <- countOverlaps(tx, get(varNames[i]))
    rawCounts[, (paste0("V",i)) := co]
  }
  mat <- as.matrix(rawCounts);colnames(mat) <- NULL
  
  colData <- DataFrame(SAMPLE = bamVarName(df, TRUE),
                       row.names=varNames)
  if (!is.null(df$rep)) colData$replicate <- df$rep
  
  res <- SummarizedExperiment(assays=list(counts=mat), rowRanges=tx, colData=colData)
  if (type == "fpkm") res <- as.data.table(scoreSummarizedExperiment(res, score = "fpkm"))
  if(!is.null(saveName)) {
    saveRDS(res, file = saveName)
  }
  return(res)
}

scoreSummarizedExperiment <- function(final, score = "transcriptNormalized", collapse = FALSE) {
  library(DESeq2)
  if (is.factor(final$SAMPLE))
    final$SAMPLE <- factor(final$SAMPLE, levels = levels(final$SAMPLE)[levels(final$SAMPLE) %in% unique(colData(final)$SAMPLE)])
  
  if (collapse) {
    collapsedAll <- collapseReplicates(final, final$SAMPLE)
    assay(collapsedAll) <- ceiling(assay(collapsedAll) / t(matrix(as.double(table(colData(final)$SAMPLE)),
                                                                  ncol = nrow(assay(collapsedAll)) ,
                                                                  nrow = length(unique(colData(final)$SAMPLE)))))
  } else collapsedAll <- final
  
  
  dds <- DESeqDataSet(collapsedAll, design = ~ SAMPLE)
  fpkmCollapsed <- DESeq2::fpkm(dds)
  if (score == "transcriptNormalized") {
    normalization <- matrix(rep(rowSums2(fpkmCollapsed), ncol(fpkmCollapsed)), ncol = ncol(fpkmCollapsed))
    fpkmTranscriptNormalized <- fpkmCollapsed / normalization
    assay(dds) <- fpkmTranscriptNormalized
    return(dds)
  } else if (score == "fpkm") {
    return(fpkmCollapsed)
  }
  else if (score == "fpkm") {
    return(fpkmCollapsed)
  } 
  else if (score == "log2fpkm") {
    return(log2(fpkmCollapsed))
  }
  else if (score == "log10fpkm") {
    return(log10(fpkmCollapsed))
  }
  return(dds)
}

#' Get difference between wild type and target
#' Per library type like RNA-seq or Ribo-seq.
#' @param dt a data.table of counts or fpkm etc.
#' @return a data.table with differences
SEdif <- function(dt, df) {
  dif <- data.table()
  bamVars <- colnames(dt)
  
  bamVarsT <- bamVarName(df, skip.type = T)
  dists <- which(!is.na(chmatch(bamVarsT, bamVarsT[1])))
  if (length(dists) != 2) stop("Wrong naming!")
  seperator <- dists[2] - dists[1]
  
  for(i in 1:(ncol(dt)/2)) {
    a <- dt[, log2((get(bamVars[(i+seperator)]) + 0.00001) / (get(bamVars[(i)]) + 0.00001))] 
    dif <- cbind(dif, a)
  }
  colnames(dif) <- bamVarsT[1:(ncol(dt)/2)]
  return(dif)
}

SESplit <- function(dif, splitCol, score = 15, dtS, df) {
  whichMatch <- rowSums(dtS > score) == ncol(dtS)
  
  vdif <- dif[whichMatch, ]
  ismiRNA <- splitCol[whichMatch]
  
  vdifm <- melt(vdif)
  vdifm$ismiRNA <- rep(ismiRNA, ncol(dif))
  vdifm$stage <- gsub(".*_", x = vdifm$variable, replacement = "")
  
  if(!is.null(df$experiment)) {
    vdifm$type <- sub("_.*", sub("..._", x = vdifm$variable,
                                 replacement = "", perl = T), replacement = "")
  } else {
    vdifm$type <- libraryTypes(vdifm$variable)
  }
  
  d <- data.table()
  libTypes <- libraryTypes(df)
  d <- vdifm[type == libTypes[1],]
  
  for(i in libTypes[-1]) {
    d <- cbind(d, vdifm$value[vdifm$type == i])
  }
  
  
  colnames(d) <- c(c("variable", libTypes[1], "ismiRNA", "stage", "type"), libTypes[-1])
  return(d)
}
