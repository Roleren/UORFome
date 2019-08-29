libraryTypes <- function(df){
  if (is(df, "data.frame")) {
    return(colnames(df)[!(colnames(df) %in% c("stage", "type", "rep", "varName"))])
  } else if (is(df, "character") | is(df, "factor")) {
    return(gsub("_.*", x = df, replacement = ""))
  } else stop("library types must be data.frame or character vector!")
}

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
makeSummarizedExperimentFromBam <- function(df, txdb, saveName = NULL, longestPerGene = TRUE, 
                                            geneOrTxNames = "gene") {
  library(SummarizedExperiment)
  libTypes <- libraryTypes(df)
  validateExperiments(df)
  
  txdb <- loadTxdb(txdb)
  tx <- loadRegion(txdb, "mrna")
  if (geneOrTxNames == "gene") names(tx) <- ORFik:::txNamesToGeneNames(names(tx), txdb)
  
  varNames <- bamVarName(df)
  outputBams(df, tx)
  
  rawCounts <- data.table(matrix(0, ncol = length(varNames), nrow = length(tx)))
  for (i in 1:length(varNames)) { # For each stage
    print(varNames[i])
    co <- countOverlaps(tx, get(varNames[i]))
    rawCounts[, (paste0("V",i)) := co]
  }
  mat <- as.matrix(rawCounts);colnames(mat) <- NULL
  
  colData <- DataFrame(SAMPLE = bamVarName(df, TRUE),
                       row.names=varNames)
  if (!is.null(df$rep)) colData$replicate <- df$rep
  
  res <- SummarizedExperiment(assays=list(counts=mat), rowRanges=tx, colData=colData)
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
  return(dds)
}

validateExperiments <- function(df) {
  libTypes <- libraryTypes(df)
  if (!is(df, "data.frame")) stop("df must be data.frame!")
  if (!all((c("stage", "type") %in% colnames(df)))) stop("stage and type must be colnames in df!")
  if (length(libTypes) == 0) stop("df have no valid sequencing libraries!")
  if (nrow(df) == 0) stop("df must have at least 1 row!")
  
  emptyFiles <- c()
  for (i in libTypes) {
    emptyFiles <- c(emptyFiles, as.numeric(sapply(as.character(df[, i]), file.size)) == 0)
  }
  if (any(emptyFiles)) {
    print(cbind(df[emptyFiles, libTypes], which(emptyFiles)))
    stop("Empty files in list, see above for which")
  }
}

#' Get difference between wild type and target
#' Per library type like RNA-seq or Ribo-seq.
#' @param dt a data.table of counts or fpkm etc.
#' @return a data.table with differences
SEdif <- function(dt) {
  dif <- data.table()
  bamVars <- colnames(dt)
  for(i in 1:nrow(df)) {
    a <- dt[, get(bamVars[(i+6)]) - get(bamVars[(i)])]
    dif <- cbind(dif, a)
  }
  colnames(dif) <- paste(libraryTypes(df), sort(df$stage), sep = "_")
  return(dif)
}

SESplit <- function(dif, splitCol, score = 15) {
  dif$ismiRNA <- splitCol
  vdif <- dif[dt$RNA_WT_2 >= score & dt$RFP_WT_2 >= score & dt$RNA_MZ_2 >= score & dt$RFP_MZ_2 >= score, ]
  ismiRNA <- vdif$ismiRNA
  vdif$ismiRNA <- NULL
  vdifm <- melt(vdif)
  vdifm$ismiRNA <- rep(ismiRNA, 6)
  vdifm$stage <- gsub(".*_", x = vdifm$variable, replacement = "")
  vdifm$type <- libraryTypes(vdifm$variable)
  d <- vdifm[type == "RNA",]
  d <- cbind(d, vdifm$value[vdifm$type == "RFP"])
  colnames(d) <- c("variable", "RNA", "ismiRNA", "stage", "type", "RFP")
  return(d)
}
