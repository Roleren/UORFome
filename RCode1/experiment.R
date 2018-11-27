# experiment validation
# hek vs hela
findGoodUORFs <- function(tissue = "all"){
  grl <- getUorfsInDb()
  RNA <- readTable("RNAByTissueMean")
  rowMeans <- rowMeans(RNA[,2:ncol(RNA)])
  rowMeans <- rowMeans[data.table::chmatch(txNames(grl), RNA$txNames)]
  cageTissuesPrediction <- readTable("tissueAtlasByCageAndPred")
  load(file = paste0("forests/prediction_", tissue, ".rdata"))
  hits <- (cageTissuesPrediction$kidney != cageTissuesPrediction$ovary) & prediction$predict == 1
  startCodonMetrics(hits)
  hits <- which(hits)
  best <- hits[which.max(rowMeans[hits])]
  bestOrder <- hits[order(rowMeans[hits], decreasing = T)]
  uorfData$StartCodons[bestOrder[1:10]]
  
  
  uorfTable <- makeUORFPredicateTable(tissue)
  uorfTable[best,]
  
  summary(uorfTable$RFPFpkm)
  summary(uorfTable$ORFScores)
  summary(uorfTable$startCodonCoverage)
  # conclusion, very good scores of rna seq and ribo seq
  # conserved in many species
  # Ribosome 
  
  # check some genes here:
  refTable <- readTable("uORFTxToGene") 
  #ATF4 confirmed
  ORFsInGeneMetrics("ATF4")
  #ATF5 Could not find
  ORFsInGeneMetrics("ATF5")
  #ABCC2 confirmed two
  ORFsInGeneMetrics("ABCC2")
  #ADH5 found none
  ORFsInGeneMetrics("ADH5")
  #ADAM10
  ORFsInGeneMetrics("ADAM10")
  #PRKCH
  ORFsInGeneMetrics("PRKCH")
  #KLF9
  ORFsInGeneMetrics("KLF9")
  #MC2R
  ORFsInGeneMetrics("MC2R")
  #BIRC3
  ORFsInGeneMetrics("BIRC3")
}

#' Get indices of uORFs from genes
#' @param hgncSymbol a character vector of gene symbols
#' @return a integer vector of indices
getORFsGeneSymbols <- function(hgncSymbol = "ATF4", refTable = refTable){
  library(biomaRt)
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  geneHits <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol'),
                      filters = 'hgnc_symbol', values = hgncSymbol, mart = ensembl)
  if(nrow(geneHits) == 0) stop(p("could not find any genes with the name: ", hgncSymbol))
  
  return(which(refTable$geneNames == geneHits[1, 1]))
}

getORFsGoTerms <- function(uORFGenes){
  library(biomartr)
  old <- uORFGenes
  uORFGenes <- unique(uORFGenes)
  Go <- biomartr::getGO(organism = "Homo sapiens", 
                            genes    = uORFGenes,
                            filters  = "ensembl_gene_id")
  desc <- Go$goslim_goa_description
  
  return(desc[data.table::chmatch(old, uORFGenes)])
}

ORFsInGeneMetrics <- function(hgncSymbol = "ATF4"){
  hits <- getORFsGeneSymbols(hgncSymbol = hgncSymbol, refTable = refTable)
  print(paste("found ", sum(finalCagePred[hits]), "/", length(hits),
              "predicted translated uorfs in", hgncSymbol))
  if(sum(finalCagePred[hits]) == 0) return(NULL)
  hits <- hits[finalCagePred[hits]]
  print(uorfData$distORFCDS[hits])
  print(as.character(uorfData$StartCodons[hits]))
  print(uorfData$lengths[hits])
  print(uorfData$rankInTx[hits])
  return(hits)
}