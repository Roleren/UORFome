# experiment validation
# hek vs hela
findGoodUORFs <- function(tissue = "all"){
  grl <- getUorfsInDb(T,T,T)
  RNA <- readTable("RNAByTissueMean")
  rowMeans <- rowMeans(RNA[,2:ncol(RNA)])
  rowMeans <- rowMeans[data.table::chmatch(txNames(grl), RNA$txNames)]
  cageTissuesPrediction <- readTable("tissueAtlasByCageAndPred")
  load(file = paste0("forests/prediction_", tissue, ".rdata"))
  hits <- which((cageTissuesPrediction$kidney != cageTissuesPrediction$ovary) & prediction$predict == 1)
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
  dt <- readTable("uORFTxToGene") 
  #ATF4
  atf4Hits <- getUORFsGeneSymbols(refTable = dt)
  finalCagePred[atf4Hits]
  
}

#' Get indices of uORFs from genes
#' @param hgncSymbol a character vector of gene symbols
#' @return a integer vector of indices
getORFsGeneSymbols <- function(hgncSymbol = "ATF4", refTable){
  library(biomaRt)
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  geneHits <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol'),
                      filters = 'hgnc_symbol', values = hgncSymbol, mart = ensembl)
  if(nrow(geneHits) == 0) stop(p("could not find any genes with name", hgncSymbol))
  
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