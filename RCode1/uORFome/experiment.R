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

getAllORFGeneSymbols <- function(geneNames){
  library(biomaRt)
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  uniqueGenes <- unique(geneNames)
  geneHits <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                    filters = 'ensembl_gene_id', values = uniqueGenes, mart = ensembl)
  group2 <- data.table::chmatch(geneNames, geneHits$ensembl_gene_id)
  
  return(data.table(geneNames = geneNames, symbol = geneHits$hgnc_symbol[group2]))
}

geneSymbolsTo <- function(geneNames, org.db = org.Dr.eg.db){
  library(clusterProfiler)
  
  geneNames <- as.character(geneNames)
  
  uniqueGenes <- unique(geneNames)
  geneHits <- bitr(uniqueGenes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Dr.eg.db)
  group2 <- data.table::chmatch(geneNames, geneHits$ENSEMBL)
  
  return(data.table(symbol = geneHits$SYMBOL[group2]))
}

getORFsGoTerms <- function(uORFGenes, organism = "Homo sapiens"){
  library(biomartr)
  old <- uORFGenes
  uORFGenes <- unique(uORFGenes)
  Go <- biomartr::getGO(organism = organism, 
                            genes    = uORFGenes,
                            filters  = "ensembl_gene_id")
  desc <- Go$goslim_goa_description
  
  return(desc[data.table::chmatch(as.character(old), as.character(uORFGenes))])
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

primersToMake <- function(){
  uorfTable <- makeUORFPredicateTable()
  uorfData <- getAllSequenceFeaturesTable()
  teTable <- readTable("TEByTissueMeanWithoutInf")
  grl <- getUorfsInDb()
  hits <- uorfData$StartCodons == "ATG" & uorfTable$RFPFpkm > 0.5 & uorfTable$startCodonCoverage > 1.0 & uorfTable$ORFScores > 0.01 & teTable$HeLa > teTable$HEK293*5
  grl <- grl[hits]
  # now check them in IGV
  #candidates are here
  gr <- GRanges(c("chr1:6026375-6026398:+", "chr1:10210717-10210740:+",
            "chr1:10430850-10431044:+","chr1:11273233-11273253:+","chr1:11980456-11980482:+",
            "chr1:153536104-153536118:-","chr1:153536171-153536233:-", "chr2:30148292-30148309:+",
            "chr3:55480956-55480997:-", "chr5:32711306-32711335:+"))
  names <- c( "ENST00000341524(HeLa, not Hek)","ENST00000377093(Both)",
              "ENST00000602296(both)","ENST00000376810(both)","ENST00000444836(both)",
              "ENST00000496817(HeLa, not HEK)1","ENST00000496817(HeLa, not HEK)",
              "ENST00000379519(Both)","ENST00000614415(HeLa, not HEK)","ENST00000265074(HeLa, not HEK)"
  )
  names(gr) <- names    
  
  control <- GRanges(c("chr22:39520648-39520659:+"))
  names(control) <- c("ATFuORF2")
  
  gr <- c(gr, control)
  getSequencesFromFasta(groupGRangesBy(gr))
  
  fileConn<-file("output.txt")
  writeLines(c("Hello","World"), fileConn)
  close(fileConn)
  write.(seqs, "primers.fa")
}

uORFAnalysis <- function(){
  uorfTable <- makeUORFPredicateTable()
  uorfData <- getAllSequenceFeaturesTable()
  teTable <- readTable("TEByTissueMeanWithoutInf")
  grl <- getUorfsInDb()
  RNA <- readTable("RNAByTissueMean")
  rowMeans <- rowMeans(RNA[,2:ncol(RNA)])
  rowMeans <- rowMeans[data.table::chmatch(txNames(grl), RNA$txNames)]
  
  hits <- which(final$Matrix == 1 & rowMeans > 10 & uorfTable$RFPFpkm > 0.5 &
                  uorfTable$startCodonCoverage > 1.0 & uorfData$eejuORF == 1 & uorfData$isOverlappingCds == F)
  
  for(i in hits[25]) {
    print(p("uORF: ", names(grl[i])))
    print(uorfData[i,])
    print(uorfTable[i,])
    print(grl[i])
  }
  notOver <- !overCDS()
  hits <- rowMeans > 10 & uorfTable$RFPFpkm > 0.5 &
                  uorfTable$startCodonCoverage > 1.0 & uorfData$eejuORF == 1 & uorfData$isOverlappingCds == F
  
  hits <- which((prediction$p1 > 0.9 & hits[notOver]))
  grl2 <- grl[notOver]
  uorfData2 <- uorfData[notOver,]
  uorfTable2 <- uorfTable[notOver,]
  
  for(i in hits[34]) {
    print(p("uORF: ", names(grl2[i])))
    print(uorfData2[i,])
    print(uorfTable2[i,])
    print(grl2[i])
  }
  
  sum(uorfData$StartCodons == "ATG")
  sum(uorfData2$StartCodons == "ATG" & prediction$predict == 1)
  
  grl2[which(uorfData2$StartCodons == "ATG" & prediction$predict == 0 & rowMeans[notOver] > 10)][1]
  cages <- readTable("tissueAtlasByCage")
  cages$uorfID <- NULL
  cage <- cages[notOver,]
  hits <- which(cage$kidney == 0 & cage$ovary == 0 & cage$brain == 0 & cage$blood == 0 & cage$prostate == 0 & prediction$predict == 1 & rowMeans[notOver] > 10)
  i = 200
  uorfData2$StartCodons[hits][i]
  grl2[hits][i]
  
  # hit matrix oo, on, no, oo
  sum(cage$kidney == 1 & cage$ovary == 1 & cage$brain == 1 & cage$blood == 1 & cage$prostate == 1 & prediction$predict == 1)
  sum(cage$kidney == 1 & cage$ovary == 1 & cage$brain == 1 & cage$blood == 1 & cage$prostate == 1 & prediction$predict == 0)
  sum(cage$kidney == 0 & cage$ovary == 0 & cage$brain == 0 & cage$blood == 0 & cage$prostate == 0 & prediction$predict == 1)
  sum(cage$kidney == 0 & cage$ovary == 0 & cage$brain == 0 & cage$blood == 0 & cage$prostate == 0 & prediction$predict == 0)
  
  
}

#' ATF4 gene, transcript isoform ENST00000404241, 6 should be at start chr22 39,520,586
#' Conclusion: It now finds all three 
ATF4Check <- function(){
    
  order(prediction$p1, decreasing = T)
  
  hits6 <- which((ORFik:::txNames(grl) == "ENST00000404241") & widthPerGroup(grl, F) == 6)
  grl[hits6] # nr 2 is correct
  hitsPred <- hits6 %in% which(prediction$predict == 1) 
  
  hits12 <- which((ORFik:::txNames(grl) == "ENST00000404241") & widthPerGroup(grl, F) == 12)
  grl[hits12] # nr 3 is correct
  hitsPred <- hits12 %in% which(prediction$predict == 1) 
  
  hitsLast <- which((ORFik:::txNames(grl) == "ENST00000404241") & startSites(grl, F, F, T) == 39520747)
  grl[hitsLast]
  hitsPred <- hitsLast %in% which(prediction$predict == 1) 
}
