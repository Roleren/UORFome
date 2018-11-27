# Gviz example
library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
library(gridExtra)
library(foreach)

# bioconductor
library(Gviz)
library(Biostrings)
library(Rsamtools)

browserPlots <- function(){
  refTable <- readTable("uORFTxToGene") 
  if(!exists("geneTrack")) geneTrack <- getGeneTrack()
  grl <- getUorfsInDb()
  # ribo seq
  riboName <- list.files(rfpFolder)[103]
  rfp <- fread.bed(p(rfpFolder, riboName))
  riboName <- gsub(x = gsub(x = riboName, pattern = ".reads.*", replacement = ""), "\\.", " ")
  uorfTable <- makeUORFPredicateTable()
  uorfData <- getAllSequenceFeaturesTable()
  
  # ATF4 uORFs
  gene <- "ATF4"
  hits <- getORFsGeneSymbols(gene, refTable)
  # verified uORFs
  verifiedUORFs <- GRanges('chr22', IRanges(c(39520586,39520648, 39520747, 39521354), 
                                 width = c(6,12,5,175)), '+', names = c('1','2','3','3'))
  verifiedUORFs <- groupGRangesBy(verifiedUORFs, verifiedUORFs$names)
  # normal
  transluORFs <- grl[hits]
  predicted <- verifiedIndices(transluORFs, verifiedUORFs)
  browseruORFRegion(transluORFs, rfp, geneTrack,riboName, gene, predicted)
  browseruORFRegion(verifiedUORFs, rfp, geneTrack,riboName, gene)
  browseruORFRegion(transluORFs[pr], rfp, geneTrack,riboName, gene)
  #ABCC2 crap found only
  gene <- "ABCC2"
  verifiedUORFs <- GRanges("chr10",IRanges(99782699, width = 48), "+")
  verifiedUORFs <- groupGRangesBy(verifiedUORFs, "1")
  hits <- getORFsGeneSymbols(gene, refTable)
  transluORFs <- grl[hits]
  predicted <- verifiedIndices(transluORFs, verifiedUORFs)
  browseruORFRegion(transluORFs, rfp, geneTrack,riboName, gene, predicted)
  # BIRC3
  gene <- "BIRC3"
  hits <- getORFsGeneSymbols(gene, refTable)
  transluORFs <- grl[hits]
  browseruORFRegion(transluORFs, rfp, geneTrack,riboName, gene)
  #ADAM10
  gene <- "ADAM10"
  hits <- getORFsGeneSymbols(gene, refTable)
  transluORFs <- grl[hits]
  #HSPA8
  gene <- "HSPA8"
  hits <- getORFsGeneSymbols(gene, refTable)
  transluORFs <- grl[hits]
  #TP63
  gene <- "TP63"
  hits <- getORFsGeneSymbols(gene, refTable)
  transluORFs <- grl[hits]
  #GABRA1
  gene <- "GABRA1"
  hits <- getORFsGeneSymbols(gene, refTable)
  transluORFs <- grl[hits]
  
  # for each rfp if needed
  foreach(i=list.files(rfpFolder), .inorder = F, .packages = c("ORFik", "data.table", "Gviz")) %dopar% {
    rfp <- fread.bed(p(rfpFolder, i))
    riboName <- gsub(x = gsub(x = i, pattern = ".reads.*", replacement = ""), "\\.", " ")
    browseruORFRegion(transluORFs, rfp, geneTrack,riboName, gene)
  }
}
#' Track plot of uORF region, with ribo seq
#' @param transluORFs the uORFs as GRanges or GRangesList
#' @param rfp character vector or GRanges of ribo seq
#' @param geneTrack track from txdb
#' @param name output name of pdf
#' @param predicted (darkred), a vector as factor for color, or a single valid color. 
#' @param ht highlight region
#' @return NULL, will save a pdf
browseruORFRegion <- function(transluORFs, rfp, geneTrack = NULL,
                              riboName = "unknown sample \nP-site aligned", gene = "unknown", predicted = "darkred",
                              ht = NULL ){
  # define region of interest:
  if(ORFik:::is.grl(transluORFs)) transluORFs <- unlistGrl(transluORFs)
  myregion_extend <- uorfRegion(transluORFs)
  mychrom <- as.character(seqnames(myregion_extend))
  exonMatching <- chmatch(transluORFs$names, unique(transluORFs$names))
  orfNames <- paste0("uORF_", exonMatching)
  if(length(predicted) > 1) predicted <- predicted[exonMatching]
  
  options(ucscChromosomeNames=FALSE)
  genomeTrack<- GenomeAxisTrack()
  
  transluORFTrack<- GeneRegionTrack(transluORFs,
                                    showId=TRUE,
                                    id=orfNames,
                                    name=p('uORF predictions ', gene),
                                    transcript=orfNames,
                                    symbol=orfNames
                                    )
  if(is.null(geneTrack)) geneTrack <- getGeneTrack()
  
  # read in ribo-seq wig files here, either shifted all, or just a few
  if(is.character(rfp)) rfp <- fread.bed(rfp)

  scores <- unlist(NumericList(ORFik:::coveragePerTiling(GRangesList(region = myregion_extend), rfp,
                                                         is.sorted = T, keep.names = F)), use.names = FALSE)
  #scores <- unlist(NumericList(ORFik:::coveragePerTiling(GRangesList(region = myregion_extend), GRanges(mychrom, rep(end(myregion_extend) - 1000, 10), "-"), is.sorted = T, keep.names = F)), use.names = FALSE)
  if(!strandBool(myregion_extend)) scores <- rev(scores)
  riboTrack<- DataTrack(chromosome = mychrom,
                         start=start(myregion_extend):end(myregion_extend),
                         width=0,
                         data=scores,
                         genome='',
                         name=p("ribo-seq\n", riboName))
  
  displayPars(riboTrack)<- list(showAxis=TRUE,
                                 type='h',
                                 col='darkblue',
                                 fontcolor.title='black',
                                 col.axis='black',
                                 cex.axis=0.5,
                                 background.title='white')
  displayPars(geneTrack)<- list(
    fontcolor.title='black',
    col.axis='black',
    cex.axis=0.5,
    background.title='white')
  
  displayPars(transluORFTrack)<- list(
    fontcolor.title='black',
    col.axis='black',
    fill=predicted,
    cex.group=0.25,
    cex.axis=0.5,
    background.title='white')
  
  mytracks<- list(genomeTrack,
                  riboTrack,
                  transluORFTrack,
                  geneTrack
  )
  #mytracks<- list(genomeTrack, riboTrack)
  showBrowserTrack(mytracks, myregion_extend, transluORFs, mychrom, gene)
  return(NULL)
}


getGeneTrack <- function(){
  if(file.exists(paste0(dataBaseFolder, "/browserTracks/", "geneTracks.rdata"))){
    load(paste0(dataBaseFolder, "/browserTracks/", "geneTracks.rdata"))
    return(geneTrack)
  }
  options(ucscChromosomeNames=FALSE)
  getGTF()
  geneTrack <- GeneRegionTrack(Gtf,name="gene models",showId=TRUE,geneSymbol=TRUE)
  # can get a txdb object
  save(geneTrack, file = paste0(dataBaseFolder, "/browserTracks/", "geneTracks.rdata"))
  return(geneTrack)
}

verifiedIndices <- function(transluORFs, verifiedUORFs){
  predicted <- factor(rep(1, length(transluORFs)), levels = c(1,2))
  pr <- which(startSites(transluORFs) %in% startSites(verifiedUORFs) &
                stopSites(transluORFs) %in% stopSites(verifiedUORFs))
  predicted[pr] <- as.factor(rep(2, length(pr)))
  return(predicted)
}

uorfRegion <- function(transluORFs, downstream = 1000, upstream = 1000){
  mychrom= as.character(seqnames(transluORFs[1]))
  mystart=min(start(transluORFs))
  myend=max(end(transluORFs))
  myregion_extend=GRanges(mychrom, IRanges(mystart-downstream, myend+upstream))
  strand(myregion_extend) = as.character(strand(transluORFs[1]))
  return(myregion_extend)
}

showBrowserTrack <- function(mytracks, myregion_extend, transluORFs, mychrom, gene, 
                             name = paste0(dataBaseFolder, "/browserTracks/", gene ,
                                           "/browserTrack ",riboName,".pdf")){
  if(!dir.exists(paste0(dataBaseFolder, "/browserTracks/", gene))){
    dir.create(paste0(dataBaseFolder, "/browserTracks/", gene))
  }
  
  if(length(transluORFs) > 250) {
    width = 14
    height = 10
  } else {
    width = 12
    height = 8
  }
  pdf(name, width=width, height=height)
  plotTracks(mytracks,
             from = start(myregion_extend) - 1300,
             to = end(myregion_extend) + 500,
             chromosome = mychrom)
  dev.off()
  return(NULL)
}
# make a cage track # find the right fibroblast cell line
#cageFile=''
# get the sequence in region
# fastaFile='/export/valenfs/data/references/GRCh38_human/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa'
# fa=FaFile(fastaFile)
# refseq<- readDNAStringSet(fastaFile)
# #refseq<- getSeq(fa, myregion_extend)
# names(refseq)<- gsub(" .*","",names(refseq))
# seqTrack<- SequenceTrack(refseq, chromosome = mychrom)
# displayPars(seqTrack)<- list(fontcolor=biovizBase::getBioColor(),noLetters=FALSE,fontsize=6)

# grl <- getUorfsInDb()
# orfs <- grl
# orfTrack<- AnnotationTrack(orfs[orfs %over% myregion_extend], name="uORFs based on sequence", col='darkgrey', fill='darkgrey')
# displayPars(orfTrack)<- list(
#   fontcolor.title='black',
#   col.axis='black',
#   fill='darkgrey',
#   cex.axis=0.5,
#   background.title='white')
# # use this to highlight
# mytracks<- list(genomeTrack,
#                 #cageTrack,
#                 riboTrackF,
#                 #riboTrackR,
#                 transluORFTrack,
#                 #orfTrack,
#                 orfTrack
#                 #geneTrack
# )
# pdf('browser_close.pdf',width=10,height=6)
# plotTracks(mytracks,
#            from = start(transluORFs['uORF_2'])-85,
#            to = end(transluORFs['uORF_2'])+10,
#            chromosome = seqnames(transluORFs['uORF_2']))
# 
# dev.off()
# ht <- HighlightTrack(trackList = mytracks,  start = start(transluORFs['uORF_2'])-80,
#                      end = end(transluORFs['uORF_2'])+20,
#                      chromosome = seqnames(transluORFs['uORF_2']))