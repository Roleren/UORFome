
library(ggplot2)
library(gridExtra)
library(GenomicAlignments)
library(GenomicFeatures)
library(ORFik)
getGTF()
getCDS()
getCageTx()
getThreeUTRs()
oldThree <- threeUTRs
threeUTRs <- threeUTRs[widthPerGroup(threeUTRs, F) > 5]
threeUTRs <- threeUTRs[countOverlaps(threeUTRs, cds) == 0]
rfp <- fread.bed(p(rfpFolder, list.files(rfpFolder)[41]))
footprints <- rfp # assign here
uorfTable <- makeUORFPredicateTable()
uorfData <- getAllSequenceFeaturesTable()

grl <- getUorfsInDb()
grl <- grl[uorfData$StartCodons == "ATG" & uorfTable$RFPFpkm > 0.5 & uorfTable$startCodonCoverage > 0.1]
grl <- grl[widthPerGroup(grl) > 20]
gr <- unlistGrl(grl)
names(gr) <- gr$names
grl <- groupGRangesBy(gr)
df <- ORFik:::riboTISCoverageProportion(cds, tx, footprints, average = TRUE,
                                        onlyProportion = FALSE)

df2 <- ORFik:::riboTISCoverageProportion(oldThree, tx, footprints, average = TRUE,
                                         onlyProportion = FALSE)
df3 <- ORFik:::riboTISCoverageProportion(threeUTRs, tx, footprints, average = TRUE,
                                         onlyProportion = FALSE)
df4 <- riboTISCoverageProportion(grl, tx, footprints, average = TRUE,
                                         onlyProportion = FALSE)
#  heatmaps
pdf("cdsHeatmap.pdf", width=13, height=9)
mLen5<-ggplot(df , aes(x=pos, y=length, fill=prop)) + geom_tile()  +
  scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(-2,2),name="Proportions") +
  xlab("Position relative to start codon") + ylab("Protected fragment length") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 12))
plot(mLen5)
dev.off()
#print(mLen5)
pdf("threeUTRsHeatmapAll.pdf", width=13, height=9)
mLen5<-ggplot(df2 , aes(x=pos, y=length, fill=prop)) + geom_tile()  +
  scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(-2,2),name="Proportions") +
  xlab("Position relative to start codon") + ylab("Protected fragment length") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 12))
plot(mLen5)
dev.off()
pdf("./DataBasePlots/threeUTRsHeatmap.pdf", width=13, height=9)
mLen5<-ggplot(df3 , aes(x=pos, y=length, fill=prop)) + geom_tile()  +
  scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(-2,2),name="Proportions") +
  xlab("Position relative to start codon") + ylab("Protected fragment length") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 12))
plot(mLen5)
dev.off()
pdf("./DataBasePlots/UORFsATGNew.pdf", width=12, height=8)
mLen5<-ggplot(df4, aes(x=pos, y=length, fill=prop)) + geom_tile()  +
  scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(-2,2),name="Proportions") +
  xlab("Position relative to start codon") + ylab("Protected fragment length") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 12))
plot(mLen5)
dev.off()
# find orfs

# faFile <- "home/roler/Desktop/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.fai"
#
# fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE)
# tx_seqs <- extractTranscriptSeqs(faFile, fiveUTRs)
# findMapORFs()
cdsProp <- split(Rle(df$prop), df$length)
names(cdsProp) <- NULL
# positive set
test <- ORFik:::riboTISCoverageProportion(cds, tx, footprints,
                                          average = FALSE)
# negative set
test2 <- ORFik:::riboTISCoverageProportion(threeUTRs, tx, footprints,
                                           average = FALSE)

dif <- lapply(seq.int(1:length(test)), function(x) abs(test[[x]] - cdsProp[x]))
dif2 <- lapply(seq.int(1:length(test)), function(x) sum(dif[[x]]))

dif3 <-  lapply(seq.int(1:length(test2)), function(x) abs(test2[[x]] - cdsProp[x]))
dif4 <- lapply(seq.int(1:length(test2)), function(x) sum(dif3[[x]]))


a <- initiationScore(cds, cds, tx, rfp)

windowPerGroup <- function(gr, tx, downstream = 0L, upstream = 0L) {
  g <- ORFik:::asTX(gr, tx)
  
  starts <- pmax(start(g) - upstream, 1L)
  indices <- chmatch(txNames(gr), names(tx))
  if (upstream != 0) {
    ends <- pmin(end(g) + downstream, widthPerGroup(tx[indices], FALSE))
    ranges(g) <- IRanges(starts, ends)
  } else {
    starts(g) <- starts
  }
  
  return(ORFik:::pmapFromTranscriptF(g, tx, indices))
}

riboTISCoverageProportion <- function(grl, tx, footprints,
                                      onlyProportion = FALSE, average = FALSE,
                                      pShifted = TRUE, keep.names = FALSE,
                                      upStart = if (pShifted) 5 else 20,
                                      downStop = if (pShifted) 20 else 5) {
  windowSize <- upStart + downStop + 1
  window <- windowPerGroup(startSites(grl, TRUE, TRUE, TRUE), tx, downStop,
                           upStart)
  noHits <- widthPerGroup(window) < windowSize
  if (all(noHits)) {
    warning("no grl had valid window size!")
    return(RleList(Rle(values = 0, lengths = windowSize)))
  }
  window <- window[!noHits]
  # fix names, find a better way to store this, should be a function.
  if (grep("_", names(grl[1]))) {
    names(window) <- names(grl[!noHits])
    g <- unlist(window, use.names = TRUE)
    names(g) <- sub("\\..*", "", names(g), perl = TRUE)
    mcols(g) <- DataFrame(row.names = names(g), names = names(g))
    window <- groupGRangesBy(g)
  }
  
  unlTile <- tile1(window, matchNaming = FALSE)
  if(length(unlTile) != length(window)) stop("Bad naming, most likely _
                                             is not used for ORF correctly")
  unlTile <- unlistGrl(unlTile)
  
  rwidth <- ORFik:::readWidths(footprints)
  footprints <- footprints[rwidth < 31 & rwidth > 26]
  rwidth <- ORFik:::readWidths(footprints)
  allLengths <- sort(unique(rwidth))
  if (!(ORFik:::is.gr_or_grl(footprints) & unique(width(footprints)) == 1)) {
    footprints <- resize(granges(footprints), 1)
  }
  
  lengthProportions <- c()
  for (l in allLengths) {
    ends_uniq <- footprints[rwidth == l]
    
    cvg <- overlapsToCoverage(unlTile, ends_uniq, FALSE, type = "within")
    
    cvg <- cvg /sum(cvg)
    cvg[is.nan(unlist(sum(runValue(cvg)), use.names = FALSE))] <-
      RleList(Rle(values = 0, lengths = windowSize))
    if (average) {
      cvg <- Reduce(`+`, cvg)
      lengthProportions <- c(lengthProportions, as.vector(cvg/sum(cvg)))
    } else {
      lengthProportions <- c(lengthProportions, cvg)
    }
  }
  if (!onlyProportion) {
    if (average) {
      lengths <- unlist(lapply(allLengths, function(x){rep.int(x,windowSize)}),
                        use.names = FALSE)
      df <- data.frame(prop = lengthProportions, length = lengths,
                       pos = rep(seq.int(-upStart, downStop),
                                 length(allLengths)))
      return(df)
    } else {
      warning("Can only return proportion when average == FALSE")
      return(lengthProportions)
    }
  } else {
    if (keep.names & !average) {
      names(lengthProportions[[1]]) <- names(window)
    }
    return(lengthProportions)
  }
}
