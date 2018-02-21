
# database plots

# plot info :
# 1. general stat, median, mean, count, max, mean width, max width, min width, etc.
# 2. Cage per tissue, difference (counts per tissue), unique per tissue, variance per uorf on cage
# 3. ribo seq, how does the cage change with fpkm support ? 
# 4. Expl.

rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./DataBaseCreator.R")


#1. general statistics:

#check max/min set of uorfs, check that everything looks good


#find how few orfs are minimum
rm(fiveUTRs)
getLeaders()
numberOfTx <- length(fiveUTRs) 

widthsMeanOriginal <- mean(ORFik:::widthPerGroup(fiveUTRs,F))
originalUorfsByTx <- getUnfilteredUORFs(fiveUTRs, assignRanges = F)

gr <- unlist(originalUorfsByTx, use.names = F)
originalUorfs <- groupGRangesBy(gr, gr$names)
originalUorfs <- removeORFsWithinCDS(originalUorfs)
gr <- unlist(originalUorfs, use.names = F)
originalUorfsByTx <- groupGRangesBy(gr, names(gr))
nTxWithUorfs <- length(originalUorfsByTx)
uorfTxRatio <- nTxWithUorfs/numberOfTx
nOrfsOriginal <- length(originalUorfs)
uorfsPerTx <- nOrfsOriginal/numberOfTx
uorfsPerTxNormalized <- uorfsPerTx/uorfTxRatio


#find how many orfs are maximum
rm(fiveUTRs)
fiveUTRsMax <- leaderAllSpanning()

widthsMeanMax <- mean(ORFik:::widthPerGroup(fiveUTRsMax,F))
maxUorfsByTx  <- getUnfilteredUORFs(fiveUTRsMax, assignRanges = F)
nTxWithUorfsMax <- length(maxUorfsByTx)
uorfTxRatioMax <- nTxWithUorfsMax/numberOfTx
gr <- unlist(maxUorfsByTx, use.names = F)
maxUorf <- groupGRangesBy(gr, gr$names)
nOrfsMax <- length(maxUorf)
uorfsPerTxMax <- nOrfsMax/numberOfTx
uorfsMaxPerTxNormalized <- uorfsPerTxMax/uorfTxRatioMax


  
#now data-base results:
grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)

#first some sanity checks
# ORfs found can not be more than maxOrfs
if (length(grl) > nOrfsMax) stop("number of orfs in atlas > than maxOrfs possible!")
# ORfs found must all be in uidMax
uidGRL <- toUniqueIDFromGR(grl)
if (length(uidGRL) != length(unique(uidGRL))) stop("uids as gr was not unique, something is wrong!")
uidMax <- toUniqueIDFromGR(maxUorf)
uidMax <- unique(uidMax)
overlap <- uidGRL %in% uidMax
uidOriginal <- toUniqueIDFromGR(originalUorfs)
nuorfsOriginalUnique <- length(unique(uidOriginal)) 
if (sum(overlap) +1  < min(length(uidGRL),length(uidMax))) stop("Not all orfs in atlas are in maxOrfs")
#whats up with that one that is not there ????

    
#width per uorf
widths <- ORFik:::widthPerGroup(grl,F)
medianWidths <- median(widths)
meanWidths <- mean(widths)
maxWidth <- max(widths)
minWidth <- min(widths)

#width per exon
widths <- width(unlist(grl, use.names = F))

exonmedianWidths <- median(widths)
exonmeanWidths <- mean(widths)
maxWidth <- max(widths)
minWidth <- min(widths)

#exon stats
numberOfexons <- length(widths)
exonsPerGroup <- GroupGRangesExonsPerGroup(grl)
medianExons <- median(exonsPerGroup)
meanExons <- mean(exonsPerGroup)

#some max checks
a <- findOverlaps(grl,fiveUTRsMax, type = "within")
index <- which.max(runLength(Rle(to(a))))
leaderWithMostUORFs <- fiveUTRsMax[index]

#check width changes from cage, the effect


df <- data.frame()

library(ggplot2)

#make table of stats



#plots: 

#1. width histogram

widthsGRL <- ORFik:::widthPerGroup(grl,F)
df <- data.frame(widths = widthsGRL)
lengthPlot <- ggplot() +
  geom_histogram(data = df, aes(x = log10(widths)),color="orange") +
  xlab("length of uorfs (log)") + ylab("number of counts")  +
  labs(title="Length of UORFs")
lengthPlot


#2. cage average hits per tissue

tissueAtlas <- readTable("tissueAtlasByCage")

colsums <- colSums(tissueAtlas[,2:length(tissueAtlas)])
badIds <- which(colsums == 0) +1 # <-no matches found for tissue in cage
names(badIds) <- NULL
tissueAtlas[,badIds] <- NULL
colsums <- colsums[colsums > 0]
df <- data.frame(colsums = colsums)
line <- data.frame(y = c(nuorfsOriginalUnique,nuorfsOriginalUnique), x = c(0, nrow(df)) )
distributionPlot <- ggplot()  +
  geom_histogram(data = df, aes(x = 1:length(colsums), y = colsums),color="orange",stat="identity") +
  xlab("tissue") + ylab("uorfs in tissue")  +
  labs(title="cage hits of UORFs")  +  geom_line(aes(x = line$x, y = line$y), color = "blue")
distributionPlot

variation <- colsums - nuorfsOriginalUnique

uorfsInAllTissues <- sum(rowSums(tissueAtlas[
  ,2:length(tissueAtlas)]) == (length(tissueAtlas)-1))
uorfsRatioTissue <- uorfsInAllTissues/nrow(tissueAtlas)

uorfsToUniqueTissue <- rowSums(tissueAtlas[
  ,2:length(tissueAtlas)]) == 1

uorfsToUniqueTissueSum <- sum(uorfsToUniqueTissue)

#3. Show that fiveUTRs acctually gets reduced in many samples
# stacked bar for 3 points, unchanged, down, up
load("leaderWidthChanges.rdata")
df <- as.data.frame(dt)
widthChangeCagePlotSame <- ggplot()  +
  geom_histogram(data = df, aes(x = 1:nrow(df), y = same),color="orange",stat="identity") +
  xlab("cage experiment") + ylab("width")  +
  labs(title="width changes per cage experiment: same")  
widthChangeCagePlotSame

widthChangeCagePlotBigger <- ggplot()  +
  geom_histogram(data = df, aes(x = 1:nrow(df), y = bigger),color="orange",stat="identity") +
  xlab("cage experiment") + ylab("width")  +
  labs(title="width changes per cage experiment: bigger") 
widthChangeCagePlotBigger

widthChangeCagePlotSmaller <- ggplot()  +
  geom_histogram(data = df, aes(x = 1:nrow(df), y = smaller),color="orange",stat="identity") +
  xlab("cage experiment") + ylab("width")  +
  labs(title="width changes per cage experiment: smaller") 
widthChangeCagePlotSmaller

widthChangeCagePlotMeanBigger <- ggplot()  +
  geom_histogram(data = df, aes(x = 1:nrow(df), y = meanBigger),color="orange",stat="identity") +
  xlab("cage experiment") + ylab("width")  +
  labs(title="width changes per cage experiment: mean of bigger") 
widthChangeCagePlotMeanBigger


# plot ribo support differences
ribo <- readTable("RiboByTissueTF")
riboColSums <- colSums(ribo[,-1])
df.ribo <- as.data.frame(riboColSums)

riboSupportTissue <- ggplot()  +
  geom_histogram(data = df.ribo, aes(x = rownames(df.ribo), y = riboColSums),color="orange",stat="identity") +
  xlab("ribo tissue") + ylab("uorfs supported")  +
  labs(title="Ribo-seq support, number of uorfs per tissue")  
riboSupportTissue

# case of brain, which of the ribo-seq supported reads
# are not in standard transcriptome ?


brainRibo <- ribo$brain
brainCage <- tissueAtlas$brain
supportedUorfs <- sum(brainRibo)
indeces <- which(brainRibo == 1)
idsSupported <- ribo$uorfID[indeces]
uniqueSupportedForBrain <- idsSupported %in% uidOriginal
idsResult <- sum(uniqueSupportedForBrain)



# get unique brain grl
# find genes to check for interesting features
uidsBrain <- idsSupported[!uniqueSupportedForBrain]
grlBrain <- toGRFromUniqueID(uidsBrain)
overlaps <-  findOverlaps(grlBrain, fiveUTRsMax, type = "within")
txname <- names(fiveUTRsMax[to(overlaps)])

getAllTranscriptLengths()
geneID <- allLengths$gene_id[allLengths$tx_name == txname[2]]


# now find a brain example to use:
# using HTR3A gene, this made no hit, so not valid

allLengths[allLengths$gene_id == "ENSG00000166736",]
findOverlaps(grl, fiveUTRsMax["ENST00000504030"], type = "within")
grl[27154] # <- verified as uorf
findOverlaps(grlBrain, fiveUTRsMax["ENST00000504030"], type = "within")

# lets try huntington gene HTT
# ENSG00000197386
allLengths[allLengths$gene_id == "ENSG00000197386",]
findOverlaps(grlBrain, fiveUTRsMax["ENST00000355072"], type = "within")
# no hit, lets try with all allowed
grlBrainAll <- toGRFromUniqueID(idsSupported)
findOverlaps(grlBrainAll, fiveUTRsMax["ENST00000355072"], type = "within")
# found a match, but is not unique to brain
indecesKid <- which(ribo$kidney == 1)
idsSupportedKid <- ribo$uorfID[indecesKid]
grlKid<- toGRFromUniqueID(idsSupportedKid)
findOverlaps(grlKid, fiveUTRsMax["ENST00000355072"], type = "within")

clustering <- function(){
  # Ribo seq clustering of 4 groups
  # gonzales 1:5, kidney, 18:25, remember +1 for uorfid!
  riboAtlas <- readTable("riboAll")
  
  usedSamples <- riboAtlas[, .(V2,V3,V4,V5,V6, V19,V20,V21,V22,V23,
                              V49,V50,V51,V52,V53, V7,V8,V9,V10,V11 )]
  names1 = paste("Brain",1:5)
  names2 = paste("Kidney",1:5)
  names3 = paste("Blood",1:5)
  names4 = paste("FibroBlast",1:5)
  colnames(usedSamples)[1:5] = names1
  colnames(usedSamples)[6:10] = names2
  colnames(usedSamples)[11:15] = names3
  colnames(usedSamples)[16:20] = names4
  binDist <- dist(x = t(usedSamples), method = "euclidean")
  hclustResult <-  hclust(binDist)
  plot(hclustResult)
  
  # cage seq of clustering of 4 groups
  load("./UORFAtlas.rdata")
  #blood: 284-288, brain: 562-566, kidney: 857-861, fibro(gum): 811,812,814, 815,817
  usedCage <- uorfAtlas[, .(`562`,`563`,`564`,`565`,`566` ,`857`,`858`,`859`,`860`,`861`
                            ,`284`,`285`,`286`,`287`,`288` ,`811`,`812`,`814`,`815`,`817`)]
  colnames(usedCage)[1:5] = names1
  colnames(usedCage)[6:10] = names2
  colnames(usedCage)[11:15] = names3
  colnames(usedCage)[16:20] = names4
  binDist <- dist(x = t(usedCage), method = "binary")
  hclustResult <-  hclust(binDist)
  plot(hclustResult)
}
