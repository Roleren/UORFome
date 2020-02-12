# write out leaders
rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/") 
source("./pipelineSetup.R")
dbDisconnect(uorfDB)
rm(uorfDB)

source("./TempScripts/tcp_pipeline.R")

setwd(p(mainFolder, "/AdamVienna/"))
plotFolder <- "/export/valenfs/projects/Håkon/AdamVienna/plots/new_plots/"

gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath);
getFasta("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa")
seqlevelsStyle(txdb)  <- seqlevelsStyle(fa)
leaders <- fiveUTRsByTranscript(txdb, use.names = T)
leaders <- leaders[widthPerGroup(leaders) > 200]
a <- extractTranscriptSeqs(fa, transcripts = leaders)
writeXStringSet(a, filepath = "leader.fasta")

df <- getTCPdf()



# Now do perl scripts
#system(p("./Adams_vienna_explain.sh -f ./leader.fasta -o ", getwd()))

########################## READ VIENNA OUTPUT #############################

hits <- fread(p(getwd(),"/csvs/am.csv"), fill = TRUE, nrow = length(leaders), header = FALSE)
if(nrow(hits) != length(leaders)) stop("did not create all")
h <- hits[,-1]
m <- setDT(melt(t(h)))
best <- m[, .(which.min(value)), by = Var2]
uniques <- m[, .(unique(value)), by = Var2]
positions <- IRanges(best$V1, width = 1)

# 50 NT
m <- setDT(melt(t(h[,1:11])))
means <- m[, .(median = median(value, na.rm = TRUE), mean = mean(value, na.rm = TRUE)), by = Var2]
means[, txNames := names(leaders)[Var2]]
means <- means[order(mean),]
means[, Var2 := NULL]
write.csv(means, "viennaLeaders_50NT.csv")

### find best per whole
means <- m[, .(median = median(value, na.rm = TRUE), mean = mean(value, na.rm = TRUE)), by = Var2]
means[, txNames := names(leaders)[Var2]]
means <- means[order(mean),]
means[, Var2 := NULL]
write.csv(means, "viennaLeaders.csv")

#means <- read.csv("viennaLeaders.csv", row.names = 1)
high <- means[mean < quantile(mean, 0.5),]
low <- means[mean >= quantile(mean, 0.5),]
highest <- means[mean < quantile(mean, 0.20),]
lowest <- means[mean >= quantile(mean, 0.80),]



########################## WORK ON STRUCTURE #############################
tx <- exonsBy(txdb, use.names = TRUE)
grl <- ORFik:::pmapFromTranscriptF(positions, tx[hits$V1], removeEmpty = TRUE)

txHigh <- high$txNames
txLow <- low$txNames
txHighest <- highest$txNames
txLowest <- lowest$txNames

lowGrl <- grl[names(grl) %in% txLowest][start(positions[names(grl) %in% txLowest]) > 76]
lowGrl <- lowGrl[!(names(lowGrl) %in% "ENSDART00000088469")] # filter out bad
highGrl <- grl[names(grl) %in% txHighest][start(positions[names(grl) %in% txHighest]) > 76]
# "ENSDART00000088469"
regionWindowAll(lowGrl, highGrl, tx, outdir = p(plotFolder, "all_"), df = df, title = "Leaders metacoverage", 
                scores = "fracPos")

df1 <- df
df1$RFP <- NULL
regionWindowAll(lowGrl, highGrl, tx, outdir = p(plotFolder, "new_all_"), df = df1, title = "Leaders metacoverage")

regionWindow(grl, tx, outdir = p(plotFolder, "whole_"), df = df1, title = "Leaders metacoverage")
# Where are the structures

leaderSSU <- scaledWindowPositions(leaders, startSites(grl, asGR = TRUE))
leaderSSU[, `:=` (fraction = "structures", feature = "leaders")]

outName <- paste("structure_sum", ".pdf", sep="")
windowCoveragePlot(leaderSSU, output = outName, scoring = "sum")
outName <- paste("structure_zscore", ".pdf", sep="")
windowCoveragePlot(leaderSSU, output = outName, scoring = "zscore")


### Sequence info

windowsStart <- startRegion(grl[start(positions) > 76], tx, TRUE, upstream = 75, downstream = 74)
seqs <- extractTranscriptSeqs(fa, windowsStart)

cons <- Biostrings::consensusMatrix(seqs, as.prob=TRUE)

dt <- melt(data.table(cons))
maxs <- dt[, .(which.max(value)), by = variable]
letts <- c("A", "C", "G", "T")[maxs$V1]
library(ggplot2)
library(ggseqlogo)

plot <- ggseqlogo(data = as.character(seqs), method = "prob") + theme_logo() + scale_x_discrete(0)
ggsave(filename = "logo.pdf", plot = plot, width = 40)

### uORF distribution
cds <- cdsBy(txdb, use.names = TRUE)[names(leaders)]

uorfRegion <- ORFik:::addCdsOnLeaderEnds(leaders, cds)
seqs1 <- ORFik:::txSeqsFromFa(uorfRegion, fa, TRUE)
orfs <- ORFik::findMapORFs(uorfRegion, seqs1, startCodon = "ATG", groupByTx = F)
ret <- filterORFs(orfs)

orfsLoc <- scaledWindowPositions(leaders, startSites(ret, asGR = T))
orfsLoc[, `:=` (fraction = "uorfLocations", feature = "leaders")]

outName <- paste(plotFolder, "test_uorfLocations", ".pdf", sep="")
windowCoveragePlot(orfsLoc, output = outName, scoring = "sum", title = "uorf (ATG) location ")

orfsLocL <- scaledWindowPositions(leaders, stopSites(ret, asGR = T))
orfsLocL[, `:=` (fraction = "uorfLocations", feature = "leaders")]
orfsLocC <- scaledWindowPositions(cds, stopSites(ret, asGR = T))
orfsLocC[, `:=` (fraction = "uorfLocations", feature = "cds")]

cov <- rbindlist(list(orfsLocL, orfsLocC))

outName <- paste(plotFolder, "test_uorfLocations_stop_tx", ".pdf", sep="")
windowCoveragePlot(cov, output = outName, scoring = "sum", title = "uorf stopSite location ")

## seperation on structure
uTx <- unique(txNames(ret))
# without uORFs
regionWindowAll(highGrl[!(names(highGrl) %in% uTx)], highGrl[names(highGrl) %in% uTx], tx,
                outdir = p(plotFolder, "t_highest_uoRFs_"), df = df1, title = "Leaders metacoverage", seperations = c("noUORFs", "uORFs"))
regionWindowAll(lowGrl[!(names(lowGrl) %in% uTx)], lowGrl[names(lowGrl) %in% uTx], tx,
                outdir = p(plotFolder, "t_lowest_uoRFs_"), df = df1, title = "Leaders metacoverage", seperations = c("noUORFs", "uORFs"))
# with uORFs
regionWindowAll(lowGrl[!(names(lowGrl) %in% uTx)], highGrl[!(names(highGrl) %in% uTx)], tx,
                outdir = p(plotFolder, "t_test_noUORFs_"), df = df1, title = "Leaders metacoverage", seperations = c("l_noUORFs", "h_noUORFs"))
regionWindowAll(lowGrl[(names(lowGrl) %in% uTx)], highGrl[(names(highGrl) %in% uTx)], tx,
                outdir = p(plotFolder, "t_withUORFs_"), df = df1, title = "Leaders metacoverage", seperations = c("l_UORFs", "h_UORFs"))

regionWindowAll(lowGrl[!(names(lowGrl) %in% uTx)], highGrl[!(names(highGrl) %in% uTx)], tx,
                outdir = p(plotFolder, "t_All_both_"), df = getTCPdf3(), title = "Leaders metacoverage", seperations = c("l_noUORFs", "h_noUORFs"))
regionWindowAll(lowGrl[(names(lowGrl) %in% uTx)], highGrl[(names(highGrl) %in% uTx)], tx,
                outdir = p(plotFolder, "t_All_both_with_"), df = getTCPdf3(), title = "Leaders metacoverage", seperations = c("l_UORFs", "h_UORFs"))
### Create initiation rate and TE

countsPerLibraryOverTranscriptPerSubunit()

### Minimum free energy meta

cov <- m[,-1]
colnames(cov) <- c("genes", "count")
cov <- cov[!is.na(cov$count),]
perGroup <- cov[,.N, by = genes]$N

cov[, ones := rep.int(1L, length(genes))]
cov[, position := cumsum(ones), by = genes]
cov$ones <- NULL


scaleTo = 100
scoring = "meanPos"

cov[, scalingFactor := (scaleTo/widthPerGroup(leaders, FALSE))[genes]]
cov[, position := ceiling(scalingFactor * position)]

cov[position > scaleTo]$position <- scaleTo
groupFPF <- quote(list(genes, position))
res <- cov[, .(score = mean(count, na.rm = TRUE)), by = eval(groupFPF)]

windowCoveragePlot(res, output = "MFE_distribution.png", scoring = "sum", title = "Minimum Free Energy distribution")

### Ribo-seq test

txnames <- filterTranscripts(txdb, 201, 100, 100)
leaders <- fiveUTRsByTranscript(txdb, use.names = T)[txnames]
cds <- cdsBy(txdb,"tx", use.names = TRUE)[txnames]
trailers <- threeUTRsByTranscript(txdb, use.names = T)[txnames]
transcriptWindow(leaders, cds, trailers, df, outdir = "/export/valenfs/projects/uORFome/AdamVienna/plots/new_plots/SSU_RIBO-SEQ_",
                 fractions = c("SSU", "Ribo-seq"))

# Ribo-seq test with structure, differentiate into groups of structure and non structure
txHigh <- txnames[txnames %in% high$txNames]
txLow <- txnames[txnames %in% low$txNames]
txHighest <- txnames[txnames %in% highest$txNames]
txLowest <- txnames[txnames %in% lowest$txNames]

transcriptWindow(leaders[txHigh], cds[txHigh], trailers[txHigh], df, outdir = "/export/valenfs/projects/uORFome/AdamVienna/plots/new_plots/HighStruct_SSU_RIBO-SEQ_",
                 fractions = c("SSU", "Ribo-seq"))
transcriptWindow(leaders[txLow], cds[txLow], trailers[txLow], df, outdir = "/export/valenfs/projects/uORFome/AdamVienna/plots/new_plots/LowStruct_SSU_RIBO-SEQ_",
                 fractions = c("SSU", "Ribo-seq"))

transcriptWindow(leaders[txHighest], cds[txHighest], trailers[txHighest], df, outdir = "/export/valenfs/projects/uORFome/AdamVienna/plots/new_plots/HighestStruct_SSU_RIBO-SEQ_",
                 fractions = c("SSU", "Ribo-seq"))
transcriptWindow(leaders[txLowest], cds[txLowest], trailers[txLowest], df, outdir = "/export/valenfs/projects/uORFome/AdamVienna/plots/new_plots/LowestStruct_SSU_RIBO-SEQ_",
                 fractions = c("SSU", "Ribo-seq"))




 
