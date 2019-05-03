# write out leaders
rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/") 
source("./pipelineSetup.R")
source("./TempScripts/tcp_pipeline.R")
gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath);
getFasta("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa")
seqlevelsStyle(txdb)  <- seqlevelsStyle(fa)
leaders <- fiveUTRsByTranscript(txdb, use.names = T)
leaders <- leaders[widthPerGroup(leaders) > 200]
a <- extractTranscriptSeqs(fa, transcripts = leaders)

setwd(p(mainFolder, "/AdamVienna/"))
writeXStringSet(a, filepath = "leader.fasta")

# Now do perl scripts
#system(p("./Adams_vienna_explain.sh -f ./leader.fasta -o ", getwd()))

# Now do final analysis

# check that all was made
if(nrow(hits) != length(leaders)) stop("did not create all")
hits <- fread(p(getwd(),"/csvs/am.csv"), fill = TRUE, nrow = length(leaders), header = FALSE)
h <- hits[,-1]
m <- setDT(melt(t(h)))
best <- m[, .(which.min(value)), by = Var2]
uniques <- m[, .(unique(value)), by = Var2]
positions <- IRanges(best$V1, width = 1)

# find best per

means <- m[, .(median = median(value, na.rm = TRUE), mean = mean(value, na.rm = TRUE)), by = Var2]
means[, txNames := names(leaders)[Var2]]
means <- means[order(mean),]
means[, Var2 := NULL]
write.csv(means, "viennaLeaders.csv")

# now create window
tx <- exonsBy(txdb, use.names = TRUE)
grl <- ORFik:::pmapFromTranscriptF(positions, tx[hits$V1], removeEmpty = TRUE)
hitMapStart <- regionWindow(grl[start(positions) > 76], tx, outdir = "frame_structure_tcp_")

regionWindowBoth(grl[start(positions) > 76], tx, outdir = "structure_tcp_")

# Where are the structures

leaderSSU <- scaledWindowPositions(leaders, startSites(grl, asGR = TRUE))
leaderSSU[, `:=` (fraction = "structures", feature = "leaders")]

outName <- paste("structure_sum", ".pdf", sep="")
windowCoveragePlot(leaderSSU, output = outName, scoring = "sum")
outName <- paste("structure_zscore", ".pdf", sep="")
windowCoveragePlot(leaderSSU, output = outName, scoring = "zscore")


# Sequence info

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

# uORF distribution
cds <- cdsBy(txdb, use.names = TRUE)

uorfRegion <- ORFik:::addCdsOnLeaderEnds(leaders, cds)
seqs1 <- ORFik:::txSeqsFromFa(uorfRegion, fa, TRUE)
orfs <- ORFik::findMapORFs(leaders, seqs1, startCodon = "ATG")
ret <- filterORFs(orfs)

orfsLoc <- scaledWindowPositions(leaders, startSites(ret, asGR = T))
orfsLoc[, `:=` (fraction = "uorfLocations", feature = "leaders")]

outName <- paste("uorfLocations", ".pdf", sep="")
windowCoveragePlot(orfsLoc, output = outName, scoring = "sum", title = "uorf (ATG) location ")

# Create initiation rate and TE

countsPerLibraryOverTranscriptPerSubunit()

# Minimum free energy meta

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

windowCoveragePlot(res, output = "MFE_distribution.pdf", scoring = "sum", title = "Minimum Free Energy distribution")

# Ribo-seq test
getTCPdf(LSU = c("/export/valenfs/projects/uORFome/AdamVienna/Ribo_seq/MZ_64/aligned/zebrafish_bazzini2012_2hAligned.sortedByCoord.out.bam",
                 "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/256Cell_trimmed.bam",
                 "/export/valenfs/projects/uORFome/AdamVienna/Ribo_seq/MZ_Shield/aligned/zebrafish_bazzini2012_6hAligned.sortedByCoord.out.bam",
                 "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Shield_trimmed.bam"))



ga <- readGAlignments()
txnames <- filterTranscripts(txdb, 201, 100, 100)
leaders <- fiveUTRsByTranscript(txdb, use.names = T)[txnames]
cds <- cdsBy(txdb,"tx", use.names = TRUE)[txnames]
trailers <- threeUTRsByTranscript(txdb, use.names = T)[txnames]
transcriptWindowPer(leaders, cds, trailers, read1 = ga, fractions = "Ribo-seq", outdir = p(mainFolder,"/ribo-test_"))


# Ribo te check

fpkms <- m[, .(median = median(value, na.rm = TRUE), mean = mean(value, na.rm = TRUE)), by = Var2]
fpkms$fpkm <- ORFik:::fpkm(leaders, ga)
cor.test(fpkms$mean, fpkms$fpkm) # No correlation found
