#STEP 1: Load libraries needed
source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")

#STEP 2: Load annotation, and filter out bad transcripts
# This is Adams CAGE for shield -> (he also uses this for most other stuff)
gtfPathShield <- "/export/valenfs/projects/adam/TCP_seq/transcript_GFF3/shield_transcripts.gff3"
txdbShield <- loadTxdb(gtfPathShield) # <- A transcript database of the gtf
validTxShield <- ORFik:::filterTranscripts(txdbShield, minFiveUTR = 100, minCDS = 75) # <- Filter

# Load mRNA transcripts, leaders, cds, and trailers seperatly
txShield <- exonsBy(txdbShield, use.names = TRUE)
leadersShield <- fiveUTRsByTranscript(txdbShield, use.names = T)
leadersShield <- leadersShield[validTxShield]
leadersShield <- leadersShield[-which(startSites(leadersShield) < 55)]
cdsShield <- ORFik:::loadRegion(txdbShield, "cds")[validTxShield]
trailersShield <- ORFik:::loadRegion(txdbShield, "trailer")
trailersShield <- trailersShield[names(trailersShield) %in% validTxShield]


#STEP 3: Load bam file and choose whole read, 5' or 3'
readsLSU4Ei <- readGAlignments("/export/valenfs/projects/uORFome/withrRNA/aligned/64_LSU_V12_4Ei.bam")
fivePrime <- ORFik:::convertToOneBasedRanges(readsLSU4Ei, method = "5prime", addSizeColumn = TRUE)
threePrime <- ORFik:::convertToOneBasedRanges(readsLSU4Ei, method = "3prime", addSizeColumn = TRUE)


#STEP 4: Make heatmap of TIS (change readsLSU4Ei with fivePrime in windowPerReadLength, to get 5prime ends etc)
scoring <- "sum" # change scoring to zscore for zscore etc
dt <- windowPerReadLength(cdsShield, txShield, readsLSU4Ei, upstream = 75, downstream = 74, 
                          zeroPosition = 75, scoring = scoring) 
dt$score <- log2(dt$score) # <- if you want log2 scores
plotTIS <- ORFik:::coverageHeatMap(coverage = dt, scoring = scoring)


#STEP 5: Make heatmap of TSS (change readsLSU4Ei with fivePrime in windowPerReadLength, to get 5prime ends etc)
scoring <- "sum" # change scoring to zscore for zscore etc
dt <- windowPerReadLength(leadersShield, extendLeaders(leadersShield, 30), readsLSU4Ei, upstream = 30, downstream = 74, 
                          scoring = scoring) 
dt$score <- log2(dt$score) # <- if you want log2 scores
plotTSS <- ORFik:::coverageHeatMap(coverage = dt, scoring = scoring)

# If you want to save ?
ggsave("TSS_heatmap_4Ei_LSU.png", plotTSS, dpi = 300)
