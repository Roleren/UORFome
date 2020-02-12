# Analysis of STOP region SSU usage with RCP seq experiments

# INFO: 
# Slides active: https://docs.google.com/presentation/d/1rTzd3OqAlRu1x6BsdKXYbYp2jFLN3a_NPMIdyWp481c/edit#slide=id.g6c1adc1b15_0_46
# Slides Eivind: https://docs.google.com/presentation/d/1tWDaN-8i_NhzEIT5yhpvTQKGDOe_P_f0ksgYfoagIDM/edit#slide=id.g469f58b04e_0_42

################################## Package loadings ##############################################################
source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")
source("./microRNA/MZdicer.R")

stop <- "/export/valenfs/projects/Hakon/mir430/Stopcodon/"
sr <- paste0(stop, "stopRegion/")
sm <- paste0(stop, "mrnaRegions/")
st <- paste0(stop, "trailerRegions/")

# Pick which dataset, split in RNA seq and rest
list.experiments()
# Add merged RNA here TODO !!!!!!!!
df <- read.experimentl("rcp-experiment_WTvsMZ") # RNA
dfr <- read.experimentl("Val19") # New RNA seq

df <- df[df$condition != "MZ" & df$libtype != "RNA",]
dfr <- dfr[dfr$condition == "",] # Pick stage 2 and 6 of RNA seq
#dfr@listData$stage <- c("2", "6")


# Get annotation
fa <- ORFik:::findFa(df)
txdb <- loadTxdb(df, seqlevelsStyle(fa))
loadRegions(txdb)

############################ FILTERS #######################################
# Length per transcript region, and longest transcript per gene ?
txNames <- filterTranscripts(txdb, 76, 76, 111, longestPerGene = FALSE)
# Remove strange scaffold, that will crash things
txNames <- txNames[txNames %in% names(startSites(mrna, keep.names = TRUE) > 75)]


# RNA filter
mrnaFPKM <- makeSummarizedExperimentFromBam(dfr, geneOrTxNames = "tx", type = "fpkm",
                                            region = mrna, saveName = p(stop, "rna_countTable.rds"))
# NOTE: LEE2013 RNA
trailerFPKM <- makeSummarizedExperimentFromBam(df, geneOrTxNames = "tx", type = "fpkm", 
                                               region = trailers, saveName = p(stop, "trailer_LSU_SSU countTable.rds"))
cdsFPKM <- makeSummarizedExperimentFromBam(df, geneOrTxNames = "tx", type = "fpkm", 
                                           region = cds, saveName = p(stop, "cds_countTable.rds"))
#keep <- (rowSums(mrnaFPKM) > 1)
minRNASeq <- 1
validMRNA <- rownames(mrnaFPKM)[rowMins(as.matrix(mrnaFPKM)) > minRNASeq]
txNames <- txNames[txNames %in% validMRNA]

# Find valid Stop contexts
stopContext <- ORFik:::startRegionString(trailers[txNames], mrna, fa, upstream = 3, 0)
table(stopContext)
stopContext <- stopContext[stopContext %in% names(table(stopContext))[table(stopContext) > 50]]
table(stopContext)
length(table(stopContext))
validStopCodons <- names(stopContext)
txNames <- txNames[txNames %in% validStopCodons]

print(paste("Will keep", length(txNames), "of the transcripts."))
print(paste("Ratio:", round(length(txNames) / length(mrna), 3)))

############################ UPDATE ANNOTATION TO FILTER ###################
splitRegions(splitList = list(txNames))
# Remove region of 3' bias in trailers 
trailers <- ORFik:::upstreamOfPerGroup(trailers, upstreamOf = stopSites(trailers, is.sorted = TRUE) + ifelse(strandBool(trailers), -35, 35))

mrnaFPKM <- mrnaFPKM[rownames(mrnaFPKM) %in% txNames,]
trailerFPKM <- trailerFPKM[rownames(trailerFPKM) %in% txNames,]
cdsFPKM <- cdsFPKM[rownames(cdsFPKM) %in% txNames,]
############################# PLOTS ########################################
### STOP CODON REGION
for (i in unique(stopContext)) {
  group <- names(stopContext[stopContext == i])
  regionWindow(trailers[group], mrna[group], df, outdir = sr, title = paste("Stop codon:", i), dfr = dfr)
}

### All mRNA regions
for (i in unique(stopContext)) {
  group <- names(stopContext[stopContext == i])
  transcriptWindow(leaders[group], cds[group], trailers[group], df, outdir = sm, 
                   title = paste("Stop codon:", i), dfr = dfr, windowSize = 76, idName = i)
}

### 3' UTR only
splitMetacoverage(trailers, stopContext, df, outdir = st, scores = c("sum", "zscore"),
                  colors = c(rep("skyblue4", 4), rep("orange", 4)), allTogether = TRUE,
                  title = "trailer metacoverage", dfr = NULL)

################################# GROUP ANALYSIS STATISTICS #########################
# Scanning effiency per stop codon
trailer.melt <- melt(trailerFPKM, variable.name = "fraction", value.name = "score")
trailer.melt <- trailer.melt[grep("SSU_", x = fraction), ]
scanning <- trailer.melt
# Normalize by RNA
rna.melt <- melt(mrnaFPKM, variable.name = "fraction", value.name = "score")
scanning <- scanning[, score := score / rna.melt$score]

scanning[, bases := rep(stopContext, length(unique(fraction)))]
scanning[, genes := rep(seq(nrow(trailerFPKM)), length(unique(fraction)))]
# puRine: R (A, G), pYrimidine: Y (C, T)
scanning[, feature := paste0(substring(bases, 1, 3),
                             ifelse(substring(bases, 4) %in% c("A", "G"),
                                    "R", "Y"))]
# Find if any stop codon group is different
difGroupAnalysis(scanning, stop, "Scanning Efficiency", "Stop Codon")

# 3' UTR SSU vs cds LSU per stop codon
trailer.melt <- melt(trailerFPKM, variable.name = "fraction", value.name = "score")
trailer.melt <- trailer.melt[grep("SSU_", x = fraction), ]
scanning <- trailer.melt
# Normalize by cds LSU
cds.melt <- melt(cdsFPKM, variable.name = "fraction", value.name = "score")
cds.melt <- cds.melt[grep("LSU_", x = fraction), ]
scanning <- scanning[, score := score / cds.melt$score]

scanning[, bases := rep(stopContext, length(unique(fraction)))]
scanning[, genes := rep(seq(nrow(trailerFPKM)), length(unique(fraction)))]
# puRine: R (A, G), pYrimidine: Y (C, T)
scanning[, feature := paste0(substring(bases, 1, 3),
                             ifelse(substring(bases, 4) %in% c("A", "G"),
                                    "R", "Y"))]
difGroupAnalysis(scanning, stop, "SSUFPKM per LSUFPKM", "Stop Codon")
