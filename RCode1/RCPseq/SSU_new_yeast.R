#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Yeast data analysis (used in article supplements)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

library(ORFikPipeline)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# DATA
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load Experiments
dfA <- read.experimentl("archer19_TCP")
dfA <- dfA[dfA$libtype == "SSU",]
dfA@expInVarName <- FALSE
dfa <- getArcher16()
dfa <- dfa[dfa$libtype == "SSU",]
dfl <- list(dfA, dfa)
names(dfl) <- c("19", "16")


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
txdb <- loadTxdb(dfA)
txNames <- filterTranscripts(txdb, 60, 100, 60, T)
loadRegions(txdb, names.keep = txNames)

# Output data
outFolder <- "/export/valenfs/projects/Hakon/RCP_SEQ/plots/"
outputLibs(dfl)

### CAGE ###
# CAGE DATA test
CAGE <- fread.bed("/export/valenfs/data/processed_data/CAGE/McManus_2017_yeast/aligned/merged_CAGE_5prime.bed")
CAGE_old <- fread.bed("/export/valenfs/data/processed_data/CAGE/wery_2015_S_cerevisiae/final_results/aligned_R64_1_1/merged_CAGE_5prime.bed")
leadersCage <- reassignTSSbyCage(leaders, CAGE, removeUnused = F, cageMcol = F)
leadersCage <- leadersCage[widthPerGroup(leadersCage) >= 60]
leadersCage_old <- reassignTSSbyCage(leaders, CAGE_old, removeUnused = F, cageMcol = F)
leadersCage_old <- leadersCage_old[widthPerGroup(leadersCage_old) >= 60]

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# PLOTS
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#


#################### Quality control
summary(readWidths(SSU_CS_r1))
#################### 1. Metacoverage
transcriptWindow(leaders, cds, trailers, df = dfl, outdir = outFolder, allTogether = TRUE)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# 2. TSS (5' and 3')
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
##################### non CAGE
interval <- 20:43
scores <- c("sum")
shifting <- c("5prime", "3prime")
up <- 70;down <- 59
upshort <- c(30, 20);downshort <- c(40, 59)
#scores <- c("log2sum")

heatMapL(leaders, extendLeaders(mrna, 76), dfl, upstream = up, downstream = down,
            outdir = outFolder, acLen = interval,
            scores = scores, shifting = "5prime", location = "TSS", skip.last = T)
heatMapL(leaders, extendLeaders(mrna, 76), dfl, upstream = 30, downstream = 74,
         outdir = outFolder, acLen = interval,
         scores = scores, shifting = "3prime", location = "TSS", skip.last = T)

#################### with CAGE
heatMapL(leadersCage, extendLeaders(leadersCage, 76), dfl, upstream = up, downstream = down,
         outdir = p(outFolder), acLen = interval,
         scores = scores, shifting = shifting, location = "TSS_CAGE_", skip.last = T)
heatMapL(leadersCage_old, extendLeaders(leadersCage_old, 76), dfl, upstream = up, downstream = down,
         outdir = p(outFolder), acLen = interval,
         scores = scores, shifting = shifting, location = "TSS_CAGE_old_", skip.last = T)
##################### filter for A16
A16_SSU_5p <- convertToOneBasedRanges(A16_SSU, method = "5prime", addSizeColumn = TRUE)
A16_SSU_3p <- convertToOneBasedRanges(A16_SSU, method = "3prime", addSizeColumn = TRUE)
l <- removeBadTxByRegion(tx = leadersCage, reads = A16_SSU_5p, upstream = -2, downstream = 40)
ll <- removeBadTxByRegion(tx = leadersCage_old, reads = A16_SSU_5p, upstream = -2, downstream = 40)

###################### Full length A16
heatMapL(l, extendLeaders(l, 76), dfl$`16`, upstream = upshort, downstream = downshort,
         outdir = p(outFolder), acLen = 18:70,
         scores = scores, shifting = shifting, location = "TSS_CAGE_fullLength_", skip.last = F)
heatMapL(ll, extendLeaders(ll, 76), dfl$`16`, upstream = upshort, downstream = downshort,
         outdir = p(outFolder), acLen = 18:70,
         scores = scores, shifting = shifting, location = "TSS_CAGE_old_fullLength_", skip.last = F, addFracPlot = FALSE)



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# 3. TIS (5' and 3')
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
heatMapL(cds, mrna, dfA, upstream = up, downstream = 74,
            outdir = outFolder, acLen = interval,
            scores = scores, shifting = shifting)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# 4. Barplots
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
####################### 16
A16_SSU_5pp <- A16_SSU_5p[readWidths(A16_SSU_5p) %in% interval]
A16_SSU_3pp <- A16_SSU_3p[readWidths(A16_SSU_3p) %in% interval]
regionBarPlot(l, extendLeaders(l, 76), A16_SSU_5p, outdir = p(outFolder, "A16_bp_A16_SSU_5prime.png"), upstream = upshort[1], downstream = downshort[1])
regionBarPlot(l, extendLeaders(l, 76), A16_SSU_3p, outdir = p(outFolder, "A16_bp_A16_SSU_3prime.png"), upstream = upshort[2], downstream = downshort[2])
regionBarPlot(ll, extendLeaders(ll, 76), A16_SSU_5p, outdir = p(outFolder, "A16_bp_A16_SSU_5prime_oldC.png"), upstream = upshort[1], downstream = downshort[1])
regionBarPlot(ll, extendLeaders(ll, 76), A16_SSU_3p, outdir = p(outFolder, "A16_bp_A16_SSU_3prime_oldC.png"), upstream = upshort[2], downstream = downshort[2])
####################### 19
SSU_CS_r1_5p <- convertToOneBasedRanges(SSU_CS_r1, method = "5prime", addSizeColumn = TRUE)
SSU_CS_r1_3p <- convertToOneBasedRanges(SSU_CS_r1, method = "3prime", addSizeColumn = TRUE)
SSU_CS_r1_5p <- SSU_CS_r1_5p[readWidths(SSU_CS_r1_5p) %in% interval]
SSU_CS_r1_3p <- SSU_CS_r1_3p[readWidths(SSU_CS_r1_3p) %in% interval]
regionBarPlot(l, extendLeaders(l, 76), SSU_CS_r1_5p, outdir = p(outFolder, "A19_bp_SSU_5prime.png"), upstream = upshort[1], downstream = downshort[1])
regionBarPlot(l, extendLeaders(l, 76), SSU_CS_r1_3p, outdir = p(outFolder, "A19_bp_SSU_3prime.png"), upstream = upshort[2], downstream = downshort[2])
regionBarPlot(ll, extendLeaders(ll, 76), SSU_CS_r1_5p, outdir = p(outFolder, "A19_bp_SSU_5prime_oldC.png"), upstream = upshort[1], downstream = downshort[1])
regionBarPlot(ll, extendLeaders(ll, 76), SSU_CS_r1_3p, outdir = p(outFolder, "A19_bp_SSU_3prime_oldC.png"), upstream = upshort[2], downstream = downshort[2])

# test direction
c <- readBam("/export/valenfs/data/processed_data/CAGE/wery_2015_S_cerevisiae/final_results/aligned_R64_1_1/SRR2048394.bam")
cc <- readBam("/export/valenfs/data/processed_data/CAGE/wery_2015_S_cerevisiae/final_results/aligned_R64_1_1/SRR2048395.bam")
regionBarPlot(a, extendLeaders(a, 76), A16_SSU_5p, outdir = p(outFolder, "A16_bp_A16_SSU_5prime_oldC.png"), upstream = 1, downstream = 1)
regionBarPlot(b, extendLeaders(b, 76), A16_SSU_5p, outdir = p(outFolder, "A16_bp_A16_SSU_5prime_oldC.png"), upstream = 1, downstream = 1)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Test new bam file Archer16
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
SSU <- readBam("/export/valenfs/projects/uORFome/archer_2016_yeast_15nt_without_3nt_trim/aligned/WT_SSU_Aligned.sortedByCoord.out.bam")
SSU_5p <- convertToOneBasedRanges(SSU, method = "5prime", addSizeColumn = TRUE)
SSU_3p <- convertToOneBasedRanges(SSU, method = "3prime", addSizeColumn = TRUE)
lll <- removeBadTxByRegion(tx = leadersCage_old, reads = SSU_5p, upstream = -2, downstream = 40)
regionBarPlot(lll, extendLeaders(lll, 76), SSU_5p, outdir = p(outFolder, "A16_bp_A16_new_SSU_5prime_oldC.pdf"), upstream = upshort[1], downstream = downshort[1])
regionBarPlot(lll, extendLeaders(lll, 76), SSU_3p, outdir = p(outFolder, "A16_bp_A16_new_SSU_3prime_oldC.pdf"), upstream = upshort[2], downstream = downshort[2])
tcpHeatMap_single(lll, extendLeaders(lll, 76), SSU_5p, upstream = upshort[1], downstream = downshort[1],
         outdir = p(outFolder, "A16_hm_new_5prime.pdf"), acLen = 18:70,
         scores = scores, shifting = shifting[1], location = "new_SSU_TSS_CAGE_old_fullLength_", skip.last = F, addFracPlot = F)
tcpHeatMap_single(lll, extendLeaders(lll, 76), SSU_3p, upstream = upshort[2], downstream = downshort[2],
                  outdir = p(outFolder, "A16_hm_new_3prime.pdf"), acLen = 18:70,
                  scores = scores, shifting = shifting[2], location = "new_SSU_TSS_CAGE_old_fullLength_", skip.last = F, addFracPlot = F)


