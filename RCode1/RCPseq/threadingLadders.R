#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Threading ladders (SSU, LSU) on TIS / TSS in zebrafish
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO:
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# These plots are for RCP-seq paper

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Set directories
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFikPipeline)
plotFolder <- "/export/valenfs/projects/Hakon/AdamVienna/plots/new_plots/"
heatMapsFolder <- p(plotFolder,"heatmaps/");finalFolder <- p(heatMapsFolder, "final/")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

gtfPath <- "/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db"
fa <- ORFik:::findFa("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa")
txdb <- loadTxdb(gtfPath, seqlevelsStyle(fa));
loadRegions(txdb, parts = c("tx", "leaders"))

# Shield CAGE annotation
gtfPathShield <- "/export/valenfs/projects/uORFome/Annotations/Zebrafish/shield_transcripts.gff3.db"
txdbShield <- loadTxdb(gtfPathShield, tx)
validTxShield <- filterTranscripts(txdbShield, 102, 50)
leadersShield <- loadRegion(txdbShield, "leaders")[validTxShield]
leadersShield <- leadersShield[!(startSites(leadersShield) < 100)]
# Filter out 2 genes
geneNames <- ORFik:::txNamesToGeneNames(names(leadersShield), txdbShield)
validTxShield <- names(leadersShield)[!(geneNames %in% c("ENSDARG00000036180", "ENSDARG0000001479"))]
# Load shield annotation with that filter
loadRegions(txdbShield, parts = c("tx", "leaders","cds", "trailers"), extension = "Shield",
            names.keep = validTxShield)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load libraries: RCP-seq mapped data
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#### For article:
# For Figure 1 heatmaps:
allSSU <- "/export/valenfs/projects/adam/TCP_seq/RCP_files/combined/stages_merged_SSU_updated.bam"
readsSSUAll <- readBam(allSSU, tx)
reads5SSUAll <- convertToOneBasedRanges(readsSSUAll, method = "5prime", addSizeColumn = TRUE, addScoreColumn = TRUE)
reads3SSUAll <- convertToOneBasedRanges(readsSSUAll, method = "3prime", addSizeColumn = TRUE, addScoreColumn = TRUE)

pathSSU <- "/export/valenfs/data/processed_data/TCP-seq/valen_all_withrRNA/aligned/shield_V15_merged_SSU.bam"
readsSSU <- readBam(pathSSU, tx)

SSUnt150Path <- "/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_16_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi/merge/SSU_merged_fractions_12_13_14.bam"
readsSSU150NT <- readBam(SSUnt150Path, tx)
readsSSU150NT5 <- convertToOneBasedRanges(readsSSU150NT, addSizeColumn = T); readsSSU150NT3 <- convertToOneBasedRanges(readsSSU150NT, method = "3prime",  addSizeColumn = T)

readsSSUGood <- c(readsSSU, readsSSU150NT)
reads5SGood <- convertToOneBasedRanges(readsSSUGood, method = "5prime", addSizeColumn = TRUE, addScoreColumn = TRUE)
reads3SGood <- convertToOneBasedRanges(readsSSUGood, method = "3prime", addSizeColumn = TRUE, addScoreColumn = TRUE)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TSS
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#### For article:

################### All stages TSS
acLen <- 17:68
leadersShield <- removeBadTxByRegion(leadersShield, reads3SSUAll, upstream = 70, downstream = 100, 200, 100)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 100), reads5SSUAll, upstream = 100, downstream = 51, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_5prime_all.pdf"),
                  scores = "log2sum", acLen = acLen, addFracPlot = FALSE)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads3SSUAll, upstream = 50, downstream = 101, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_3prime_shield_all.pdf"),
                  scores = "log2sum", acLen = acLen, addFracPlot = FALSE)
regionBarPlot(leadersShield, extendLeaders(leadersShield, extension = 100), reads5SSUAll, upstream = 100, downstream = 51, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_5prime_shield_all_bar.pdf"))
regionBarPlot(leadersShield, extendLeaders(leadersShield, extension = 50), reads3SSUAll, upstream = 50, downstream = 101, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_3prime_shield_all_bar.pdf"))

################### 4Ei 10 uM
df <- read.experimentl("val_4Ei64", expInVarName = TRUE)[2,]
outputLibs(df, leadersShield, type = "bedo")
SSU_4Ei64 <- convertToOneBasedRanges(val_4Ei64_SSU, method = "3prime", addSizeColumn = TRUE)
leadersShield4Ei64 <- removeBadTxByRegion(leadersShield, SSU_4Ei64, upstream = -77, downstream = 82, 10, min_cutoff = 10)
tcpHeatMap_single(leadersShield4Ei64, extendLeaders(leadersShield4Ei64, extension = 100), val_4Ei64_SSU, upstream = 100, downstream = 51, outdir = paste0(heatMapsFolder, "final/SSU_4Ei64_TSS_heatmap_5prime.pdf"),
                  scores = "log2sum", acLen = acLen, addFracPlot = FALSE, shifting = "5prime")
tcpHeatMap_single(leadersShield4Ei64, extendLeaders(leadersShield4Ei64, extension = 100), val_4Ei64_SSU, upstream = 50, downstream = 101, outdir = paste0(heatMapsFolder, "final/SSU_4Ei64_TSS_heatmap_3prime.pdf"),
                  scores = "log2sum", acLen = acLen, addFracPlot = FALSE, shifting = "3prime")
regionBarPlotAll(leadersShield4Ei64, extendLeaders(leadersShield4Ei64, extension = 100), df, upstream = c(100, 50), downstream = c(51, 101), outdir = paste0(heatMapsFolder, "final/"), location = "TSS", format = ".pdf")

#################### 4Ei 0.1 uM
df <- read.experimentl("val_10")
acLen <- 18:68
#ORFik:::simpleLibs(df);ORFik:::remove.experiments(df)
outputLibs(df, leadersShield, type = "bedo")
plot5 <- tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 100), SSU_4Ei, upstream = 100, downstream = 51, outdir = p(finalFolder, "SSU_4Ei0.1um64_TSS_heatmap_5prime.pdf"),
                           scores = "log2sum", acLen = acLen, addFracPlot = FALSE, shifting = "5prime")
ggslackR(plot5)
plot3 <- tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 100), SSU_4Ei, upstream = 50, downstream = 101, outdir = p(finalFolder, "SSU_4Ei0.1um64_TSS_heatmap_3prime.pdf"),
                           scores = "log2sum", acLen = acLen, addFracPlot = FALSE, shifting = "3prime")
ggslackR(plot3)
regionBarPlotAll(leadersShield, extendLeaders(leadersShield, extension = 100), df, upstream = c(100, 50), downstream = c(51, 101), outdir = finalFolder, location = "TSS", format = ".pdf", acLen = acLen)
slackrUpload("/export/valenfs/projects/Hakon/AdamVienna/plots/new_plots/heatmaps/final/val_10_bp_TSS_SSU_4Ei_5prime_sum.pdf", channels = "#visualizations")
slackrUpload("/export/valenfs/projects/Hakon/AdamVienna/plots/new_plots/heatmaps/final/val_10_bp_TSS_SSU_4Ei_3prime_sum.pdf", channels = "#visualizations")

################################### 150 NT valen 16
# Yamilla 150 NT lib TSS  & region
counts <- fpkm(leadersShield, readsSSU150NT5)
hitMap5 <-  tcpHeatMap_single(leadersShield[counts > 2], extendLeaders(leadersShield[counts > 2], extension = 100), readsSSU150NT5, upstream = 100, downstream = 50, outdir = paste0(heatMapsFolder, "150NT_SSU_5prime_TSS_sum.pdf"), scores = "sum", acLen = 18:150)
hitMap3 <- tcpHeatMap_single(leadersShield[counts > 2], extendLeaders(leadersShield[counts > 2], extension = 50), readsSSU150NT3, upstream = 50, downstream = 100, outdir = paste0(heatMapsFolder, "150NT_SSU_3prime_TSS_sum.pdf"), scores = "sum", acLen = 18:150)
regionBarPlot(leadersShield[counts > 5], extendLeaders(leadersShield[counts > 5], extension = 100), readsSSU150NT5, upstream = 100, downstream = 51, outdir = paste0(heatMapsFolder, "150NT_SSU_5prime_TSS_bar.pdf"))
regionBarPlot(leadersShield[counts > 5], extendLeaders(leadersShield[counts > 5], extension = 50), readsSSU150NT3, upstream = 50, downstream = 101, outdir = paste0(heatMapsFolder, "150NT_SSU_3prime_TSS_bar.pdf"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Loss over leaders (for figure 2)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# For article:
leaders200 <- leadersShield[widthPerGroup(leadersShield) >= 202]
windowsOne <- startRegion(leaders200, extendLeaders(leaders200, 100), upstream = 0, downstream = 201, is.sorted = TRUE)
hitMapTSS <- metaWindow(readsSSUAll, windowsOne, scoring = "sumPos", withFrames = FALSE,
                        fraction = "SSU", feature = "TSS", zeroPosition = 0)
windowCoveragePlot(hitMapTSS, scoring = "sum", feature = "TSS",type = "relation to transcript start (nt)",
                   output = paste0(heatMapsFolder, "final/leaderLoss_TSS.pdf"))

cds200 <- cdsShield[names(leaders200)]
windowsOne <- startRegion(cds200, txShield, upstream = 200, downstream = 101, is.sorted = TRUE)
windowsOnetest <- removeBadTxByRegion(windowsOne, reads5SSUAll, upstream = 250, downstream = 300, median_multiplier = 20, min_cutoff = 200)
hitMapTIS <- metaWindow(readsSSUAll, windowsOnetest, scoring = "sumPos", withFrames = FALSE,
                        fraction = "SSU", feature = "TIS", zeroPosition = 200)
windowCoveragePlot(hitMapTIS, scoring = "sum", type = "relation to start codon (nt)",
                   output = paste0(heatMapsFolder, "final/leaderLoss_TIS.pdf"))
gg <- windowCoveragePlot(rbindlist(list(hitMapTSS, hitMapTIS)), scoring = "sum", type = "relation to start codon (nt)", setMinToZero = TRUE)
gg <- gg + xlab("") + theme_classic() + theme(panel.spacing = unit(2, "lines")) + theme(strip.background = element_blank(), strip.text.x = element_blank(), title = element_blank(), legend.position = "none", strip.text = element_blank()) + ylab("Sum") +  facet_grid(fraction ~ feature, scales = "free_x",space = "free")
ggsave(gg, filename = paste0(heatMapsFolder, "final/leaderLoss_TISandTSS.pdf"), dpi = 300, width = 300, height = 75, units = "mm")
ggslackR(width = 300, height = 75)


