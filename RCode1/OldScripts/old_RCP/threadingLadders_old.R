#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Threading ladders (SSU, LSU) on TIS / TSS in zebrafish
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
############

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

##### For sanity checks:
df2 <- getTCPNew()
pathSSU4Ei <- "/export/valenfs/data/processed_data/TCP-seq/valen_all_withrRNA/aligned/64_SSU_V12_4Ei.bam"
readsSSU4Ei <- ORFik:::readBam(pathSSU4Ei, tx)
readsLSUV7 <- ORFik:::readBam(df2$LSU[df2$stage == "64"][1], tx)
readsLSUV8 <- ORFik:::readBam(df2$LSU[df2$stage == "64"][2], tx)
readsLSUV5 <- ORFik:::readBam(df2$LSU[df2$stage == "shield"][1], tx)
readsLSUV6 <- ORFik:::readBam(df2$LSU[df2$stage == "shield"][2], tx)
readsLSU <- ORFik:::readBam(df2$LSU[df2$stage == "shield"][3], tx)
readsLSUV7S1 <- ORFik:::readBam(df2$LSU[df2$stage == "sphere"][1], tx)
readsLSUV7S2 <- ORFik:::readBam(df2$LSU[df2$stage == "sphere"][2], tx)
readsLSUV7S3 <- ORFik:::readBam(df2$LSU[df2$stage == "sphere"][3], tx)
readsLSU4Ei <- ORFik:::readBam(df2$LSU[df2$type == "4Ei"], tx)
readsL <- c(readsLSUV7, readsLSU, readsLSUV5)


LSUnt150Path <- "/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_16_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi/merge/LSU_merged_fractions_18_19.bam"
readsLSU150NT <- ORFik:::readBam(LSUnt150Path, leadersShield)
readsLSU150NT5 <- convertToOneBasedRanges(readsLSU150NT, addSizeColumn = T); readsLSU150NT3 <- convertToOneBasedRanges(readsLSU150NT, method = "3prime",  addSizeColumn = T)
readsLSUGood <- c(readsLSU, readsLSUV8, readsLSU150NT)

reads5LGood <- convertToOneBasedRanges(readsLSUGood, method = "5prime", addSizeColumn = TRUE)
reads3LGood <- convertToOneBasedRanges(readsLSUGood, method = "3prime", addSizeColumn = TRUE)

reads5L <- convertToOneBasedRanges(readsL, method = "5prime", addSizeColumn = TRUE)
reads3L <- convertToOneBasedRanges(readsL, method = "3prime", addSizeColumn = TRUE)
reads5LSU <- ORFik:::convertToOneBasedRanges(readsLSU, method = "5prime", addSizeColumn = TRUE)
reads3LSU <- ORFik:::convertToOneBasedRanges(readsLSU, method = "3prime", addSizeColumn = TRUE)


allLSU <- "/export/valenfs/projects/adam/TCP_seq/RCP_files/combined/stages_merged_LSU_updated.bam"
readsLSUAll <- ORFik:::readBam(allLSU, tx)
reads5LSUAll <- ORFik:::convertToOneBasedRanges(readsLSUAll, method = "5prime", addSizeColumn = TRUE)
reads3LSUAll <- ORFik:::convertToOneBasedRanges(readsLSUAll, method = "3prime", addSizeColumn = TRUE)
# Conclusion: 64 cell 7, shield 15 and 8


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TSS
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#### Sanity tests

################################### 150 NT valen 16
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), readsLSU150NT5, upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "150NT_LSU_5prime_TSS.png"), scores = "sum", logIt = F, acLen = 16:75)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), readsLSU150NT3, upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "150NT_LSU_3prime_TSS.png"), scores = "sum", logIt = F, acLen = 16:75)

################################### Best stages
acLen <- 17:68
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads3SSUAll, upstream = 50, downstream = 96, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_3prime_shield_all_rawsum_layout_left.pdf"),
                  scores = "sum", acLen = acLen, legendPos = c(.14, .66))
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads5SSUAll, upstream = 50, downstream = 96, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_5prime_shield_all_rawsum_layout_left.pdf"),
                  scores = "sum", acLen = acLen, legendPos = c(.14, .66))


tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 100), reads5SGood, upstream = 100, downstream = 51, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_5prime_best.png"),
                  scores = "log2sum", acLen = acLen)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads3SGood, upstream = 50, downstream = 101, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_3prime_shield_best.png"),
                  scores = "log2sum", acLen = acLen)
regionBarPlot(leadersShield, extendLeaders(leadersShield, extension = 100), reads5SGood, upstream = 100, downstream = 51, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_5prime_shield_best_bar.png"))
regionBarPlot(leadersShield, extendLeaders(leadersShield, extension = 50), reads3SGood, upstream = 50, downstream = 101, outdir = paste0(heatMapsFolder, "final/TSS_heatmap_3prime_shield_best_bar.png"))


tcpHeatMap_single(leaders, extendLeaders(leaders, extension = 50), reads5LSU, upstream = 30, downstream = 99,
                  outdir = paste0(plotFolder, "/heatmaps/final/TSS_heatmap_5prime_shield_LSU_nonDepleted.png"),
                  scores = "sum")
tcpHeatMap_single(leaders, extendLeaders(leaders, extension = 50), reads3LSU, upstream = 30, downstream = 99,
                  outdir = paste0(plotFolder, "/heatmaps/final/TSS_heatmap_3prime_shield_LSU_nonDepleted.png"),
                  scores = "sum")

df <- getTCPdfAll()
df <- rbind(df, data.frame(SSU = "/export/valenfs/projects/adam/TCP_seq/valen_16/processed_29_03_18/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
                           NA, NA, LSU = "/export/valenfs/projects/adam/TCP_seq/valen_16/processed_29_03_18/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
                           stage = "shield", type = "150"))
df$RNA <- NULL
df$RFP <- NULL

tcpHeatMap_int(leadersShield, extendLeaders(leadersShield, extension = 50), df = df[2,], upstream = 30, downstream = 99,
               outdir = paste0(heatMapsFolder, "TSS_heatmap_"),
               shifting = "5prime", scores = "sum")

tcpHeatMap_int(leadersShield, extendLeaders(leadersShield, extension = 50), df = df[2,], upstream = 30, downstream = 99,
               outdir = paste0(heatMapsFolder, "TSS_heatmap_"),
               shifting = "3prime", scores = "sum")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TIS
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
reads5LFilt <- reads5L[readWidths(reads5L) > 24 & readWidths(reads5L) < 36]
region <- startRegion(cdsShield, txShield, upstream = 73, downstream = 72)
hitMapStart <- metaWindow(reads5L, region, withFrames = TRUE, zeroPosition = 73, feature = "TIS",
                          fraction = "LSU Good")
hitMapStartAll <- metaWindow(readsL, region, withFrames = TRUE, zeroPosition = 73, feature = "TIS",
                             fraction = "LSU Good")
plot <- ORFik:::pSitePlot(hitMapStart, region = "region", type = "TIS CAGE of LSU 5'", length = "all")
ggsave(paste0(heatMapsFolder, "TIS_barplot_LSU_5'.png"), plot, dpi = 300)
plot <- ORFik:::pSitePlot(hitMapStartAll, region = "region", type = "TIS CAGE of LSU whole read",
                          length = "all", scoring = "sum")
ggsave(paste0(heatMapsFolder, "TIS_barplot_LSU_whole.png"), plot, dpi = 300)

tcpHeatMap_single(cdsShield, txShield, reads5LFilt, upstream = 50, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TIS_heatmap_5prime_all_good_LSU_25-35_nonDepleted_zscore.png"),
                  scores = "zscore", logIt = FALSE)
tcpHeatMap_single(cdsShield, txShield, reads5L, upstream = 50, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TIS_heatmap_5prime_all_good_LSU_nonDepleted.png"),
                  scores = "sum")
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads3L, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TIS_heatmap_3prime_all_good_LSU_nonDepleted.png"),
                  scores = "sum")

####################################### FRACTION #########################################
filePath <- "/export/valenfs/data/processed_data/TCP-seq/valen_2018_zebrafish_5_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi"
bamFiles <- grep(".bam",list.files(filePath), value = T)
bamFiles <- paste0(filePath, "/",bamFiles)
df <- data.frame(FRACTION = bamFiles, stage = "shield", type = c(10,12,13,14,15,16,17,18,19,20,9), stringsAsFactors = F)

tcpHeatMap_int(leaders, extendLeaders(leaders, extension = 50), df = df, upstream = 30, downstream = 99,
               outdir = paste0(plotFolder, "/heatmaps/final/TSS_heatmap_"),
               shifting = "5prime", scores = "sum")

tcpHeatMap_int(leaders, extendLeaders(leaders, extension = 50), df = df, upstream = 30, downstream = 99,
               outdir = paste0(plotFolder, "/heatmaps/final/TSS_heatmap_"),
               shifting = "3prime", scores = "sum")

# ALL merged fractions
allSSU <- "/export/valenfs/projects/adam/TCP_seq/RCP_files/combined/stages_merged_SSU_updated.bam"
readsSSUAll <- readGAlignments(allSSU); seqlevelsStyle(readsSSUAll) <- seqlevelsStyle(tx)[1]
reads5SSUAll <- ORFik:::convertToOneBasedRanges(readsSSUAll, method = "5prime", addSizeColumn = TRUE)
reads3SSUAll <- ORFik:::convertToOneBasedRanges(readsSSUAll, method = "3prime", addSizeColumn = TRUE)
allLSU <- "/export/valenfs/projects/adam/TCP_seq/RCP_files/combined/stages_merged_LSU_updated.bam"
readsLSUAll <- ORFik:::readBam(allLSU, tx)
reads5LSUAll <- ORFik:::convertToOneBasedRanges(readsLSUAll, method = "5prime", addSizeColumn = TRUE)
reads3LSUAll <- ORFik:::convertToOneBasedRanges(readsLSUAll, method = "3prime", addSizeColumn = TRUE)


tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads5LSUAll, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TSS_heatmap_5prime_shield_LSU_nonDepleted.png"),
                  scores = "sum")
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads3LSUAll, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TSS_heatmap_3prime_shield_LSU_nonDepleted.png"),
                  scores = "sum")
# barplots

region <- startRegion(leadersShield, extendLeaders(leadersShield, extension = 50), upstream = 30, downstream = 99)
hitMapStart <- metaWindow(reads5LSUAll, region, withFrames = TRUE, zeroPosition = 30, feature = "TSS",
                          fraction = "LSU merged")
plot <- ORFik:::pSitePlot(hitMapStart, region = "region", type = "TSS CAGE of LSU 5'", length = "all")
ggsave(paste0(heatMapsFolder, "TSS_barplot_LSU_5'.png"), plot, dpi = 300)



tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads5SSUAll, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TSS_heatmap_5prime_shield_SSU_nonDepleted.png"),
                  scores = "sum")
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads3SSUAll, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TSS_heatmap_3prime_shield_SSU_nonDepleted.png"),
                  scores = "sum")

region <- startRegion(leaders, extendLeaders(leaders, extension = 50), upstream = 30, downstream = 99)
hitMapStart <- metaWindow(reads5SSUAll, region, withFrames = TRUE, zeroPosition = 30, feature = "TSS",
                          fraction = "SSU merged")
plot <- ORFik:::pSitePlot(hitMapStart, region = "region", type = "TSS CAGE of SSU 5'", length = "all")
ggsave(paste0(heatMapsFolder, "TSS_barplot_SSU_5'.png"), plot, dpi = 300)

# Transcript coverage
rfp <- ORFik:::readBam(df$RFP[3], tx)
# V15
dtSSU <- ORFik:::windowPerTranscript(txdbShield, readsSSUAll, fraction = "SSU")
dtLSU <- ORFik:::windowPerTranscript(txdbShield, readsLSUAll, fraction = "LSU")
dtRFP <- ORFik:::windowPerTranscript(txdbShield, rfp, fraction = "RFP")
merged <- rbindlist(list(dtSSU, dtLSU, dtRFP))
merged <- merged[feature != "trailers",]
plot <- windowCoveragePlot(merged, scoring = "sum", title = "Transcript metacoverage CAGE")
ggsave(paste0(heatMapsFolder, "metacoverage_sum.png"), plot, dpi = 300)
plot <- windowCoveragePlot(merged, scoring = "zscore", title = "Transcript metacoverage CAGE")
ggsave(paste0(heatMapsFolder, "metacoverage_zscore.png"), plot, dpi = 300)

dtSSU <- ORFik:::windowPerTranscript(txdbShield, readsSSU, fraction = "SSU")
dtLSUV15 <- ORFik:::windowPerTranscript(txdbShield, readsLSU, fraction = "LSU_V15")
dtLSUV6 <- ORFik:::windowPerTranscript(txdbShield, readsLSUV6, fraction = "LSU_V6")
dtLSUV5 <- ORFik:::windowPerTranscript(txdbShield, readsLSUV5, fraction = "LSU_V5")
dtLSUV7 <- ORFik:::windowPerTranscript(txdbShield, readsLSUV7, fraction = "LSU_V7")
dtLSUV8 <- ORFik:::windowPerTranscript(txdbShield, readsLSUV8, fraction = "LSU_V8")
dtLSUV7_S1 <- ORFik:::windowPerTranscript(txdbShield, readsLSUV7S1, fraction = "LSU_V7S1")
dtLSUV7_S2 <- ORFik:::windowPerTranscript(txdbShield, readsLSUV7S2, fraction = "LSU_V7S2")
dtLSUV7_S3 <- ORFik:::windowPerTranscript(txdbShield, readsLSUV7S3, fraction = "LSU_V7S3")

merged <- rbindlist(list(dtSSU, dtLSUV7, dtLSUV8, dtLSUV15, dtLSUV6, dtLSUV5, dtLSUV7_S1, dtLSUV7_S2, dtLSUV7_S3, dtRFP))
merged <- merged[feature != "trailers",]
cols <- c("skyblue4", rep("orange", 8), "green")
plot <- windowCoveragePlot(merged, scoring = "sum", title = "Transcript metacoverage CAGE", colors = cols)
ggsave(paste0(heatMapsFolder, "metacoverage_sum_my.png"), plot, dpi = 300)
plot <- windowCoveragePlot(merged, scoring = "zscore", title = "Transcript metacoverage CAGE", colors = cols)
ggsave(paste0(heatMapsFolder, "metacoverage_zscore_my.png"), plot, dpi = 300)
plot <- windowCoveragePlot(merged, scoring = "transcriptNormalized", title = "Transcript metacoverage CAGE", colors = cols)
ggsave(paste0(heatMapsFolder, "metacoverage_transcriptNormalized_my.png"), plot, dpi = 300)
# Do all good
dtLSUGood <- ORFik:::windowPerTranscript(txdbShield, readsL, fraction = "LSU_G")
plot <- windowCoveragePlot(dtLSUGood, scoring = "sum", title = "Transcript metacoverage CAGE", colors = rep("orange", 1))

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads5L, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TSS_heatmap_5prime_shield_LSU_nonDepleted.png"),
                  scores = "sum")
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads3L, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "TSS_heatmap_3prime_shield_LSU_nonDepleted.png"),
                  scores = "sum")

# SSU
scoring <- "sum"
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsSSU, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "SSU_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsSSU, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "SSU_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSU, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU15_shield_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSU, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU15_shield_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV6, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU6_shield_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV6, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU6_shield_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV5, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU5_shield_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV5, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU5_shield_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_64_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_64_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV8, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU8_64_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV8, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU8_64_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7S1, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_sphere1_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7S1, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_sphere1_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7S2, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_sphere2_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7S2, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_sphere2_TSS_heatmap_3prime.png"), scores = scoring)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7S3, addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_sphere3_TSS_heatmap_5prime.png"), scores = scoring)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), convertToOneBasedRanges(readsLSUV7S3, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_sphere3_TSS_heatmap_3prime.png"), scores = scoring)

################ INTRONS #########################
introns <- intronsByTranscript(txdb, use.names = TRUE)
ints <- unlistGrl(introns)
int <- ints[width(ints) > 101]
introns <- groupGRangesBy(int, 1:length(int))

regionWindow(introns, NULL, outdir = p(plotFolder, "sanityCheck_introns_"), df = df,
             title = "introns metacoverage", downstream = 101, upstream = 0, scores = "sum")

tcpHeatMap_int(introns, NULL, outdir = p(plotFolder, "sanityCheck_introns_heatmap_"), df = df,
               downstream = 101, upstream = 0, scores = "sum")

df <- getTCPdf()
df$LSU <- NULL
df <- df[2,]
tcpHeatMap_int(introns, NULL, outdir = p(plotFolder, "aaa_sanityCheck_introns_heatmap_"), df = df,
               downstream = 101, upstream = 0, scores = "sum", shifting = "5prime", zeroPosition = 1)

coverage <- tcpHeatMap_int(introns, NULL, outdir = p(plotFolder, "aaa_sanityCheck_introns_heatmap_"), df = df,
                           downstream = 101, upstream = 0, scores = "sum", shifting = "3prime", zeroPosition = 1,
                           returnCoverage = TRUE)
plot <- ORFik:::coverageHeatMap(coverage = coverage, scoring = "sum")
plot <- plot + xlim(1, 100) + theme_bw()
ggsave(paste0(plotFolder, "aaa_sanityCheck_introns_heatmap_3prime_64_SSU.png"), plot = plot, width = 350, height = 180, units = "mm",
       dpi = 300, limitsize = FALSE)

########################## rRNA ###############################

gtfPathOri <- "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81_chr.gtf"
tt <- import(gtfPathOri)
valids <- tt[grep(x = tt$transcript_biotype, pattern = "rRNA")]
if (any(is.na(valids$transcript_id))) stop("some rRNAs have NA tx names")
rRNA <- tx[unique(valids$transcript_id)]
rRNA <- rRNA[widthPerGroup(rRNA) > 99]
tcpHeatMap_single(rRNA, NULL, reads5, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA_heatmap_5prime_shield_SSU.png"))
tcpHeatMap_single(rRNA, NULL, reads3, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA_heatmap_3prime_shield_SSU.png"))
tcpHeatMap_single(rRNA, NULL, reads5LSU, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA_heatmap_5prime_shield_LSU.png"))
tcpHeatMap_single(rRNA, NULL, reads3LSU, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA_heatmap_3prime_shield_LSU.png"))


dt <- ORFik:::metaWindow(reads, rRNA, zeroPosition = 1, fraction = "SSU", feature = "rRNA", forceUniqueEven = F, scoring = "sum",
                         returnAs = "data.table")
windowCoveragePlot(coverage = dt, output = "rRNA_metaCoverage.png", scoring = "sum", title = "rRNA metaplot", type = "rRNA")
# Split 5S and 5.8S
tt <- import("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81_chr.gtf")
valids <- tt[grep(x = tt$transcript_biotype, pattern = "rRNA")]
valids5S <- valids[grep(x = valids$transcript_name, pattern = "5S_rRNA")]
valids5_8S <- valids[grep(x = valids$transcript_name, pattern = "8S_rRNA")]
if (any(is.na(valids$transcript_id))) stop("some rRNAs have NA tx names")
rRNA5S <- tx[unique(valids5S$transcript_id)]
rRNA5S <- rRNA5S[widthPerGroup(rRNA5S) > 99]
rRNA5.8S <- tx[unique(valids5_8S$transcript_id)]
rRNA5.8S <- rRNA5.8S[widthPerGroup(rRNA5.8S) > 99]

tcpHeatMap_single(rRNA5S, NULL, reads5, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA5S_heatmap_5prime_shield_SSU.png"),
                  scores = "sum")
tcpHeatMap_single(rRNA5S, NULL, reads3, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA5S_heatmap_3prime_shield_SSU.png"),
                  scores = "sum")

tcpHeatMap_single(rRNA5.8S, NULL, reads5, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA5.8S_heatmap_5prime_shield_SSU.png"),
                  scores = "sum")
tcpHeatMap_single(rRNA5.8S, NULL, reads3, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA5.8S_heatmap_3prime_shield_SSU.png"),
                  scores = "sum")
# Test repeatmasker
rRNAPath <- p(dataFolder, "/Zebrafish/danrer10.rRNA.gtf")
tr <- import(rRNAPath); seqlevelsStyle(tr) <- seqlevelsStyle(tx)[1]

validsLSU <- tr[grep(x = tr$gene_id, pattern = "LSU-rRNA")]
rRNA28S <- groupGRangesBy(validsLSU, seq.int(length(validsLSU)))
rRNA28S <- rRNA28S[widthPerGroup(rRNA28S) > 99]

tcpHeatMap_single(rRNA28S[1:10], NULL, reads5[1:10000], upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA28S_heatmap_5prime_shield_SSU.png"))
tcpHeatMap_single(rRNA28S, NULL, reads3, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA28S_heatmap_3prime_shield_SSU.png"))

validsSSU <- tr[grep(x = tr$gene_id, pattern = "SSU-rRNA")]
rRNA18S <- groupGRangesBy(validsSSU, seq.int(length(validsSSU)))
rRNA18S <- rRNA18S[widthPerGroup(rRNA18S) > 99]

tcpHeatMap_single(rRNA18S, NULL, reads5, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA18S_heatmap_5prime_shield_SSU.png"),
                  scores = "sum")
tcpHeatMap_single(rRNA18S, NULL, reads3, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/rRNA18S_heatmap_3prime_shield_SSU.png"),
                  scores = "sum")


########################## snoRNA ###############################

valids <- tt[grep(x = tt$transcript_biotype, pattern = "snoRNA")] # <- snoRNAs
if (any(is.na(valids$transcript_id))) stop("some rRNAs have NA tx names")
snoRNA <- tx[unique(valids$transcript_id)]
snoRNA <- snoRNA[widthPerGroup(snoRNA) > 99]

tcpHeatMap_single(snoRNA, NULL, reads5, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/snoRNA_heatmap_5prime_shield_SSU.png"),
                  scores = "sum")
tcpHeatMap_single(snoRNA, NULL, reads3, upstream = 0, downstream = 99, zeroPosition = 1,
                  outdir = paste0(plotFolder, "/heatmaps/snoRNA_heatmap_3prime_shield_SSU.png"),
                  scores = "sum")

#################################### YEAST #######################################################
yeastPath <- "/export/valenfs/data/references/R64_1_1_yeast/Saccharomyces_cerevisiae.R64-1-1.79_with_UTRs.gtf"
txdbYeast <- loadTxdb(yeastPath)
yeastLeaders <- loadRegion(txdbYeast, "leader")
yeastCDS <- loadRegion(txdbYeast, "cds")
yeastTx <- loadRegion(txdbYeast)
# old version
yeastCage <- ORFik:::readBam("/export/valenfs/data/processed_data/CAGE/wery_2015_S_cerevisiae/final_results/aligned_R64_1_1/SRR2048394.bam", yeastLeaders)
yeastCageScore <- convertToOneBasedRanges(yeastCage, addScoreColumn = T)
yeastCage5 <- convertToOneBasedRanges(yeastCage, addSizeColumn = T)
# new CAGE
# yeastCageScore <- fread.bed("/export/valenfs/data/processed_data/CAGE/wery_2015_S_cerevisiae/final_results/aligned_R64_1_1/merged_CAGE_5prime.bed", yeastLeaders)
yeastLeadersCage <- reassignTSSbyCage(fiveUTRs = yeastLeaders, yeastCageScore, filterValue = 10, removeUnused = F, preCleanup = F)
yeastLeadersCage <- yeastLeadersCage[widthPerGroup(yeastLeadersCage) > 69]
yeastCDS <- yeastCDS[widthPerGroup(yeastCDS) > 69]


yeastSSU <- ORFik:::readBam("/export/valenfs/data/processed_data/TCP-seq/archer_2016_yeast_15nt_without_3nt_trim/final_results/aligned_R64_1_1_tidy_tRNA/WT_SSU.bam", yeastLeaders)
yeastSSU5 <- convertToOneBasedRanges(yeastSSU, addSizeColumn = TRUE); yeastSSU3 <- convertToOneBasedRanges(yeastSSU, method = "3prime",addSizeColumn = TRUE)
# filter bad
region <- ORFik:::startRegion(yeastLeadersCage, yeastLeadersCage, upstream = -2, downstream = 14)
counts <- countOverlaps(region, yeastSSU5);summary(counts); sum(counts > 200)
yeastLeadersCage <- yeastLeadersCage[!(counts > 200)]
counts <- fpkm(yeastLeadersCage, yeastSSU5);summary(counts)

tcpHeatMap_single(yeastLeadersCage[counts > 5], extendLeaders(yeastLeadersCage[counts > 5], extension = 50), yeastSSU5, upstream = 30, downstream = 57, outdir = paste0(heatMapsFolder, "yeast_SSU_5prime_TSS_trans.pdf"), scores = "transcriptNormalized", acLen = 16:70)
tcpHeatMap_single(yeastLeadersCage[counts > 5], extendLeaders(yeastLeadersCage[counts > 5], extension = 50), yeastSSU3, upstream = 20, downstream = 67, outdir = paste0(heatMapsFolder, "yeast_SSU_3prime_TSS_trans.pdf"), scores = "transcriptNormalized", acLen = 16:70)
tcpHeatMap_single(yeastLeadersCage[counts > 5], extendLeaders(yeastLeadersCage[counts > 5], extension = 50), yeastSSU5, upstream = 30, downstream = 57, outdir = paste0(heatMapsFolder, "yeast_SSU_5prime_TSS_sum.pdf"), scores = "sum", acLen = 16:70)
tcpHeatMap_single(yeastLeadersCage[counts > 5], extendLeaders(yeastLeadersCage[counts > 5], extension = 50), yeastSSU3, upstream = 20, downstream = 67, outdir = paste0(heatMapsFolder, "yeast_SSU_3prime_TSS_sum.pdf"), scores = "sum", acLen = 16:70)
regionBarPlot(yeastLeadersCage, extendLeaders(yeastLeadersCage, extension = 50), yeastSSU5, upstream = 30, downstream = 57, outdir = paste0(heatMapsFolder, "yeast_SSU_5prime_TSS_bar.pdf"))
regionBarPlot(yeastLeadersCage, extendLeaders(yeastLeadersCage, extension = 50), yeastSSU3, upstream = 20, downstream = 67, outdir = paste0(heatMapsFolder, "yeast_SSU_3prime_TSS_bar.pdf"))

yeastLSU <- ORFik:::readBam("/export/valenfs/data/processed_data/TCP-seq/archer_2016_yeast_15nt_without_3nt_trim/final_results/aligned_R64_1_1_tidy_tRNA/WT_RS.bam", yeastLeaders)
tcpHeatMap_single(yeastLeadersCage, extendLeaders(yeastLeadersCage, extension = 50), convertToOneBasedRanges(yeastLSU, addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "yeast_LSU_5prime_TSS.png"))
tcpHeatMap_single(yeastLeadersCage, extendLeaders(yeastLeadersCage, extension = 50), convertToOneBasedRanges(yeastLSU, method = "3prime",addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "yeast_LSU_3prime_TSS.png"))
#TIS
tcpHeatMap_single(yeastCDS, yeastTx, yeastSSU5, upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "yeast_SSU_5prime_TIS.png"), scores = "transcriptNormalized", logIt = F, acLen = 1:80)
tcpHeatMap_single(yeastCDS, yeastTx, convertToOneBasedRanges(yeastSSU, method = "3prime",addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "yeast_SSU_3prime_TIS.png"), scores = "zscore", logIt = F)
tcpHeatMap_single(yeastCDS, yeastTx, convertToOneBasedRanges(yeastLSU, addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "yeast_SSU_5prime_TIS.png"), scores = "zscore", logIt = F)
tcpHeatMap_single(yeastCDS, yeastTx, convertToOneBasedRanges(yeastLSU, method = "3prime",addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "yeast_SSU_3prime_TIS.png"), scores = "zscore", logIt = F)

#################################### RFP ##########################################################
# From yamila hela test 2
txNames <- filterTranscripts(Gtf, 60, 60, 60)
getLeaders()
fiveUTRs <- fiveUTRs[txNames]
getCDS()
cds <- cds[txNames]
tx <- getTx()
tx <- tx[txNames]
#cageLeaders <- reassignTSSbyCage(fiveUTRs = fiveUTRs, yeastCageScore,  filterValue = 10, removeUnused = F)
rfpHeLa <-  ORFik:::readBam("/export/valenfs/data/processed_data/Ribo-seq/yamilla_test2_2017_HeLa/aligned_GRCh38/80S.bam", fiveUTRs)

#  AND TIS
tcpHeatMap_single(fiveUTRs, extendLeaders(fiveUTRs, extension = 50), convertToOneBasedRanges(rfpHeLa, addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "hela_rfp_5prime_TSS.png"))
tcpHeatMap_single(fiveUTRs, extendLeaders(fiveUTRs, extension = 50), convertToOneBasedRanges(rfpHeLa, method = "3prime",addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "hela_rfp_3prime_TSS.png"))

tcpHeatMap_single(cds, tx, convertToOneBasedRanges(rfpHeLa, addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "hela_rfp_5prime_TIS.png"))
tcpHeatMap_single(cds, tx, convertToOneBasedRanges(rfpHeLa, method = "3prime", addSizeColumn = TRUE), upstream = 30, downstream = 59, outdir = paste0(heatMapsFolder, "hela_rfp_3prime_TIS.png"))

# metaplot leader, cds, trailer
dt <- ORFik:::windowPerTranscript(Gtf, rfpHeLa, windowSize = 60, fraction = "rfp")
windowCoveragePlot(dt, output = paste0(heatMapsFolder, "hela_rfp_metacoverage.png"))
# barplot of TIS
region <- ORFik:::startRegion(cds, tx, TRUE, 50, 49)
hitMap <- metaWindow(rfpHeLa, region, zeroPosition = 50, fraction = "rfp", feature = "TIS", forceUniqueEven = T, scaleTo = 100, withFrames = T)
ORFik:::pSitePlot(hitMap, length = "all", region = "start", output = paste0(heatMapsFolder, "hela_rfp_TIS_barplot.png"))

# Isodio..
rfpIso <- "/export/valenfs/data/processed_data/Ribo-seq/yamilla_2018_Isodiametra/aligned_ISO/Isodiametra_S12_R1_001.bam"
isoPath <- "/export/valenfs/data/references/Isodiametra_pulchra/bowtie_index/"
txdbIso <- loadTxdb(yeastPath)
isoLeaders <- loadRegion(txdbYeast, "leader")
isoLeaders <- yeastLeaders[widthPerGroup(yeastLeaders) > 59]

####################### 150 NT rRNA ###############################
tcpHeatMap_single(rRNA5S, NULL, readsLSU150NT5, upstream = 0, downstream = 99, zeroPosition = 1, outdir = paste0(heatMapsFolder, "150NT_rRNA5S_heatmap_5prime_SSU.png"), scores = "zscore", logIt = FALSE)
tcpHeatMap_single(rRNA5S, NULL, readsLSU150NT3, upstream = 0, downstream = 99, zeroPosition = 1, outdir = paste0(heatMapsFolder, "150NT_rRNA5S_heatmap_3prime_SSU.png"), scores = "zscore", logIt = TRUE)



################################### lncRNA ######################################
lincRNAs <- ORFik:::loadTranscriptType(gtfPathOri, part = "lincRNA")
txBad <- filterTranscripts(gtfPathOri, 0, 1, 0, longestPerGene = FALSE)
tx <- loadRegion(gtfPathOri)
tx <- tx[txBad]
badLincs <- countOverlaps(lincRNAs, tx)
lincRNAs <- lincRNAs[!(badLincs > 0)]
seqlevelsStyle(lincRNAs) <- seqlevelsStyle(leadersShield)[1]
lincRNAs <- lincRNAs[widthPerGroup(lincRNAs) > 150]

tcpHeatMap_single(lincRNAs, extendLeaders(lincRNAs, 30), reads5SGood, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "lincRNA_heatmap_5prime_SSU_merged.png"),
                  scores = "sum", logIt = TRUE)
tcpHeatMap_single(lincRNAs, extendLeaders(lincRNAs, 30), reads5LGood, upstream = 30, downstream = 99,
                  outdir = paste0(heatMapsFolder, "lincRNA_heatmap_5prime_LSU_merged.png"),
                  scores = "sum", logIt = TRUE)

################################## meta coverage ##################################
readsLSUGoodTrans <- readsLSUGood[readWidths(readsLSUGood) < 35 & readWidths(readsLSUGood) > 25]

dtSSUGood <- ORFik:::windowPerTranscript(txdbShield, readsSSUGood, fraction = "SSU")
dtLSUGood <- ORFik:::windowPerTranscript(txdbShield, readsLSUGoodTrans, fraction = "LSU")
merged <- rbindlist(list(dtSSUGood, dtLSUGood))
plot <- windowCoveragePlot(merged, scoring = "zscore", title = "Transcript metacoverage CAGE")
ggsave(paste0(heatMapsFolder, "metaCoverage_cage_zscore.png"), plot)
merged <- merged[feature != "trailers",]
plot <- windowCoveragePlot(merged, scoring = "log2sum", title = "Transcript metacoverage CAGE")
ggsave(paste0(heatMapsFolder, "metaCoverage_cage_sum.png"), plot)


dtSSUHeat <- tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads5SSUAll, upstream = 30, downstream = 74, outdir = paste0(heatMapsFolder, "TSS_SSU_merged_5prime.png"), scores = "log2sum", logIt = F, acLen = 16:70)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), readsLSUGood, upstream = 30, downstream = 74, outdir = paste0(heatMapsFolder, "TSS_LSU_merged_5prime.png"), scores = "zscore", logIt = F, acLen = 1:70)

adamsFormatHeatMap <- function(plot) {
  p <- plot + scale_y_continuous(breaks = c(50, 100, 150))
  return(p)
}

################################## Silvesterol ##################################
dfSiv <- data.frame(SSU = "/export/valenfs/data/processed_data/TCP-seq/valen_2018_zebrafish_10_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi/merge/SSU_fractions_12_13_14_silvestrol.bam",
                    LSU = "/export/valenfs/data/processed_data/TCP-seq/valen_2018_zebrafish_10_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi/merge/LSU_fractions_18_19_silvestrol.bam",
                    stage = "64", type = "silvesterol")
readsSivSSU <- ORFik:::readBam(dfSiv$SSU, leadersShield)
readsSivLSU <- ORFik:::readBam(dfSiv$LSU, leadersShield)

readsSivSSU5 <- convertToOneBasedRanges(readsSivSSU, addSizeColumn = TRUE)
readsSivSSU3 <- convertToOneBasedRanges(readsSivSSU, method = "3prime", addSizeColumn = TRUE)
readsSivLSU5 <- convertToOneBasedRanges(readsSivLSU, addSizeColumn = TRUE)
readsSivLSU3 <- convertToOneBasedRanges(readsSivLSU, method = "3prime", addSizeColumn = TRUE)

tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), readsSivSSU5, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "SSU7_silvesterol_TSS_heatmap_5prime.png"), scores = "sum", acLen = 19:72)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), readsSivSSU3, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "SSU7_silvesterol_TSS_heatmap_3prime.png"), scores = "sum", acLen = 19:72)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), readsSivLSU5, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_silvesterol_TSS_heatmap_5prime.png"), scores = "sum", acLen = 19:72)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), readsSivLSU3, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_silvesterol_TSS_heatmap_3prime.png"), scores = "sum", acLen = 19:72)

tcpHeatMap_single(cdsShield, txShield, readsSivSSU5, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "SSU7_silvesterol_TIS_heatmap_5prime.png"), scores = "sum", acLen = 19:72)
tcpHeatMap_single(cdsShield, txShield, readsSivSSU3, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "SSU7_silvesterol_TIS_heatmap_3prime.png"), scores = "sum", acLen = 19:72)
tcpHeatMap_single(cdsShield, txShield, readsSivLSU5, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_silvesterol_TIS_heatmap_5prime.png"), scores = "sum", acLen = 19:72)
tcpHeatMap_single(cdsShield, txShield, readsSivLSU3, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "LSU7_silvesterol_TIS_heatmap_3prime.png"), scores = "sum", acLen = 19:72)

# Check WT
reads64SSU <- ORFik:::readBam(df$SSU[1], leadersShield)
reads64SSU5 <- convertToOneBasedRanges(reads64SSU, addSizeColumn = TRUE)
reads64LSU <- ORFik:::readBam("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_LSU_reps_1_2_peaks_removed_translating_filter.bam", leadersShield)
tcpHeatMap_single(leadersShield, extendLeaders(leadersShield, extension = 50), reads64SSU5, upstream = 30, downstream = 99, outdir = paste0(heatMapsFolder, "SSU64_WT_TSS_heatmap_5prime.png"), scores = "sum", acLen = 19:72)

# Coverage
dtSSU <- ORFik:::splitIn3Tx(leadersShield, cdsShield, trailersShield, readsSivSSU, fraction = "SSU_silv")
dtLSU <- ORFik:::splitIn3Tx(leadersShield, cdsShield, trailersShield, readsSivLSU, fraction = "LSU_silv")
dt <- rbindlist(list(dtSSU, dtLSU))
dt <- dt[feature != "trailers",]
windowCoveragePlot(dt, output = paste0(heatMapsFolder,"coveragePlot_silvesterol.png"), scoring = "sum", equalMax = T)

# With structure
means <- read.csv("viennaLeaders_50NT.csv", row.names = 1)
high <- means[mean < quantile(mean, 0.5),]
low <- means[mean >= quantile(mean, 0.5),]
highest <- means[mean < quantile(mean, 0.20),]
lowest <- means[mean >= quantile(mean, 0.80),]

txHigh <- high$txNames
txLow <- low$txNames
txHighest <- highest$txNames
txLowest <- lowest$txNames
dt[genes %in% which(names(cdsShield) %in% txHigh),]$fraction <- paste0(dt[genes %in% which(names(cdsShield) %in% txHigh),]$fraction,"_high")
dt[genes %in% which(names(cdsShield) %in% txHighest),]$fraction <- paste0(dt[genes %in% which(names(cdsShield) %in% txHighest),]$fraction,"_max")
dt[genes %in% which(names(cdsShield) %in% txLow),]$fraction <- paste0(dt[genes %in% which(names(cdsShield) %in% txLow),]$fraction,"_low")
dt[genes %in% which(names(cdsShield) %in% txLowest),]$fraction <- paste0(dt[genes %in% which(names(cdsShield) %in% txLowest),]$fraction,"_min")
dt <- dt[!c(fraction %in% c("SSU_silv" , "LSU_silv")), ]
levels(dt$fraction) <- c("SSU_silv_high_max", "SSU_silv_high", "SSU_silv_low", "SSU_silv_low_min", "LSU_silv_high_max", "LSU_silv_high",
                         "LSU_silv_low", "LSU_silv_low_min")

pp <- windowCoveragePlot(dt, scoring = "zscore", scaleEqual = T, colors = c(rep("skyblue4", 4), rep("orange", 4)))
ggsave(paste0(heatMapsFolder,"structure_coveragePlot_silvesterol.png"), height = 15)
# WT
dtSSUWT <- ORFik:::splitIn3Tx(leadersShield, cdsShield, trailersShield, reads64SSU, fraction = "SSU")
dtLSUWT <- ORFik:::splitIn3Tx(leadersShield, cdsShield, trailersShield, reads64LSU, fraction = "LSU")
dtWT <- rbindlist(list(dtSSUWT, dtLSUWT))
dtWT <- dtWT[feature != "trailers",]

dtWT[genes %in% which(names(cdsShield) %in% txHighest),]$fraction <- paste0(dtWT[genes %in% which(names(cdsShield) %in% txHighest),]$fraction,"_max")
dtWT[genes %in% which(names(cdsShield) %in% txHigh),]$fraction <- paste0(dtWT[genes %in% which(names(cdsShield) %in% txHigh),]$fraction,"_high")
dtWT[genes %in% which(names(cdsShield) %in% txLow),]$fraction <- paste0(dtWT[genes %in% which(names(cdsShield) %in% txLow),]$fraction,"_low")
dtWT[genes %in% which(names(cdsShield) %in% txLowest),]$fraction <- paste0(dtWT[genes %in% which(names(cdsShield) %in% txLowest),]$fraction,"_min")
dtWT <- dtWT[!c(fraction %in% c("SSU" , "LSU")), ]
pp <- windowCoveragePlot(dtWT, scoring = "zscore", scaleEqual = T, colors = c(rep("skyblue4", 4), rep("orange", 4)))
ggsave(paste0(heatMapsFolder,"structure_coveragePlot_WT.png"), height = 15)
################################## RNA ##################################

df <- getTCPdfAll()
readsRNA <- ORFik:::readBam(df$RNA[3], leadersShield)
readsRNA5 <- convertToOneBasedRanges(readsRNA, addSizeColumn = TRUE)
readsRNA3 <- convertToOneBasedRanges(readsRNA, method = "3prime", addSizeColumn = TRUE)

windowsOne <- startRegion(leadersShield, extendLeaders(leadersShield, extension = 75), TRUE, upstream = 75, downstream = 74)
hitMap <- metaWindow(readsRNA, windowsOne, scoring = "meanPos", withFrames = FALSE, fraction = "RNA", feature = "shield")
windowCoveragePlot(hitMap, output = paste0(heatMapsFolder,"RNA_TSS_5prime.png"))

seqlevelsStyle(leaders)  <- seqlevelsStyle(leadersShield)[1]
leaders <- leaders[startSites(leaders) > 100 & widthPerGroup(leaders) > 100]
windowsOne <- startRegion(leaders, extendLeaders(leaders, extension = 75), TRUE, upstream = 75, downstream = 74)
hitMap <- metaWindow(readsRNA, windowsOne, scoring = "meanPos", withFrames = FALSE, fraction = "RNA", feature = "shield")
windowCoveragePlot(hitMap, output = paste0(heatMapsFolder,"RNA_TSS_5prime_nonCAGE.png"))

seqlevelsStyle(tx)  <- seqlevelsStyle(leadersShield)[1]
tx <- tx[startSites(tx) > 100 & widthPerGroup(tx) > 100]
#windowsOne <- startRegion(tx, extendLeaders(leaders, extension = 75), TRUE, upstream = 75, downstream = 74)
hitMap <- metaWindow(readsRNA, tx, scoring = NULL, withFrames = FALSE, fraction = "RNA", feature = "shield")


windowCoveragePlot(hitMap, output = paste0(heatMapsFolder,"RNA_TSS_5prime_nonCAGE_tx.png"))

dtSSU <- ORFik:::splitIn3Tx(leadersShield, cdsShield, trailersShield, readsSivSSU, fraction = "SSU_silv")
windowCoveragePlot(dt, output = paste0(heatMapsFolder,"coveragePlot_silvesterol.png"), scoring = "zscore", equalMax = F)


############################ Ribo-seq chew 2013 heatmaps ########################################
gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath);
seqlevelsStyle(txdb) <- "UCSC"
leaders <- loadRegion(txdb, "leaders")
leaderss <- leaders[seqnamesPerGroup(leaders, F) == grep("M",seqlevels(leaders), value = T)]
trailers <- loadRegion(txdb, "trailers")
trailerss <- trailers[seqnamesPerGroup(trailers, F) == grep("M",seqlevels(trailers), value = T)]
a <- readBam("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Shield_trimmed.bam", "UCSC")
b <- readBam("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Dome_trimmed.bam", "UCSC")
d <- readBam("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Bud_trimmed.bam", "UCSC")
e <- readBam("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/256Cell_trimmed.bam", "UCSC")
"/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/"
a_5p <- convertToOneBasedRanges(a, method = "5prime", addSizeColumn = TRUE)
b_5p <- convertToOneBasedRanges(b, method = "5prime", addSizeColumn = TRUE)
d_5p <- convertToOneBasedRanges(d, method = "5prime", addSizeColumn = TRUE)
e_5p <- convertToOneBasedRanges(e, method = "5prime", addSizeColumn = TRUE)
merged <- c(a_5p, b_5p, d_5p, e_5p)
tcpHeatMap_single(cdss, extendLeaders(cdss, 100), merged, upstream = -60, downstream = 164, outdir = "test.png", scores = "zscore",
                  acLen = 25:33)
tcpHeatMap_single(cdss, extendLeaders(cdss, 100), merged, upstream = -60, downstream = 164, outdir = "test.png", scores = "transcriptNormalized",
                  acLen = 25:33)
# 25:30 has periodicity
windowPerReadLength(grl = cdss, tx = cdss, reads = merged,
                    pShifted = FALSE, upstream = -40, downstream = 159, scoring = "periodic")

tcpHeatMap_single(cdss, cdss, merged, upstream = -40, downstream = 159, outdir = p(heatMapsFolder,"test.png"), scores = "sum",
                  acLen = c(25,27,28,29,30))

####################################### new 4Ei sanity checks (2020) ##################################################
#
source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")
# 64
create.experimentl("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_12_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi/merged", exper = "val_4Ei64",txdb = "/export/valenfs/projects/adam/TCP_seq/transcript_GFF3/shield_transcripts.gff3")
df <- read.experimentl("val_4Ei64")[2,]
ORFik:::simpleLibs(df)
remove.experiments(df)
# shield
create.experimentl("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_13_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi/merge/", exper = "val_4EiShi",txdb = "/export/valenfs/projects/adam/TCP_seq/transcript_GFF3/shield_transcripts.gff3")
df <- read.experimentl("val_4EiShi")[2,]
ORFik:::simpleLibs(df)
remove.experiments(df)

rm(list=ls())
df <- read.experimentl("val_4Ei64")[2,]
heatMapRegion(df, "TSS", scores = "log2sum")
remove.experiments(df)
df <- read.experimentl("val_4EiShi")[2,]
heatMapRegion(df, "TSS", scores = "log2sum")
df@txdb <- "/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db"
heatMapRegion(df, "TSS", scores = "log2sum", outdir =  paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/othertxdb/"))
remove.experiments(df)

df <- read.experimentl("Val_withrRNA")[9,]
outputLibs(df, leadersShield, type = "bedo")
SSU_valRRNA5 <- convertToOneBasedRanges(SSU, method = "5prime", addSizeColumn = TRUE)
leadersShield4Ei <- removeBadTxByRegion(leadersShield, SSU_valRRNA5, upstream = 100, downstream = 51, 10, min_cutoff = 10)
leadersShield4Ei <- removeBadTxByRegion(leadersShield, SSU_valRRNA5, upstream = 80, downstream = -60, 10, min_cutoff = 10)
tcpHeatMap_single(leadersShield4Ei, extendLeaders(leadersShield4Ei, extension = 100), SSU, upstream = 100, downstream = 51, outdir = paste0(heatMapsFolder, "final/SSU_4Ei_TSS_heatmap_5prime.pdf"),
                  scores = "log2sum", acLen = acLen, addFracPlot = FALSE, shifting = "5prime")
tcpHeatMap_single(leadersShield4Ei, extendLeaders(leadersShield4Ei, extension = 100), SSU, upstream = 50, downstream = 101, outdir = paste0(heatMapsFolder, "final/SSU_4Ei_TSS_heatmap_3prime.pdf"),
                  scores = "log2sum", acLen = acLen, addFracPlot = FALSE, shifting = "3prime")
regionBarPlotAll(leadersShield4Ei, extendLeaders(leadersShield4Ei, extension = 100), df, upstream = c(100, 50), downstream = c(51, 101), outdir = paste0(heatMapsFolder, "final/"), location = "TSS", format = ".pdf")

