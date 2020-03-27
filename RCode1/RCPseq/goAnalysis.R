#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO (GO analysis) (This is not used in article)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# If just rerun, start on plots section!!!!!!!!!!!!!!
# This first part depends on a lot of packages, if you only want data before ploting, go step 3

# NOTES:
# This GO analysis uses Gorilla
# It uses adams filtered list of genes with IR scores, using RFP.
# Later a go score of LSU was used, since we started to believe this was made from a bias in
# RCP seq protocol.

### LOAD THIS FIRST PART OF SOURCES FIRST!
library(ORFikPipeline)
library(ggpubr)
plotFolder <- "/export/valenfs/projects/Hakon/AdamVienna/plots/new_plots/"
setwd("/export/valenfs/projects/Hakon/AdamVienna/")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# DATA
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Anotation
gtfPath <- p("/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadTxdb(gtfPath)


# 2: Add go terms column
dt <- setDT(readRDS("kozakTx.rds")) # <- instead you can load this file
orfGo <- getORFsGoTerms(dt$X.gene_id, organism = "Danio rerio")
dt$go <- orfGo
dt[stage == "64cell", "stage"] <- "64"
levels(dt$stage)[1] <- "64"
dt$transcript_id <- as.character(dt$transcript_id)

saveRDS(object = dt, file = "kozakTxWithGo.rds")
dt <- setDT(readRDS("kozakTxWithGo.rds")) # <- final adams

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Expression
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

txNames <- filterTranscripts(txdb,minFiveUTR =  100, 150, 0)
leaders <- ORFik:::loadRegion(txdb, "leader")[txNames]
tx <- ORFik:::loadRegion(txdb)[txNames]
cds <- ORFik:::loadRegion(txdb, "cds")[txNames]
trailers <- ORFik:::loadRegion(txdb, "trailer")


df <- getTCPdfAll()

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Make my matrix
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

dd <- data.table()
for ( s in df$stage) {
  df_stage <- df[df$stage == s,]

  rna <- readGAlignments(df_stage$RNA); seqlevelsStyle(rna) <- seqlevelsStyle(tx)[1]; strand(rna) <- "*"
  rfp <- readGAlignments(df_stage$RFP); seqlevelsStyle(rfp) <- seqlevelsStyle(tx)[1]
  ssu <- readGAlignments(df_stage$SSU); seqlevelsStyle(ssu) <- seqlevelsStyle(tx)[1]

  rfpFPKM <- fpkm(cds, rfp)
  ssuFPKM <- fpkm(leaders, ssu)
  rnaTx <- fpkm(tx, rna)
  rnaCDS <- fpkm(cds, rna)
  dd <- rbindlist(list(dd, data.table(txNames = names(cds), rfp = rfpFPKM, ssu = ssuFPKM, rnaTx, rnaCDS = rnaCDS, stage = s)))
}

dd$stage <- as.factor(dd$stage)
saveRDS(object = dd, file = "expression_me.rds")
dd <- setDT(readRDS("expression_me.rds"))

# make combined
# dd <- dd[rfp > 0. & rnaTx > 0. & ssu > 0.,]
ddd <- dplyr:::left_join(dt, dd, by = c("transcript_id"="txNames", "stage"="stage"))
View(ddd[, c("transcript_id", "rnaCDS", "CDS_totalRNA_FPKM", "stage")])
View(dt[, c("transcript_id", "CDS_totalRNA_FPKM", "stage")])
saveRDS(object = ddd, file = "expression_both.rds")
ddd <- setDT(readRDS("expression_both.rds"))


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# PLOTS (START HERE IF DONE ->!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
setwd("/export/valenfs/projects/Hakon/AdamVienna/")
dt <- setDT(readRDS("kozakTxWithGo.rds")) # <- final adams
dd <- setDT(readRDS("expression_me.rds"))
ddd <- setDT(readRDS("expression_both.rds"))

############ IR PLOTS ######################
# dd <- dd[rfp > 0. & rnaTx > 0. & ssu > 0.,]
# me
ggplot(dd, aes(rfp / (rnaTx + 1), ssu / (rnaTx + 1))) +
  geom_point(size = 0.5) +
  xlim(0, 2000) + ylim(0, 2000) +
  geom_smooth(method='lm', fill = NA) +
  xlab("TE") + ylab("SE") +
  scale_x_log10() + scale_y_log10() +
  facet_grid( ~ stage, scales = "free") +
  theme_bw() +
  stat_cor(method = "pearson") +
  ggtitle(label = "TE vs SE", subtitle = "TE: Translational Efficiency \nSE: Scanning Efficiency")

# adam
ggplot(dt, aes(complete_CDS_RFP_FPKM / (complete_CDS_RNA_FPKM + 1), complete_leader_SSU_FPKM / (complete_CDS_RNA_FPKM  + 1))) +
  geom_point(size = 0.5) +
  xlim(0, 2000) + ylim(0, 2000) +
  geom_smooth(method='lm', fill = NA) +
  xlab("TE") + ylab("SE") +
  scale_x_log10() + scale_y_log10() +
  facet_grid( ~ stage, scales = "free") +
  theme_bw() +
  stat_cor(method = "pearson") +
  ggtitle(label = "TE vs SE (Adam's)", subtitle = "TE: Translational Efficiency \nSE: Scanning Efficiency")




plot(log10(ddd$CDS_totalRNA_FPKM), log10(ddd$rnaCDS))
ggplot(ddd, aes(ssu / (rnaTx + 1), color = TSS_start_sequence)) +
  stat_ecdf() +
  scale_x_log10()

ggplot(ddd, aes(rfp, rnaTx))  +
  geom_point()+
  scale_x_log10() + scale_y_log10() +
  facet_grid( ~ stage, scales = "free") +
  stat_cor(method = "pearson")

# changed by CAGE
tt <- ddd$transcript_id[ddd$cage_tags_at_leader_peak > 5]
tt <- tt[tt %in% names(leaders)]
ggplot(ddd, aes(CDS_totalRNA_FPKM, rnaCDS))  +
  geom_point()+
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method='lm', fill = NA) +
  facet_grid( ~ stage, scales = "free") +
  stat_cor(method = "pearson")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Gorilla
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# why are the high IR high  and the low IR low? (from GOrilla)
# All
fwrite(data.table(dt[order(Initiation_rate_RPF, decreasing = TRUE),]$X.gene_id), file = "/export/valenfs/projects/Hakon/AdamVienna/GO/ir_2019_orderered_RFP.csv")
fwrite(data.table(dt[order(Initiation_rate, decreasing = TRUE),]$X.gene_id), file = "/export/valenfs/projects/Hakon/AdamVienna/GO/ir_2019_orderered_LSU.csv")
fwrite(data.table(dt[order(Initiation_rate, decreasing = FALSE),]$X.gene_id), file = "/export/valenfs/projects/Hakon/AdamVienna/GO/ir_2019_orderered_LSU_increasing.csv")

fwrite(data.table(dt[order(Initiation_rate_RPF, decreasing = TRUE),][stage == "shield",]$X.gene_id), file = "/export/valenfs/projects/Hakon/AdamVienna/GO/i_2019_shield_orderered_RFP.csv")
fwrite(data.table(dt[order(Initiation_rate, decreasing = TRUE),][stage == "shield",]$X.gene_id), file = "/export/valenfs/projects/Hakon/AdamVienna/GO/ir_2019_shield_orderered_LSU.csv")

fwrite(data.table(dt[order(Initiation_rate_RPF, decreasing = TRUE),][stage == "64cell",]$X.gene_id), file = "/export/valenfs/projects/Hakon/AdamVienna/GO/i_2019_64cell_orderered_RFP.csv")
fwrite(data.table(dt[order(Initiation_rate, decreasing = TRUE),][stage == "64cell",]$X.gene_id), file = "/export/valenfs/projects/Hakon/AdamVienna/GO/ir_2019_64cell_orderered_LSU.csv")
# high / low

#' Get output genes from Gorilla
#'
#' How this is made:
#' 1. Find some genes and rank them by a score
#' 2. Send as single list to Gorilla with correct species
#' 3. Pick GO group you want
#' 4. List genes in that group and copy them to a csv file
#' 5. Load that csv file with this function
#' @return a list of transcript names
getGOrilla <- function(fileName, txdb) {
  d <- read.csv(fileName, header = FALSE, skip = 0, sep = "-")
  d <- d$V1
  d <- as.character(d)
  d <- sub(pattern = " ", x = d, replacement = "")

  library(clusterProfiler)
  library(org.Dr.eg.db)
  geneHits <- bitr(d, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Dr.eg.db)

  genes <- geneHits$ENSEMBL

  g <- GenomicFeatures::transcriptLengths(txdb)
  gg <- g[g$gene_id %in% genes,]
  return(gg$tx_name)
}

gg_high <- getGOrilla("./GO/go_high.csv", txdb) # ER
gg_low <- getGOrilla("./GO/go_low.csv", txdb) # Nucleus


#gg_extraC <- dt$transcript_id[dt$go == "extracellular region" | dt$go == "extracellular space"]

#ddG <-  dd[txNames %in% c(gg_high, gg_low, gg_extraC),]
ddG <-  dd[txNames %in% c(gg_high, gg_low),]
ddG[, rfpTPM := (rfp / sum(rfp))*1e6]
ddG[, ssuTPM := (ssu / sum(ssu))*1e6]
ddG$fraction <- "GO_nuc"
ddG$fraction[ddG$txNames %in% c(gg_high)] <- "GO_ER"
#ddG$fraction[ddG$txNames %in% c(gg_extraC)] <- "GO_ER"
ddG$fraction <- as.factor(ddG$fraction)
ddG <- ddG[rfp > 0. & ssu > 0.,]
saveRDS(object = ddG, file = "expression_both_withGorilla_filteredRFPSSU.rds")

ddG$hasUORFs <- ddG$txNames %in% c(txNames(uORF_high), txNames(uORF_low))
ggplot(ddG, aes(rfpTPM, ssuTPM, color = fraction)) +
  geom_point() +
  xlim(0, 3000) + ylim(0, 3000) +
  geom_smooth(method='lm', fill = NA) +
  labs(x = "tpm RFP", y = "tpm SSU") + stat_cor(method = "pearson") +
  facet_wrap( ~ hasUORFs)

ggplot(ddG, aes(rfp / (rnaCDS + 1), ssu / (rnaCDS + 1), color = fraction)) +
  geom_point() +
  geom_smooth(method='lm', fill = NA) +
  scale_x_log10() + scale_y_log10()


# Load whatever SSU; rfp and rna you want. Here I use shield
ddG <- ddG[stage == "shield",]
tNames <- ddG$txNames[ddG$fraction == "GO_ER"]
tNames <- tNames[tNames %in% names(trailers)]
rna <- ORFik:::readBam(df$RNA[3], tx); strand(rna) <- "*"
rfp <- ORFik:::readBam(df$RFP[3], tx)
#ssu <- ORFik:::readBam(df_stage$SSU[3], tx)
ssu <- readsSSUAll

high_ssu <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = ssu, fraction = "GO_ER_ssu")
high_rfp <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rfp, fraction = "GO_ER_rfp")
high_rna <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rna, fraction = "GO_ER_rna")

tNames <- ddG$txNames[ddG$fraction == "GO_nuc"]
tNames <- tNames[tNames %in% names(trailers)]
low_ssu <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = ssu, fraction = "GO_nuc_ssu")
low_rfp <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rfp, fraction = "GO_nuc_rfp")
low_rna <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rna, fraction = "GO_nuc_rna")

# tNames <- ddG$txNames[ddG$fraction == "GO_ER"]
# tNames <- tNames[tNames %in% names(trailers)]
# er_ssu <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = ssu, fraction = "GO_low_ssu")
# er_rfp <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rfp, fraction = "GO_low_rfp")
# er_rna <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rna, fraction = "GO_low_rna")

mergedHiLo <- rbindlist(list(high_ssu, high_rfp, low_ssu, low_rfp, er_ssu, er_rfp))

p <- windowCoveragePlot(dt, scoring = "sum", title = "GO - IR expression normalized by RNA", )

# with RNA
dt[grep(x = dt$fraction, pattern = "high")]$score <- dt[grep(x = dt$fraction, pattern = "high")]$score  / (high_rna$score + 1)
dt[grep(x = dt$fraction, pattern = "low")]$score <- dt[grep(x = dt$fraction, pattern = "low")]$score  / (low_rna$score + 1)
pp <- windowCoveragePlot(dt, scoring = "sum", title = "GO - IR expression normalized by RNA")

# TIS
tNames <- dd$txNames[dd$fraction == "GO_high"]
window <- ORFik:::startRegion(cds[tNames], tx[tNames], upstream = 75, downstream = 74)
window <- window[widthPerGroup(window) == 150]
high_ssu <- metaWindow(x = ssu, windows = window, fraction = "GO_high_ssu", scoring = NULL, scaleTo = 150, feature = "TIS")
high_rfp <- metaWindow(x = rfp, windows = window,  fraction = "GO_high_rfp", scoring = NULL, scaleTo = 150, feature = "TIS")

tNames <- dd$txNames[dd$fraction == "GO_low"]
window <- ORFik:::startRegion(cds[tNames], tx[tNames], upstream = 75, downstream = 74)
window <- window[widthPerGroup(window) == 150]
low_ssu <- metaWindow(x = ssu, windows = window, fraction = "GO_low_ssu", scoring = NULL, scaleTo = 150, feature = "TIS")
low_rfp <- metaWindow(x = rfp, windows = window, fraction = "GO_low_rfp", scoring = NULL, scaleTo = 150, feature = "TIS")
dt <- rbindlist(list(high_ssu, high_rfp, low_ssu, low_rfp))
p <- windowCoveragePlot(dt, scoring = "sum", title = "GO - IR expression", type = "TIS region")

# TSS
tNames <- dd$txNames[dd$fraction == "GO_high"]
window <- ORFik:::startRegion(leaders[tNames], tx[tNames], upstream = 0, downstream = 99)
window <- window[widthPerGroup(window) == 100]
high_ssu <- metaWindow(x = ssu, windows = window, fraction = "GO_high_ssu", zeroPosition = 1, scoring = NULL, feature = "TSS")
high_rfp <- metaWindow(x = rfp, windows = window,  fraction = "GO_high_rfp", zeroPosition = 1, scoring = NULL, feature = "TSS")

tNames <- dd$txNames[dd$fraction == "GO_low"]
window <- ORFik:::startRegion(leaders[tNames], tx[tNames], upstream = 0, downstream = 99)
window <- window[widthPerGroup(window) == 100]
low_ssu <- metaWindow(x = ssu, windows = window, fraction = "GO_low_ssu", zeroPosition = 1, scoring = NULL, feature = "TSS")
low_rfp <- metaWindow(x = rfp, windows = window, fraction = "GO_low_rfp", zeroPosition = 1, scoring = NULL, feature = "TSS")
dt <- rbindlist(list(high_ssu, high_rfp, low_ssu, low_rfp))
p <- windowCoveragePlot(dt, scoring = "sum", title = "GO - IR expression", type = "TSS region")

# Try to recreate dotPlot
gorillaToDotplot <- function(goDf1, outName=NULL) {
  library(ggplot2)
  scaler <- max(nchar(as.character(goDf1$Description)))/39
  # Plot
  p <- ggplot(goDf1, aes(x=GeneRatio, y=Description, color = -log10(qvalue), size = Count)) +
    geom_point(stat='identity')  +
    theme_bw(base_size=13) +
    scale_color_gradientn(colours = c("navy", "purple","red", "red2"), name = "-log10(q-value)") +
    labs(title="GO enrichment",
         subtitle="Onthology: Component") +
    xlim(0.0, 0.3)


  if (!is.null(outName)) {

    ggsave(filename = outName, plot = p, width = 160, height = 140, units = "mm",
           dpi = 300, limitsize = FALSE)
  }
  return(p)
}

# ER (values copied from Gorilla output)
goDf <- data.table(t(matrix(c(3.64E-14, 	2.96E-11, 634,
                 2.18E-13, 	8.84E-11, 544,
                 2.34E-13, 	6.33E-11, 542,
                 9.83E-10, 	1.99E-7, 713,
                 3.1E-7, 5.03E-5, 98,
                 5.48E-6, 7.42E-4, 48,
                 8.13E-6, 9.43E-4, 77,
                 6.37E-5, 6.47E-3, 649,
                 5.45E-4, 4.92E-2, 6),
               nrow = 3)))
colnames(goDf) <- c("pvalue", 	"qvalue", "Count")
goDf$GeneRatio <- goDf$Count/2489
goDf$ID <- c("GO:0044425", "GO:0031224", "GO:0016021", "GO:0016020", "GO:0005783", "GO:0005789", "GO:0044432", "GO:0044444", "GO:0005793")
goDf$Description <- c("membrane part", "intrinsic component of membrane", "integral component of membrane",
                 "membrane", "endoplasmic reticulum", "endoplasmic reticulum membrane", "endoplasmic reticulum part",
                 "cytoplasmic part", "endoplasmic reticulum-Golgi compartment")

goDf$Description <- factor(goDf$Description, levels = goDf$Description[order(goDf$GeneRatio, decreasing = FALSE)])
goDf <- goDf[order(GeneRatio, decreasing = FALSE),]
a <- gorillaToDotplot(goDf, paste0(plotFolder, "/heatmaps/final/GO_ER.pdf"))
# Nuc
goDfNuc <- data.table(t(matrix(c("GO:0005634", 	"nucleus", 	1.24E-11,	1.01E-8, 	669,
                                 "GO:0044428", 	"nuclear part", 	3.37E-6, 	1.37E-3, 	302,
                                 "GO:0044451", 	"nucleoplasm part", 	2.16E-4, 	5.86E-2, 	121,
                                 "GO:0005667", 	"transcription factor complex", 	4.38E-4, 	8.9E-2, 50),
                               nrow = 5)))
colnames(goDfNuc) <- c("ID", "Description", "pvalue","qvalue", "Count")
class(goDfNuc$pvalue) <- "numeric";class(goDfNuc$qvalue) <- "numeric";class(goDfNuc$Count) <- "numeric"
goDfNuc$GeneRatio <- goDfNuc$Count/2489
goDfNuc$Description <- factor(goDfNuc$Description, levels = goDfNuc$Description[order(goDfNuc$GeneRatio, decreasing = FALSE)])
goDfNuc <- goDfNuc[order(GeneRatio, decreasing = FALSE),]
b <- gorillaToDotplot(goDfNuc, paste0(plotFolder, "/heatmaps/final/GO_nuc.pdf"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Length and uORF analysis
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Length
summary(widthPerGroup(leaders[dd$txNames[ddG$fraction == "GO_ER"]]))
summary(widthPerGroup(leaders[dd$txNames[ddG$fraction == "GO_nuc"]]))
t.test(widthPerGroup(leaders[ddG$txNames[ddG$fraction == "GO_nuc"]]), widthPerGroup(leaders[ddG$txNames[ddG$fraction == "GO_ER"]]))
wilcox.test(widthPerGroup(leaders[ddG$txNames[ddG$fraction == "GO_ER"]]), widthPerGroup(leaders[!(names(leaders) %in% ddG$txNames)]))
t.test(widthPerGroup(leaders[ddG$txNames[ddG$fraction == "GO_ER"]]), widthPerGroup(leaders[!(names(leaders) %in% ddG$txNames)]))
# uORFs
getFasta("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa")
cds <- cdsBy(txdb, use.names = TRUE)[names(leaders)]

uorfRegion <- ORFik:::addCdsOnLeaderEnds(leaders, cds)
seqs1 <- ORFik:::txSeqsFromFa(uorfRegion, fa, TRUE)
orfs <- ORFik::findMapORFs(uorfRegion, seqs1, startCodon = "ATG", groupByTx = F)
ret <- filterORFs(orfs)

uORF_high <- ret[which(txNames(ret) %in% gg_high)]
uORF_low <- ret[which(txNames(ret) %in% gg_low)]

te_High <- ORFik:::translationalEff(uORF_high, RNA = rna, RFP = rfp, tx = tx, pseudoCount = TRUE)
te_Low <- ORFik:::translationalEff(uORF_low, RNA = rna, RFP = rfp, tx = tx, pseudoCount = TRUE)


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# 4Ei vs WT
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#<- load from threadingLadders
readsLSU4Ei <- ORFik:::readBam(df2$LSU[df2$type == "4Ei"], tx)
readsSSU4Ei <- ORFik:::readBam(pathSSU4Ei, tx)
readsLSUGoodTrans <- readsLSUGood[readWidths(readsLSUGood) < 35 & readWidths(readsLSUGood) > 25]
readsLSU4EiTrans <- readsLSU4Ei[readWidths(readsLSU4Ei) < 35 & readWidths(readsLSU4Ei) > 25]
rna <- ORFik:::readBam(df$RNA[1], leadersShield); strand(rna) <- "*"

se4Ei <- translationalEff(leadersShield, rna, readsSSU4Ei, cdsShield, pseudoCount = 1)
te4Ei <- translationalEff(cdsShield, rna, readsLSU4EiTrans, cdsShield, pseudoCount = 1)
seWT <- translationalEff(leadersShield, rna, readsSSUGood, cdsShield, pseudoCount = 1)
teWT <- translationalEff(cdsShield, rna, readsLSUGoodTrans, cdsShield, pseudoCount = 1)
#teWTAll <- translationalEff(cdsShield, rna, readsLSUGoodTrans, cdsShield, pseudoCount = 1, with.fpkm = T)
dt <- data.table(se4Ei, te4Ei, seWT, teWT)
dtMelt <- melt(dt)
dtMelt$value <- log10(dtMelt$value)
dtMelt <- dtMelt[dtMelt$value > 0.1 | dtMelt$value < -0.5,]

g <- ggplot(dtMelt, aes(value))
g + geom_density(aes(fill=factor(variable)), alpha=0.2) +
  labs(title="Density plot",
       subtitle="rates of scanning and translation",
       caption="Source: yamilla",
       x="log10 (TE or SE)",
       fill="# type:")

ggplot(dt, aes(se4Ei, seWT)) + geom_point() + geom_smooth(method='lm', fill = NA) + xlim(0, 3) + ylim(0, 3) + stat_cor(method = "pearson") + labs(title = "SE (4Ei vs WT)")
ggplot(dt, aes(te4Ei, teWT)) + geom_point() + geom_smooth(method='lm', fill = NA) + xlim(0, 3) + ylim(0, 3) + labs(title = "TE (4Ei vs WT)") + stat_cor(method = "pearson")
ggplot(dtMelt, aes(variable, value, color = variable)) + geom_violin() + labs(title = "SE (4Ei vs WT)")

# Prepare for GOrilla (make csv and copy genes to browser, not the file!)
ratioSE <- se4Ei / seWT
names(ratioSE) <- names(leadersShield)
oSE <- sort(ratioSE, decreasing = T)
geneNames <- txNamesToGeneNames(names(oSE), txdb)
write.table(data.frame(genes = geneNames, stringsAsFactors = F), file = "4EiSE.csv")
a <- clusterProfiler::enrichGO(gene = geneNames,
                               OrgDb = org.Dr.eg.db,
                               keyType = 'ENSEMBL',
                               ont = "ALL")

ratioTE <- te4Ei / teWT
names(ratioTE) <- names(leadersShield)
oTE <- sort(ratioTE, decreasing = T)
geneNames <- txNamesToGeneNames(names(oTE), txdb)
write.table(data.frame(genes = geneNames, stringsAsFactors = F), file = "4EiTE.csv")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# NEW december 2019, LSU DECREASING IR
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
gg_LSU <- getGOrilla("./GO/go_2019_shield_IR_LSU_decreasing.csv", txdb) # ER
df <- read.experimentl("rcp-experiment_WTvsMZ")
df <- df[c(4,8),]
# New rna seq
dfr <- read.experimentl("Val19")
dfr <- dfr[c(3),]
txdb <- loadTxdb(df)
loadRegions(txdb)
txNames <- filterTranscripts(txdb, 70, 70, 70, longestPerGene = FALSE)
splitRegions(splitList = list(txNames))
gg_LSU <- gg_LSU[gg_LSU %in% txNames]
backGround <- txNames[!(txNames %in% gg_LSU)]
outputLibs(list(df, dfr), seqlevelsStyle(txdb))

region <- startRegion(cds, mrna, T, 50, 49)
splitList <- rep("Others", length(txNames))
splitList[txNames %in% gg_LSU] <- "GO"
splitMetacoverage(region, splitList, df, outdir = "./GO/leader_region_LSU_GO_2019.png", dfr = dfr, title = "TIS metacoverage")
