# Go analysis
# This first part depends on a lot of packages, if you only want data before ploting, go step 3
# 1: get adams data (only run if you have all packages)
### NOTE: LOAD THIS FIRST PART OF SOURCES FIRST!
rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/") 
library(ORFik)
source("./DataBaseSetup.R")
dbDisconnect(uorfDB)
rm(uorfDB)
setwd("/export/valenfs/projects/uORFome/RCode1/") 
source("./TempScripts/tcp_pipeline.R")
plotFolder <- "/export/valenfs/projects/uORFome/AdamVienna/plots/new_plots/"
library(ggpubr)
library(ORFik)
setwd(p(mainFolder, "/AdamVienna/"))

gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath);

# df <- combined.dat.ranked.filter # <- this is output from adams file, so can no run this from here.
# Location: 
#dt <- setDT(df)
# txNames <- dt$transcript_id
#saveRDS(df, file = "kozakTx.rds") 


# 2: Add go terms column
dt <- setDT(readRDS("kozakTx.rds")) # <- instead you can load this file
orfGo <- getORFsGoTerms(dt$X.gene_id, organism = "Danio rerio")
dt$go <- orfGo
dt[stage == "64cell", "stage"] <- "64"
levels(dt$stage)[1] <- "64"
dt$transcript_id <- as.character(dt$transcript_id)

saveRDS(object = dt, file = "kozakTxWithGo.rds")
dt <- setDT(readRDS("kozakTxWithGo.rds")) # <- final adams

######################## Expression ###############################

txNames <- filterTranscripts(txdb,minFiveUTR =  100, 150, 0)
leaders <- ORFik:::loadRegion(txdb, "leader")[txNames]
tx <- ORFik:::loadRegion(txdb)[txNames]
cds <- ORFik:::loadRegion(txdb, "cds")[txNames]
trailers <- ORFik:::loadRegion(txdb, "trailer")


df <- getTCPdfAll()

################################ make mine ##################################

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
######################## PLOTS #########################################
######################## STEP 3 ########################################
# START HERE IF DONE ->
setwd(p(mainFolder, "/AdamVienna/"))
dt <- setDT(readRDS("kozakTxWithGo.rds")) # <- final adams
dd <- setDT(readRDS("expression_me.rds"))
ddd <- setDT(readRDS("expression_both.rds"))

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

######################## GORILLA #######################
# why are the high high  and the low low? (from GOrilla) 
# All

# high / low
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

gg_high <- getGOrilla("./go_high.csv", txdb) # ER
gg_low <- getGOrilla("./go_low.csv", txdb) # Nucleus
#gg_extraC <- dt$transcript_id[dt$go == "extracellular region" | dt$go == "extracellular space"]

#ddG <-  dd[txNames %in% c(gg_high, gg_low, gg_extraC),]
ddG <-  dd[txNames %in% c(gg_high, gg_low),]

ddG$fraction <- "GO_nuc"
ddG$fraction[ddG$txNames %in% c(gg_high)] <- "GO_ER"
#ddG$fraction[ddG$txNames %in% c(gg_extraC)] <- "GO_ER"
ddG$fraction <- as.factor(ddG$fraction)
ddG <- ddG[rfp > 0. & ssu > 0.,]

ggplot(ddG, aes((rfp / sum(rfp)), (ssu / sum(ssu)), color = fraction)) + 
  geom_point() + 
  xlim(10, 2000) + ylim(10, 2000) + 
  geom_smooth(method='lm', fill = NA) + 
  scale_x_log10() + scale_y_log10()

ggplot(ddG, aes(rfp / (rnaCDS + 1), ssu / (rnaCDS + 1), color = fraction)) + 
  geom_point() + 
  xlim(10, 2000) + ylim(10, 2000) + 
  geom_smooth(method='lm', fill = NA) + 
  scale_x_log10() + scale_y_log10()

# Transcript
ddG <-  dd[txNames %in% c(gg_high, gg_low, gg_extraC),]

ddG$fraction <- "GO_low"
ddG$fraction[ddG$txNames %in% c(gg_high)] <- "GO_high"
ddG$fraction[ddG$txNames %in% c(gg_extraC)] <- "GO_ER"
ddG$fraction <- as.factor(ddG$fraction)
ddG <- ddG[stage == "shield",]
tNames <- ddG$txNames[ddG$fraction == "GO_high"] 
tNames <- tNames[tNames %in% names(trailers)]
high_ssu <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = ssu, fraction = "GO_high_ssu")
high_rfp <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rfp, fraction = "GO_high_rfp")
high_rna <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rna, fraction = "GO_high_rna")

tNames <- ddG$txNames[ddG$fraction == "GO_low"] 
tNames <- tNames[tNames %in% names(trailers)]
low_ssu <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = ssu, fraction = "GO_low_ssu")
low_rfp <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rfp, fraction = "GO_low_rfp")
low_rna <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rna, fraction = "GO_low_rna")

tNames <- ddG$txNames[ddG$fraction == "GO_ER"] 
tNames <- tNames[tNames %in% names(trailers)]
er_ssu <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = ssu, fraction = "GO_low_ssu")
er_rfp <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rfp, fraction = "GO_low_rfp")
er_rna <- ORFik:::splitIn3Tx(leaders[tNames], cds[tNames], trailers[tNames], reads = rna, fraction = "GO_low_rna")

dt <- rbindlist(list(high_ssu, high_rfp, low_ssu, low_rfp, er_ssu, er_rfp))

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

###################################### Length and uORF analysis ################
# Length
summary(widthPerGroup(leaders[dd$txNames[dd$fraction == "GO_high"]]))
summary(widthPerGroup(leaders[dd$txNames[dd$fraction == "GO_low"]]))
t.test(widthPerGroup(leaders[dd$txNames[dd$fraction == "GO_low"]]), widthPerGroup(leaders[dd$txNames[dd$fraction == "GO_high"]]))

# uORFs

cds <- cdsBy(txdb, use.names = TRUE)[names(leaders)]

uorfRegion <- ORFik:::addCdsOnLeaderEnds(leaders, cds)
seqs1 <- ORFik:::txSeqsFromFa(uorfRegion, fa, TRUE)
orfs <- ORFik::findMapORFs(uorfRegion, seqs1, startCodon = "ATG", groupByTx = F)
ret <- filterORFs(orfs)

uORF_high <- ret[which(txNames(ret) %in% gg_high)]
uORF_low <- ret[which(txNames(ret) %in% gg_low)]

te_High <- ORFik:::translationalEff(uORF_high, RNA = rna, RFP = rfp, tx = tx, pseudoCount = TRUE)
te_Low <- ORFik:::translationalEff(uORF_low, RNA = rna, RFP = rfp, tx = tx, pseudoCount = TRUE)


#################################### 4Ei vs WT ####################################
#<- load from threadingLadders
readsLSU4Ei <- ORFik:::readBam(df2$LSU[df2$type == "4Ei"], tx)
readsSSU4Ei <- ORFik:::readBam(pathSSU4Ei, tx)
readsLSUGoodTrans <- readsLSUGood[readWidths(readsLSUGood) < 35 & readWidths(readsLSUGood) > 25]
readsLSU4EiTrans <- readsLSU4Ei[readWidths(readsLSU4Ei) < 35 & readWidths(readsLSU4Ei) > 25]
rna <- ORFik:::readBam(df$RFP[3], leadersShield)

se4Ei <- translationalEff(leadersShield, rna, readsSSU4Ei, cdsShield, pseudoCount = 1)
te4Ei <- translationalEff(cdsShield, rna, readsLSU4EiTrans, cdsShield, pseudoCount = 1)
seWT <- translationalEff(leadersShield, rna, readsSSUGood, cdsShield, pseudoCount = 1)
teWT <- translationalEff(cdsShield, rna, readsLSUGoodTrans, cdsShield, pseudoCount = 1)
teWTAll <- translationalEff(cdsShield, rna, readsLSUGoodTrans, cdsShield, pseudoCount = 1, with.fpkm = T)
dt <- data.table(se4Ei, te4Ei, seWT, teWT)
dt <- melt(dt)
dt$value <- log10(dt$value)
dt <- dt[dt$value > 0.1 | dt$value < -0.5,]

g <- ggplot(dt, aes(value))
g + geom_density(aes(fill=factor(variable)), alpha=0.2) + 
  labs(title="Density plot", 
       subtitle="rates of scanning and translation",
       caption="Source: yamilla",
       x="log10 (TE or SE)",
       fill="# type:")

ggplot() + geom_point(aes(se4Ei, seWT)) + xlim(0, 3) + ylim(0, 3) + labs(title = "SE (4Ei vs WT)")
ggplot() + geom_point(aes(te4Ei, teWT)) + xlim(0, 3) + ylim(0, 3) + labs(title = "TE (4Ei vs WT)")
ggplot(dt, aes(variable, value, color = variable)) + geom_violin() + labs(title = "SE (4Ei vs WT)")
