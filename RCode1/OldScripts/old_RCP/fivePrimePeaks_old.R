library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(cowplot)
library(data.table)
library(ORFikPipeline)
Palette1 <- c('skyblue4', 'orange')


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create TISU definition from Adams tables
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

region.64cell.TISU           <- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/TSS_matrix/64cell_WT_TISU_TSS_100_200_500/TSS_5prime_100.csv", header=T, check.names = F)
region.64cell.4Ei.100nm.TISU <- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/TSS_matrix/64cell_4Ei_100nm_TISU_TSS_100_200_500/TSS_5prime_100.csv", header=T, check.names = F)
region.64cell.4Ei.10um.TISU  <- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/TSS_matrix/64cell_4Ei_10um_TISU_TSS_100_200_500/TSS_5prime_100.csv", header=T, check.names = F)


region.64cell.TISU           <- filter(region.64cell.TISU, gene_overlaps_another_gene == 0, leader_overlap_an_upstream_gene == 0, gene_overlaps_ncRNA == 0, leader_length <=30, tisu_score >=-12) #TISU
region.64cell.4Ei.100nm.TISU <- filter(region.64cell.4Ei.100nm.TISU, gene_overlaps_another_gene == 0, leader_overlap_an_upstream_gene == 0, gene_overlaps_ncRNA == 0, leader_length <=30, tisu_score >=-12) #TISU
region.64cell.4Ei.10um.TISU  <- filter(region.64cell.4Ei.10um.TISU,  gene_overlaps_another_gene == 0, leader_overlap_an_upstream_gene == 0, gene_overlaps_ncRNA == 0, leader_length <=30, tisu_score >=-12) #TISU

region.64cell.TISU$type           <- "WT TISU"
region.64cell.4Ei.100nm.TISU$type <- "4Ei 100nM TISU"
region.64cell.4Ei.10um.TISU$type  <- "4Ei 10μM TISU"


region.64cell.TISU           <- select (region.64cell.TISU,      "gene_id", "type", "fraction", "position_in_region", "count")
region.64cell.4Ei.100nm.TISU <- select (region.64cell.4Ei.100nm.TISU, "gene_id", "type", "fraction", "position_in_region", "count")
region.64cell.4Ei.10um.TISU  <- select (region.64cell.4Ei.10um.TISU,  "gene_id", "type", "fraction", "position_in_region", "count")

region <- rbind(region.64cell.TISU, region.64cell.4Ei.100nm.TISU, region.64cell.4Ei.10um.TISU)

region$fraction <- factor(region$fraction, levels=c("SSU", "LSU", "40S", "80S"))
region$fraction[region$fraction == "LSU" ] <- "80S"
region$fraction[region$fraction == "SSU" ] <- "40S"
region$fraction <- factor(region$fraction, levels=c("40S", "80S"))

region$type <- factor(region$type, levels=c("WT TISU", "4Ei 100nM TISU", "4Ei 10μM TISU"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Use TISU definition to find SE difference in 4Ei vs WT
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# Gene definition
genes <- unique(region$gene_id)
# WT (used for all)
SSUWT <- readBam("/export/valenfs/projects/adam/TCP_seq/RCP_files/combined/stages_merged_SSU_updated.bam")

################ 64 CELL
# Load Transcript table
adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
adams <- adams[stage == "64cell",]
if (!all(genes %in% adams$X.gene_id)) stop("Not all genes found!")
adams <- adams[adams$X.gene_id %in% genes, ];dim(adams)


# Load rna-seq for 10mM 4Ei
df <- read.experimentl("Val19", expInVarName = FALSE)
l <- makeSummarizedExperimentFromBam(saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_leaders.rds")
txNames <- rownames(assay(l))
Ei <- assay(l)[,1]
names(Ei) <- ORFik:::txNamesToGeneNames(names(Ei), df)
Ei <- Ei[names(Ei) %in% genes]
dt <- data.table(RNA_4Ei = Ei, X.gene_id = names(Ei))
adams <- adams[adams$X.gene_id %in% dt$X.gene_id, ]
final <- data.table::merge.data.table(adams, dt, all = FALSE, by = c('X.gene_id'))
final <- final[, .(X.gene_id, complete_CDS_totalRNA_FPKM_new, RNA_4Ei, leader_SSU_FPKM)]
# 4Ei 10 uM
df <- read.experimentl("val_4Ei64", expInVarName = TRUE)[2,]; outputLibs(df)
# 4Ei 0.1 uM
df <- read.experimentl("val_10"); outputLibs(df)


loadRegions(df, "leaders", extension = "4Ei", names.keep = adams$transcript_id)
seqlevelsStyle(leaders4Ei) <- seqlevelsStyle(SSU_4Ei)
SSU0.1 <- fpkm(leaders4Ei, SSU_4Ei); names(SSU0.1) <- NULL
SSU10 <- fpkm(leaders4Ei, val_4Ei64_SSU); names(SSU10) <- NULL
SSUWTFPKM <- fpkm(leaders4Ei, SSUWT); names(SSUWTFPKM) <- NULL
dt2 <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leaders4Ei), df),
                  SSU0.1FPKM = SSU0.1, SSU10FPKM = SSU10, SSUWTFPKM)
final2 <- data.table::merge.data.table(final, dt2, all = FALSE, by = c('X.gene_id'))
res <- data.table(RNA = c(final2$complete_CDS_totalRNA_FPKM_new, final2$RNA_4Ei), SSU = c(final2$SSUWTFPKM, final2$SSU10FPKM),
                  type = c(rep("WT", nrow(final2)), rep("4Ei10mM", nrow(final2))))
res$RNA <- res$RNA + 1
res$SE <- res$SSU / res$RNA

summary(res)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
master <- data.table(res, type2 = paste(res$type, "TISU", "64cell"))

# With SE of WT removed if 0:
res2 <- copy(res)
res2$SE <- res2$SE + 0.000001
dt3 <- copy(dt2)
dt3 <- dt3[dt3$SSUWTFPKM > 0,]
final3 <- data.table::merge.data.table(final, dt3, all = FALSE, by = c('X.gene_id'))
res3 <- data.table(RNA = c(final3$complete_CDS_totalRNA_FPKM_new, final3$RNA_4Ei), SSU = c(final3$SSUWTFPKM, final3$SSU10FPKM),
                  type = c(rep("WT", nrow(final3)), rep("4Ei10mM", nrow(final3))))
res3$RNA <- res3$RNA + 1
res3$SE <- res3$SSU / res3$RNA

t.test(log2(res3[type == "WT",]$SE + 0.000001), log2(res3[type != "WT",]$SE + 0.000001), paired = TRUE, alternative = "two.sided")
wilcox.test(res3[type == "WT",]$SE, res3[type != "WT",]$SE, paired = TRUE)
wilcox.test(res2[type == "WT",]$SE, res2[type != "WT",]$SE, paired = TRUE)

################ SHIELD

adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
adams <- adams[stage == "shield",]
if (!all(genes %in% adams$X.gene_id)) stop("Not all genes found!")
adams <- adams[adams$X.gene_id %in% genes, ];dim(adams)


# Load rna-seq for 10mM 4Ei
l <- makeSummarizedExperimentFromBam(saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_leaders.rds")
txNames <- rownames(assay(l))
Ei <- assay(l)[,4]
names(Ei) <- ORFik:::txNamesToGeneNames(names(Ei), df)
Ei <- Ei[names(Ei) %in% genes]
dt <- data.table(RNA_4Ei = Ei, X.gene_id = names(Ei))
adams <- adams[adams$X.gene_id %in% dt$X.gene_id, ]
final <- data.table::merge.data.table(adams, dt, all = FALSE, by = c('X.gene_id'))
final <- final[, .(X.gene_id, complete_CDS_totalRNA_FPKM_new, RNA_4Ei, leader_SSU_FPKM)]
# 4Ei 10 uM
df <- read.experimentl("val_4EiShi.csv", expInVarName = TRUE)[2,]; outputLibs(df)

loadRegions(df, "leaders", extension = "4Ei", names.keep = adams$transcript_id)
seqlevelsStyle(leaders4Ei) <- seqlevelsStyle(SSU_4Ei)
SSU10 <- fpkm(leaders4Ei, val_4EiShi_SSU); names(SSU10) <- NULL
SSUWTFPKM <- fpkm(leaders4Ei, SSUWT); names(SSUWTFPKM) <- NULL
dt2 <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leaders4Ei), df),
                  SSU10FPKM = SSU10, SSUWTFPKM)
final2 <- data.table::merge.data.table(final, dt2, all = FALSE, by = c('X.gene_id'))
res <- data.table(RNA = c(final2$complete_CDS_totalRNA_FPKM_new, final2$RNA_4Ei), SSU = c(final2$SSUWTFPKM, final2$SSU10FPKM),
                  type = c(rep("WT", nrow(final2)), rep("4Ei10mM", nrow(final2))))
res$RNA <- res$RNA + 1
res$SE <- res$SSU / res$RNA

summary(res)

coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
master <- rbindlist(list(master, data.table(res, type2 = paste(res$type, "TISU", "shield"))))

# With SE of WT removed if 0:
res2 <- copy(res)
res2$SE <- res2$SE + 0.000001
dt3 <- copy(dt2)
dt3 <- dt3[dt3$SSUWTFPKM > 0,]
final3 <- data.table::merge.data.table(final, dt3, all = FALSE, by = c('X.gene_id'))
res3 <- data.table(RNA = c(final3$complete_CDS_totalRNA_FPKM_new, final3$RNA_4Ei), SSU = c(final3$SSUWTFPKM, final3$SSU10FPKM),
                   type = c(rep("WT", nrow(final3)), rep("4Ei10mM", nrow(final3))))
res3$RNA <- res3$RNA + 1
res3$SE <- res3$SSU / res3$RNA

t.test(log2(res3[type == "WT",]$SE + 0.000001), log2(res3[type != "WT",]$SE + 0.000001), paired = TRUE, alternative = "two.sided")
wilcox.test(res3[type == "WT",]$SE, res3[type != "WT",]$SE, paired = TRUE)
wilcox.test(res2[type == "WT",]$SE, res2[type != "WT",]$SE, paired = TRUE)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Use All - TISU = Background definition to find SE difference in 4Ei vs WT
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

############ Shield
adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
adams <- adams[stage == "shield",]
if (!all(genes %in% adams$X.gene_id)) stop("Not all genes found!")
adams <- adams[!(adams$X.gene_id %in% genes), ];dim(adams)
backgenes <- adams$X.gene_id

# Load rna-seq for 10mM 4Ei
l <- makeSummarizedExperimentFromBam(saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_leaders.rds")
txNames <- rownames(assay(l))
Ei <- assay(l)[,4]
names(Ei) <- ORFik:::txNamesToGeneNames(names(Ei), df)
Ei <- Ei[names(Ei) %in% backgenes]
dt <- data.table(RNA_4Ei = Ei, X.gene_id = names(Ei))
adams <- adams[adams$X.gene_id %in% dt$X.gene_id, ]
final <- data.table::merge.data.table(adams, dt, all = FALSE, by = c('X.gene_id'))
final <- final[, .(X.gene_id, complete_CDS_totalRNA_FPKM_new, RNA_4Ei, leader_SSU_FPKM)]
# 4Ei 10 uM
df <- read.experimentl("val_4EiShi.csv", expInVarName = TRUE)[2,]; outputLibs(df)

loadRegions(df, "leaders", extension = "4Ei", names.keep = adams$transcript_id)
seqlevelsStyle(leaders4Ei) <- seqlevelsStyle(SSU_4Ei)
SSU10 <- fpkm(leaders4Ei, val_4EiShi_SSU); names(SSU10) <- NULL
SSUWTFPKM <- fpkm(leaders4Ei, SSUWT); names(SSUWTFPKM) <- NULL
dt2 <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leaders4Ei), df),
                  SSU10FPKM = SSU10, SSUWTFPKM)
final2 <- data.table::merge.data.table(final, dt2, all = FALSE, by = c('X.gene_id'))
res <- data.table(RNA = c(final2$complete_CDS_totalRNA_FPKM_new, final2$RNA_4Ei), SSU = c(final2$SSUWTFPKM, final2$SSU10FPKM),
                  type = c(rep("WT", nrow(final2)), rep("4Ei10mM", nrow(final2))))
res$RNA <- res$RNA + 1
res$SE <- res$SSU / res$RNA

summary(res)

coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
master <- rbindlist(list(master, data.table(res, type2 = paste(res$type, "shield"))))

# With SE of WT removed if 0:
res2 <- copy(res)
res2$SE <- res2$SE + 0.000001
dt3 <- copy(dt2)
dt3 <- dt3[dt3$SSUWTFPKM > 0,]
final3 <- data.table::merge.data.table(final, dt3, all = FALSE, by = c('X.gene_id'))
res3 <- data.table(RNA = c(final3$complete_CDS_totalRNA_FPKM_new, final3$RNA_4Ei), SSU = c(final3$SSUWTFPKM, final3$SSU10FPKM),
                   type = c(rep("WT", nrow(final3)), rep("4Ei10mM", nrow(final3))))
res3$RNA <- res3$RNA + 1
res3$SE <- res3$SSU / res3$RNA

t.test(log2(res3[type == "WT",]$SE + 0.000001), log2(res3[type != "WT",]$SE + 0.000001), paired = TRUE, alternative = "two.sided")
wilcox.test(res3[type == "WT",]$SE, res3[type != "WT",]$SE, paired = TRUE)
wilcox.test(res2[type == "WT",]$SE, res2[type != "WT",]$SE, paired = TRUE)

############ 64 cell
adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
adams <- adams[stage == "64cell",]
if (!all(genes %in% adams$X.gene_id)) stop("Not all genes found!")
adams <- adams[!(adams$X.gene_id %in% genes), ];dim(adams)
backgenes <- adams$X.gene_id

# Load rna-seq for 10mM 4Ei
l <- makeSummarizedExperimentFromBam(saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_leaders.rds")
txNames <- rownames(assay(l))
Ei <- assay(l)[,1]
names(Ei) <- ORFik:::txNamesToGeneNames(names(Ei), df)
Ei <- Ei[names(Ei) %in% backgenes]
dt <- data.table(RNA_4Ei = Ei, X.gene_id = names(Ei))
adams <- adams[adams$X.gene_id %in% dt$X.gene_id, ]
final <- data.table::merge.data.table(adams, dt, all = FALSE, by = c('X.gene_id'))
final <- final[, .(X.gene_id, complete_CDS_totalRNA_FPKM_new, RNA_4Ei, leader_SSU_FPKM)]
# 4Ei 10 uM
df <- read.experimentl("val_4Ei64", expInVarName = TRUE)[2,]; outputLibs(df)

loadRegions(df, "leaders", extension = "4Ei", names.keep = adams$transcript_id)
seqlevelsStyle(leaders4Ei) <- seqlevelsStyle(SSU_4Ei)
SSU10 <- fpkm(leaders4Ei, val_4Ei64_SSU); names(SSU10) <- NULL
SSUWTFPKM <- fpkm(leaders4Ei, SSUWT); names(SSUWTFPKM) <- NULL
dt2 <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leaders4Ei), df),
                  SSU10FPKM = SSU10, SSUWTFPKM)
final2 <- data.table::merge.data.table(final, dt2, all = FALSE, by = c('X.gene_id'))
res <- data.table(RNA = c(final2$complete_CDS_totalRNA_FPKM_new, final2$RNA_4Ei), SSU = c(final2$SSUWTFPKM, final2$SSU10FPKM),
                  type = c(rep("WT", nrow(final2)), rep("4Ei10mM", nrow(final2))))
res$RNA <- res$RNA + 1
res$SE <- res$SSU / res$RNA

summary(res)

coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
master <- rbindlist(list(master, data.table(res, type2 = paste(res$type, "64cell"))))

# With SE of WT removed if 0:
res2 <- copy(res)
res2$SE <- res2$SE + 0.000001
dt3 <- copy(dt2)
dt3 <- dt3[dt3$SSUWTFPKM > 0,]
final3 <- data.table::merge.data.table(final, dt3, all = FALSE, by = c('X.gene_id'))
res3 <- data.table(RNA = c(final3$complete_CDS_totalRNA_FPKM_new, final3$RNA_4Ei), SSU = c(final3$SSUWTFPKM, final3$SSU10FPKM),
                   type = c(rep("WT", nrow(final3)), rep("4Ei10mM", nrow(final3))))
res3$RNA <- res3$RNA + 1
res3$SE <- res3$SSU / res3$RNA

t.test(log2(res3[type == "WT",]$SE + 0.000001), log2(res3[type != "WT",]$SE + 0.000001), paired = TRUE, alternative = "two.sided")
wilcox.test(res3[type == "WT",]$SE, res3[type != "WT",]$SE, paired = TRUE)
wilcox.test(res2[type == "WT",]$SE, res2[type != "WT",]$SE, paired = TRUE)


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# MASTER LIST, all varianters, definition to find SE difference in 4Ei vs WT
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

master64 <- master[grep("64cell", type2),]
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master64), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master64), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()

masterShield <- master[grep("shield", type2),]
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(masterShield), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(masterShield), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()

master$type3 <- gsub("64cell|shield", master$type2, replacement = "")
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master), aes(x = log2(SE), colour= as.factor(type3))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()
master$type3 <- gsub("64cell|shield", master$type2, replacement = "")
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master), aes(x = log10(SE), colour= as.factor(type3))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-2.5, 2.5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()

masterTISU <- master[grep(pattern = "TISU", master$type3),]
masterTISUWT <- masterTISU[type == "WT",]
masterTISUWT64 <- masterTISUWT[grep(pattern = "64cell", type2),]
summary(masterTISUWT64$SE)

masterTISU4Ei <- masterTISU[type != "WT",]
masterTISU4Ei64 <- masterTISU4Ei[grep(pattern = "64cell", type2),]
summary(masterTISU4Ei64$SE)

box <- data.table(SE = c(masterTISUWT64$SE, masterTISU4Ei64$SE),
                  type = c(rep("WT", length(masterTISUWT64$SE)),rep("4Ei", length(masterTISU4Ei64$SE))))
ggplot(box, aes(y = log2(SE), fill = type)) + geom_boxplot()


# PSEUDO count
masterp <- copy(master)
masterp$SE <- masterp$SE + 0.0001
master64 <- masterp[grep("64cell", type2),]
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master64), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master64), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()

masterShieldp <- masterp[grep("shield", type2),]
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(masterShieldp), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()

# Not pseudo
masterShield <- master[grep("shield", type2),]
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(masterShield), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(masterShield), aes(x = log2(SE), colour= as.factor(type2))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()

master$type3 <- gsub("64cell|shield", master$type2, replacement = "")
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master), aes(x = log2(SE), colour= as.factor(type3))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()
master$type3 <- gsub("64cell|shield", master$type2, replacement = "")
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master), aes(x = log10(SE), colour= as.factor(type3))) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-2.5, 2.5))
coverage_plot_transcript_raw_100_ecdf
ggslackR()

masterTISU <- masterp[grep(pattern = "TISU", masterp$type3),]
masterTISUWT <- masterTISU[type == "WT",]
masterTISUWT64 <- masterTISUWT[grep(pattern = "64cell", type2),]
summary(masterTISUWT64$SE)

masterTISU4Ei <- masterTISU[type != "WT",]
masterTISU4Ei64 <- masterTISU4Ei[grep(pattern = "64cell", type2),]
summary(masterTISU4Ei64$SE)

box <- data.table(SE = c(masterTISUWT64$SE, masterTISU4Ei64$SE),
                  type = c(rep("WT", length(masterTISUWT64$SE)),rep("4Ei", length(masterTISU4Ei64$SE))))
ggplot(box, aes(y = log2(SE), fill = type)) + geom_boxplot()
ggslackR()
ggplot(box, aes(y = (SE), fill = type)) + geom_boxplot()
ggslackR()
ggplot(box, aes(y = (SE), fill = type)) + geom_boxplot() + coord_cartesian(ylim = c(0, 10))
ggslackR()

# Change in SE WT vs 4Ei (all tx)

coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(master64), aes(x = log2(SE), color = as.factor(type))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()

coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(masterShieldp), aes(x = log2(SE), color = as.factor(type))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf
ggslackR()


# Other data-set sanity test
