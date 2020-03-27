#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Sanity tests (NOT USED!!!!!!!!!!!!!!)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

############################## CREATE DATA ###################################################
library(ORFikPipeline)
library(dplyr)
library(ggpubr)
library(GGally)
plotFolder <- "/export/valenfs/projects/Hakon/RCP_SEQ/plots/new_IR_plots/"

adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")

rnaFilt <- 10
rnaFiltNew <- 6.14

#combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM_new >=rnaFiltNew, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate_RPF > 0, !is.na(Initiation_rate_RPF), is.finite(Initiation_rate_RPF), leader_length >= 100)
dt <- adams
dt$Scanning_efficiency_new <-   dt$leader_SSU_FPKM / dt$complete_CDS_totalRNA_FPKM_new

corPlot <- ggplot(data = dt, aes(x = Scanning_efficiency_new, y = Scanning_efficiency)) +
  geom_point() +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_classic() +
  facet_grid( ~ stage)
corPlot
ggsave(p(plotFolder, "corPlot_kozak_IR_Ranking.pdf"), corPlot)

ad <- copy(dt)
# OLD SE
scan <- ad[, c("Scanning_efficiency", "stage")]
setDT(scan)
scan[, .N, by = stage]
dt <- data.table("64cell" = scan$Scanning_efficiency[scan$stage == "64cell"],
                 "sphere" = scan$Scanning_efficiency[scan$stage == "sphere"],
                 "shield" = scan$Scanning_efficiency[scan$stage == "shield"])
dt <- log2(dt)
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "old_SE.pdf"), corPlot)
dt <- dt[is.finite(rowSums(dt)),]
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "old_SE_finite.pdf"), corPlot)

# NEW SE
scan <- ad[, c("Scanning_efficiency_new", "stage")]
setDT(scan)
scan[, .N, by = stage]
dt <- data.table("64cell" = scan$Scanning_efficiency_new[scan$stage == "64cell"],
                 "sphere" = scan$Scanning_efficiency_new[scan$stage == "sphere"],
                 "shield" = scan$Scanning_efficiency_new[scan$stage == "shield"])
dt <- log2(dt)
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "new_SE.pdf"), corPlot)
dt <- dt[is.finite(rowSums(dt)),]
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "new_SE_finite.pdf"), corPlot)

############################# TRANSLATIONAL EFFICIENCY #############################################
ad$TE_new <-   ad$CDS_RFP_FPKM / ad$complete_CDS_totalRNA_FPKM_new
# OLD TE
scan <- ad[, c("TE", "stage")]
setDT(scan)
scan[, .N, by = stage]
dt <- data.table("64cell" = scan$TE[scan$stage == "64cell"],
                 "sphere" = scan$TE[scan$stage == "sphere"],
                 "shield" = scan$TE[scan$stage == "shield"])
dt <- log2(dt)
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "old_TE.pdf"), corPlot)
dt <- dt[is.finite(rowSums(dt)),]
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "old_TE_finite.pdf"), corPlot)
# NEW TE
scan <- ad[, c("TE_new", "stage")]
setDT(scan)
scan[, .N, by = stage]
dt <- data.table("64cell" = scan$TE_new[scan$stage == "64cell"],
                 "sphere" = scan$TE_new[scan$stage == "sphere"],
                 "shield" = scan$TE_new[scan$stage == "shield"])
dt <- log2(dt)
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "new_TE.pdf"), corPlot)
dt <- dt[is.finite(rowSums(dt)),]
corPlot <- ggpairs(dt); corPlot
ggsave(p(plotFolder, "new_TE_finite.pdf"), corPlot)
