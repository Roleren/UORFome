# Kozak heatmap
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Initiation rate by nucleotide and position
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFikPipeline)
library(ggpubr)
library(ggthemes)
library(dplyr)
plotFolder <- "/export/valenfs/projects/Hakon/AdamVienna/plots/new_plots/"
setwd("/export/valenfs/projects/Hakon/AdamVienna/")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load kozak grouping
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
ready <- readRDS("expression_both_withGorilla_filteredRFPSSU.rds")
ready_ER_filtered <- ready[fraction == "GO_ER",]
dim(ready_ER_filtered)
# filter is this: ready_ER_filtered <- ready_ER[complete_CDS_totalRNA_FPKM >=10& complete_leader_SSU_FPKM >=0& complete_CDS_RFP_FPKM >=0& overlapping_gene==FALSE& leader_potentially_overlaps_upstream_gene==FALSE& gene_potentially_overlaps_downstream_leader==FALSE& gene_overlaps_non_coding_transcript==FALSE& histone==FALSE& Initiation_rate_RPF > 0& !is.na(Initiation_rate_RPF)& is.finite(Initiation_rate_RPF)& leader_length >= 100, ]

seqs <- ready_ER_filtered$initiation_sequence
rate <- ready_ER_filtered$Initiation_rate_RPF
start <- 6; stop <- 14;  center <- ceiling((stop - start + 1)/2)
min.observations <- ">q1"; skip.startCodon = T; type = "IR"
plot_matrix2_log <- kozakHeatmap(seqs, rate,
                                 start, stop, center,
                                 min.observations, skip.startCodon, type = type)
plot_matrix2_log
ggsave(paste0(plotFolder, "IR_by_nucleotide_10FPKM_100nt_leaders_count1000_median.png"), plot_matrix2_log, height=100, width=250, units = 'mm', limitsize = F)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# New RNA-seq
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Rerun from /export/valenfs/projects/adam/TCP_seq/RCP_files/RCP_plots/initiation_plots/initiation_sequences_4_upstream_context.R
ready_adam <- setDT(readRDS("kozakTxWithGo.rds"))
dim(ready_adam) # <- this should be final size
adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
adams <- adams[, .(complete_CDS_totalRNA_FPKM_new, stage, transcript_id)]
ready_adam <- data.table::merge.data.table(adams, ready_adam, all.y = TRUE, all.x = FALSE, by = c('stage','transcript_id'))

dim(ready_adam)
summary(ready_adam$complete_CDS_totalRNA_FPKM_new)
summary(ready_adam$complete_CDS_totalRNA_FPKM)
rnaFilt <- 10
rnaFiltNew <- 1
ready_filtered <- ready_adam[complete_CDS_totalRNA_FPKM_new >=rnaFiltNew& complete_leader_SSU_FPKM >=0& complete_CDS_RFP_FPKM >=0& overlapping_gene==FALSE& leader_potentially_overlaps_upstream_gene==FALSE& gene_potentially_overlaps_downstream_leader==FALSE& gene_overlaps_non_coding_transcript==FALSE& histone==FALSE& Initiation_rate_RPF > 0& !is.na(Initiation_rate_RPF)& is.finite(Initiation_rate_RPF)& leader_length >= 100, ]
dim(ready_filtered)

plot_matrix2_log <- kozakHeatmap(ready_filtered$initiation_sequence,
                                 ready_filtered$Initiation_rate_RPF,
                                 start, stop, center,
                                 min.observations, skip.startCodon, type = type)
plot_matrix2_log
ggsave(paste0(plotFolder, "IR_original_by_nucleotide_10FPKM_100nt_leaders_count1000_median.png"), plot_matrix2_log, height=100, width=250, units = 'mm', limitsize = F)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Boxplot difference in groups
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# General
ready_filtered$fraction <- "Other"
ready_filtered[ready_filtered$transcript_id %in% ready_ER_filtered$transcript_id,]$fraction <- "ER"
ready_filtered <- ready_filtered[complete_leader_SSU_FPKM >0 & complete_CDS_RFP_FPKM > 0, ]
dt <- ready_filtered
dt$Scanning_efficiency_new <- dt$complete_leader_SSU_FPKM / dt$complete_CDS_totalRNA_FPKM_new
dt$TE_new <- dt$complete_CDS_RFP_FPKM    / dt$complete_CDS_totalRNA_FPKM_new
dt$TE_LSU_new <- dt$complete_CDS_LSU_FPKM    / dt$complete_CDS_totalRNA_FPKM_new
dtt <- dt %>% group_by(initiation_sequence_sub, perfect_kozak, upstream_kozak_strength) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), mean_IR_RFP = mean(Initiation_rate_RPF), median_IR_LSU = median(Initiation_rate), sequence_count = n(), median_TE=median(TE), median_TE_leader = median(Leader_RFP_translating_efficiency), median_TE_new = median(TE_new), median_TE_LSU_new = median(TE_LSU_new), median_SE =median(Scanning_efficiency), median_SE_new = median(Scanning_efficiency_new), median_RFP = median(CDS_RFP_FPKM), median_RNA_NEW = median(complete_CDS_totalRNA_FPKM_new), median_RFP_leader = median(leader_RFP_FPKM), median_SSU = median(CDS_SSU_FPKM), median_SSU_leader = median(leader_SSU_FPKM), median_LSU = median(CDS_LSU_FPKM), median_LSU_leader = median(leader_LSU_FPKM), IR_SD=sd(Initiation_rate_RPF))
dttt <- merge(dt, dtt, by = "initiation_sequence_sub")
dttt$IR_dif_ER_Other <- dttt$Initiation_rate_RPF - dttt$median_IR_RFP

ggplot(data = ready_filtered, aes(x = fraction, y = Initiation_rate_RPF)) +
  geom_boxplot() +
  ylim(0, 20)

# Relative to kozak median
tile_IE_median_20 <- ggplot(data=dtt, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR_RFP))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  #scale_fill_viridis(discrete=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
tile_IE_median_20

# box
ggplot(data = dttt, aes(x = fraction, y = IR_dif_ER_Other)) +
  geom_boxplot() +
  ylim(-5, 5) +
  ylab("IR: singles - median ")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Other plots (violin, ecdf..)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# violin
violin <- ggplot(data = dttt, aes(x = fraction, y = IR_dif_ER_Other, color = fraction)) +
  geom_boxplot(alpha = 0.05) +
  geom_violin(alpha = 0.9) +
  ylim(-5, 5) +
  ylab("IR: singles - median ")
violin
ggsave(paste0(plotFolder, "IR_relative_difference.png"), violin, height=100, width=250, units = 'mm', limitsize = F)
# dotplot
dot <- ggplot(data = dttt, aes(x = median_IR_RFP, y = Initiation_rate_RPF, color = fraction)) +
  geom_point(alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ fraction) +
  stat_regline_equation()
dot
ggsave(paste0(plotFolder, "IR_relative_difference.png"), violin, height=100, width=250, units = 'mm', limitsize = F)
# ecdf
dot <- ggplot(data = dttt, aes(x = log2(Initiation_rate_RPF/median_IR_RFP), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed IR / Expected IR)") +
  ylab("Cumulative frequency")

dot
ggsave(paste0(plotFolder, "IR_ECDF_LOG2.pdf"), dot, height=100, width=250, units = 'mm', limitsize = F, dpi = 300)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Statistical tests
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
other <- dttt[fraction == "Other",]
ER <- dttt[fraction == "ER",]
ER <- ER$Initiation_rate_RPF/ER$median_IR_RFP
other <- other$Initiation_rate_RPF/other$median_IR_RFP
t.test(ER, other)
mean(ER) / mean(other)
median(ER) / median(other)

library(Bolstad2)
sintegral(seq(440), ecdf(ER)(1:440), n.pts = 440)

library(pROC)
# Compute roc
res.roc <- roc(observed.classes, )
plot.roc(res.roc, print.auc = TRUE)

# with cyto

candidates <- unique(other$go)[grep(pattern = "cyto", x = unique(other$go))]
candidates <- candidates[7]
cytosol <- other[go %in% candidates,]
cytosol$fraction <- "Cytosol"
other2 <- other[!(go %in% candidates),]
final <- rbindlist(list(other2, ER, cytosol))

dot <- ggplot(data = final, aes(x = log2(Initiation_rate_RPF/median_IR_RFP), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed IR / Expected IR)") +
  ylab("Cumulative frequency")
dot

# IR test ER vs median by kozak
# ER

irTest <- ggscatter(ER[sequence_count > 20,], x = "median_IR_RFP", y = "Initiation_rate_RPF",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Median ER per Kozak", ylab = "Observer IR")
irTest <- irTest +
  scale_x_log10() +
  scale_y_log10()

irTest
ggsave(paste0(plotFolder, "IR_dotplot_ER.pdf"), irTest, height=100, width=250, units = 'mm', limitsize = F, dpi = 300)

cor.test(ER[sequence_count > 20,]$Initiation_rate_RPF, ER[sequence_count > 20,]$median_IR_RFP, method = "spearman")

# Other
irTest <- ggscatter(other[sequence_count > 20,], x = "median_IR_RFP", y = "Initiation_rate_RPF",
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = TRUE, cor.method = "spearman",
                    xlab = "Median ER per Kozak", ylab = "Observer IR", )
irTest <- irTest +
  scale_x_log10() +
  scale_y_log10()
irTest

ggscatter(ER, aes(x = log10(Initiation_rate_RPF), y =log10(median_IR_RFP))) +
  geom_point()

ggsave(paste0(plotFolder, "IR_dotplot_other.pdf"), irTest, height=100, width=250, units = 'mm', limitsize = F, dpi = 300)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# New tests ecdf plots (sanity tests only)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# WITH TE, for figure 4e check
# IR with RFP:
dot <- ggplot(data = dttt, aes(x = log2(Initiation_rate_RPF/median_IR_RFP), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed IR / Expected IR)") +
  ylab("Cumulative frequency")
dot
# IR with LSU:
dot <- ggplot(data = dttt, aes(x = log2(Initiation_rate/median_IR_LSU), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed IR LSU / Expected IR LSU)") +
  ylab("Cumulative frequency")
dot

# TE defined as: complete_CDS_RFP_FPKM    / complete_CDS_totalRNA_FPKM)
dot <- ggplot(data = dttt, aes(x = log2(TE/median_TE), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed TE / Expected TE)") +
  ylab("Cumulative frequency")
dot
# TE new RNA:
dot <- ggplot(data = dttt, aes(x = log2(TE_new/median_TE_new), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed new TE / Expected TE)") +
  ylab("Cumulative frequency")
dot
# TE new RNA and LSU:
dot <- ggplot(data = dttt, aes(x = log2(TE_LSU_new/median_TE_LSU_new), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed new TE by LSU / Expected TE by LSU)") +
  ylab("Cumulative frequency")
dot

# leader TE defined as:
dot <- ggplot(data = dttt, aes(x = log2(Leader_RFP_translating_efficiency/median_TE_leader), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed TE leader / Expected TE leader)") +
  ylab("Cumulative frequency")
dot
#ggsave(paste0(plotFolder, "IR_ECDF_LOG2.pdf"), dot, height=100, width=250, units = 'mm', limitsize = F, dpi = 300)

# SE:
dot <- ggplot(data = dttt, aes(x = log2(Scanning_efficiency/median_SE), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed SE / Expected SE)") +
  ylab("Cumulative frequency")
dot
# SE new:
dot <- ggplot(data = dttt, aes(x = log2(Scanning_efficiency_new/median_SE_new), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed SE / Expected SE)") +
  ylab("Cumulative frequency")
dot

# cds RFP:
dot <- ggplot(data = dttt, aes(x = log2(CDS_RFP_FPKM/median_RFP), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed RFP / Expected RFP)") +
  ylab("Cumulative frequency")

dot

# leader RFP:
dot <- ggplot(data = dttt, aes(x = log2(leader_RFP_FPKM/median_RFP_leader), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed leader RFP / Expected leader RFP)") +
  ylab("Cumulative frequency")
dot

# RNA:
dot <- ggplot(data = dttt, aes(x = log2(complete_CDS_totalRNA_FPKM_new/median_RNA_NEW), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed RNA / Expected RNA)") +
  ylab("Cumulative frequency")
dot
# SSU:
dot <- ggplot(data = dttt, aes(x = log2(CDS_SSU_FPKM/median_SSU), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed RNA / Expected RNA)") +
  ylab("Cumulative frequency")
dot
# SSU leader:
dot <- ggplot(data = dttt, aes(x = log2(leader_SSU_FPKM/median_SSU_leader), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed SSU leader / Expected SSU leader)") +
  ylab("Cumulative frequency")
dot
# LSU:
dot <- ggplot(data = dttt, aes(x = log2(CDS_LSU_FPKM/median_LSU), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed LSU / Expected LSU)") +
  ylab("Cumulative frequency")
dot
# LSU leader:
dot <- ggplot(data = dttt, aes(x = log2(leader_LSU_FPKM/median_LSU_leader), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (Observed LSU leader / Expected LSU leader)") +
  ylab("Cumulative frequency")
dot
# RNA:
dot <- ggplot(data = dttt, aes(x = log2(complete_CDS_totalRNA_FPKM_new/complete_CDS_totalRNA_FPKM), color = fraction)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_base() +
  xlab("Log2 (new RNA / old RNA)") +
  ylab("Cumulative frequency")
dot
#ggsave(paste0(plotFolder, "IR_ECDF_LOG2.pdf"), dot, height=100, width=250, units = 'mm', limitsize = F, dpi = 300)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Correlation plots for ER
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#
library(GGally)
df <- dttt[, .(complete_CDS_totalRNA_FPKM, complete_CDS_totalRNA_FPKM_new, fraction)]
d <- df[fraction == "Other",]
d$fraction <- NULL
ggpairs(d)
d <- df[fraction == "ER",]
d$fraction <- NULL
ggpairs(d)
