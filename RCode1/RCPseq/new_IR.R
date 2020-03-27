#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO (used in article)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Used to create new IR figures

#/export/valenfs/projects/adam/TCP_seq/RCP_files/RCP_plots/initiation_plots/initiation_sequences_4_upstream_context.R
# We are now using LSU for article
library(ORFikPipeline)
library(dplyr)
library(ggpubr)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Run analysis: Start here, if making rna-seq is done!
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#
adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
# statistics
cor.test(adams$complete_CDS_RNA_FPKM, adams$complete_CDS_totalRNA_FPKM_new)
summary(adams$complete_CDS_totalRNA_FPKM)
summary(adams$complete_CDS_totalRNA_FPKM_new)
any(is.na(adams$complete_CDS_totalRNA_FPKM_new))
# Justification for 10 as filter
quantile(adams$complete_CDS_totalRNA_FPKM, 0.645)
quantile(adams$complete_CDS_totalRNA_FPKM_new, 0.645)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# With LSU (THIS way is used now in article)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# With RFP you get dim: 103 groups using 20, with LSU  85, we use 20 to get 102 groups
combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM_new >=rnaFiltNew, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate > 0, !is.na(Initiation_rate), is.finite(Initiation_rate), leader_length >= 100)
combined.dat.ranked.filter$initiation_sequence_sub <- substring(combined.dat.ranked.filter$initiation_sequence, 6, 9)

combined.dat.ranked.filter$perfect_kozak <- NA
combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$perfect_kozak <- 1
kozak_x_lables <-combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$initiation_sequence_sub

combined.dat.ranked.filter.sequence.group <- combined.dat.ranked.filter %>% group_by(initiation_sequence_sub, perfect_kozak, upstream_kozak_strength) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), median_IR = median(Initiation_rate), mean_IR_RFP = mean(Initiation_rate_RPF), mean_IR = mean(Initiation_rate), sequence_count = n(), median_TE=median(TE), IR_SD=sd(Initiation_rate_RPF))
a <- (combined.dat.ranked.filter %>% group_by(initiation_sequence_sub, stage) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), median_IR = median(Initiation_rate), mean_IR = mean(Initiation_rate), sequence_count = n())); a[a$initiation_sequence_sub %in% c("AAGC", "AAAC", "TGGA"),]

combined.dat.ranked.filter.sequence.group <- filter (combined.dat.ranked.filter.sequence.group, sequence_count >= 20)
dim(combined.dat.ranked.filter.sequence.group)
combined.dat.ranked.filter.sequence.group$initiation_sequence_sub <- factor(combined.dat.ranked.filter.sequence.group$initiation_sequence_sub, levels = combined.dat.ranked.filter.sequence.group$initiation_sequence_sub[order(combined.dat.ranked.filter.sequence.group$median_IR)])
LSU_ranking_old <- combined.dat.ranked.filter.sequence.group[combined.dat.ranked.filter.sequence.group$initiation_sequence_sub %in% c("AAGC", "AAAC", "TGGA"),][c(2,1,3),];LSU_ranking_old[,c("initiation_sequence_sub", "mean_IR", "median_IR", "sequence_count")]

median_IR_RFP_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR_RFP))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
median_IR_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
median_IR_RFP_20
median_IR_20


tile_KS_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=upstream_kozak_strength)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("#146DA8", "white", "#A87225")) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
ggsave(p(plotFolder, "bar_median_IR_20.pdf"), median_IR_20, dpi = 300, width = 20, heigh = 5)
ggsave(p(plotFolder, "bar_tile_KS_20.pdf"), tile_KS_20, dpi = 300, width = 20, heigh = 5)
# mean
combined.dat.ranked.filter.sequence.group$initiation_sequence_sub <- factor(combined.dat.ranked.filter.sequence.group$initiation_sequence_sub, levels = combined.dat.ranked.filter.sequence.group$initiation_sequence_sub[order(combined.dat.ranked.filter.sequence.group$mean_IR)])
mean_IR_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=log2(mean_IR))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "blue1", "white", "red1", "red2", "red3")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "IR")
tile_KS_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=upstream_kozak_strength)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("#146DA8", "white", "#A87225")) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Kozak")
mean_IR_20
ggslackR(width = 300, height = 75)
tile_KS_20
ggslackR(width = 300, height = 75)
ggsave(p(plotFolder, "LSU_bar_tile_by_mean_IR_20.pdf"), mean_IR_20, dpi = 300, width = 20, heigh = 5)
ggsave(p(plotFolder, "LSU_bar_tile_by_mean_KS_20.pdf"), tile_KS_20, dpi = 300, width = 20, heigh = 5)

# New ordering
b <- as.character(old[order(old$median_IR_RFP),]$initiation_sequence_sub)
d <- as.character(combined.dat.ranked.filter.sequence.group[order(combined.dat.ranked.filter.sequence.group$median_IR),]$initiation_sequence_sub)
length(d);length(b)
b <- b[b %in% d]
all(d %in% b)
sum(b == d)
all(b == d)
# The 3 experimentally verified ones
"AAGC" %in% b
"AAAC" %in% b
"TGGA" %in% b
print("low: AAGC");paste("RFP:",which(b == "AAGC"));paste("LSU:",which(d == "AAGC"))
print("middle: AAAC");paste("RFP:",which(b == "AAAC"));paste("LSU:",which(d == "AAAC"))
print("high TGGA");paste("RFP:",which(b == "TGGA"));paste("LSU:",which(d == "TGGA"))

# Correlation LSU vs RFP
mat <- combined.dat.ranked.filter.sequence.group[, c("median_IR_RFP", "median_IR")]
corPlot <- ggpairs(mat); corPlot
ggsave(p(plotFolder, "corPlot_IR_Ranking.pdf"), corPlot)
# Kozak strength cor rfp
mat <- combined.dat.ranked.filter.sequence.group[, c("upstream_kozak_strength", "median_IR_RFP")]
mat$median_IR_RFP <- log2(mat$median_IR_RFP)
corPlot <- ggpairs(mat); corPlot
ggsave(p(plotFolder, "corPlot_kozakvsIR_Old_Ranking.pdf"), corPlot)
# Kozak strength cor LSU median
mat <- combined.dat.ranked.filter.sequence.group[, c("upstream_kozak_strength", "median_IR")]
corPlot <- ggplot(data = mat, aes(x = upstream_kozak_strength, y = log2(median_IR))) +
  geom_point() +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_classic()
corPlot
ggsave(p(plotFolder, "corPlot_kozak_IR_Ranking.pdf"), corPlot)
# Kozak strength cor LSU mean (used in paper!)
mat <- combined.dat.ranked.filter.sequence.group[, c("initiation_sequence_sub", "upstream_kozak_strength", "mean_IR")]
colors <- rep("gray", nrow(mat))
colors[which(mat$initiation_sequence_sub == "TGGA")] <- "red"; colors[which(mat$initiation_sequence_sub == "AAAC")] <- "black";colors[which(mat$initiation_sequence_sub == "AGCC")] <- "blue"

corPlot <- ggplot(data = mat, aes(x = upstream_kozak_strength, y = log2(mean_IR))) +
  geom_point(size = 3, color = colors) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -5.5) +
  theme_classic() +
  xlab("Kozak Strength") +
  ylab("log2(mean IR)") +
  theme(axis.text=element_text(size=12))
corPlot
ggsave(p(plotFolder, "corPlot_kozak_IR_Ranking_mean.pdf"), corPlot)
ggslackR(width = 200, height = 100)

