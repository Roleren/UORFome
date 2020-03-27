#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Old ways of doing it
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# LOAD DATA
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load new RNA-seq (This is old RNA_seq fpkm, not used anymore)
df <- read.experimentl("Val19", expInVarName = FALSE)
txdb <- loadTxdb(df)
loadRegions(txdb, parts = c("mrna", "leaders", "cds", "trailers", "tx"));
outputLibs(df, leaders)

data_for_pairs <- makeSummarizedExperimentFromBam(df, region = mrna, geneOrTxNames = "tx", longestPerGene = FALSE,
                                                  type = "fpkm", saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_mrna.rds")
leader_fpkms <- makeSummarizedExperimentFromBam(df, region = leaders, geneOrTxNames = "tx", longestPerGene = FALSE,
                                                type = "fpkm", saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_leaders.rds")
cds_fpkms <- makeSummarizedExperimentFromBam(df, region = cds, geneOrTxNames = "tx", longestPerGene = FALSE,
                                             type = "fpkm", saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_mrna.rds")
trailer_fpkms <- makeSummarizedExperimentFromBam(df, region = trailers, geneOrTxNames = "tx", longestPerGene = FALSE,
                                                 type = "fpkm", saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_trailers.rds")
leader_fpkms[,transcript_id := names(leaders)]
cds_fpkms[,transcript_id := names(cds)]
trailer_fpkms[,transcript_id := names(trailers)]

cds_fpkmsNew <- rbindlist(list(cds_fpkms[,.(transcript_id, complete_CDS_totalRNA_FPKM_new = RNA__64cell_f1)],
                               cds_fpkms[,.(transcript_id, complete_CDS_totalRNA_FPKM_new = RNA__Sphere_f1)],
                               cds_fpkms[,.(transcript_id, complete_CDS_totalRNA_FPKM_new = RNA__Shield_f1)]))
cds_fpkmsNew$stage <- as.factor(c(rep("64cell", nrow(cds_fpkms)), rep("sphere", nrow(cds_fpkms)), rep("shield", nrow(cds_fpkms))))


# Adams list, it filters on all 3 stages at the same time.
adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered.rds")
setDT(adams)
a <- copy(adams)
a$transcript_id <- as.character(a$transcript_id)
adams <- data.table::merge.data.table(adams, cds_fpkmsNew, all = TRUE, by = c('stage','transcript_id'))
adams <- adams[!is.na(complete_leader_RNA_FPKM),]
#saveRDS(object = adams,"/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")


################################ WITH LSU PER STAGE ##########################################
# With RFP you get dim: 103 groups using 20, with LSU  85, we use 13 to get 102 groups (USED IN ARTICLE!)
combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM_new >=rnaFiltNew, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate > 0, !is.na(Initiation_rate), is.finite(Initiation_rate), leader_length >= 100)
combined.dat.ranked.filter$initiation_sequence_sub <- substring(combined.dat.ranked.filter$initiation_sequence, 6, 9)

combined.dat.ranked.filter$perfect_kozak <- NA
combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$perfect_kozak <- 1
kozak_x_lables <-combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$initiation_sequence_sub

combined.dat.ranked.filter.sequence.group <- combined.dat.ranked.filter %>% group_by(initiation_sequence_sub, perfect_kozak, upstream_kozak_strength, stage) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), median_IR = median(Initiation_rate), mean_IR_RFP = mean(Initiation_rate_RPF), mean_IR = mean(Initiation_rate), sequence_count = n(), median_TE=median(TE), IR_SD=sd(Initiation_rate_RPF))

combined.dat.ranked.filter.sequence.group <- filter (combined.dat.ranked.filter.sequence.group, sequence_count >= 5)
dim(combined.dat.ranked.filter.sequence.group)
combined.dat.ranked.filter.sequence.group$initiation_sequence_sub <- factor(combined.dat.ranked.filter.sequence.group$initiation_sequence_sub, levels = combined.dat.ranked.filter.sequence.group$initiation_sequence_sub[order(combined.dat.ranked.filter.sequence.group$median_IR)])


median_IR_RFP_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR_RFP))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid( ~ stage) +
  NULL
median_IR_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid( ~ stage) +
  NULL
median_IR_RFP_20
median_IR_20
ggsave(p(plotFolder, "RFP_PERSTAGE_bar_median_IR_20.pdf"), median_IR_RFP_20, dpi = 300, width = 30, heigh = 5)
ggsave(p(plotFolder, "LSU_PERSTAGE_bar_median_IR_20.pdf"), median_IR_20, dpi = 300, width = 30, heigh = 5)
combined.dat.ranked.filter.sequence.group[combined.dat.ranked.filter.sequence.group$initiation_sequence_sub %in% c("AAGC", "AAAC", "TGGA"),]

combined.dat.ranked.filter.sequence.group[combined.dat.ranked.filter.sequence.group$initiation_sequence_sub == c("AAGC", "AAAC", "TGGA"),]
tile_KS_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=upstream_kozak_strength)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("#146DA8", "white", "#A87225")) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
ggsave(p(plotFolder, "bar_tile_KS_20_new.pdf"), tile_KS_20, dpi = 300, width = 30, heigh = 5)

################################ WITH LSU PER STAGE grouped ##########################################
combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM_new >=rnaFiltNew, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate > 0, !is.na(Initiation_rate), is.finite(Initiation_rate), leader_length >= 100)
combined.dat.ranked.filter$initiation_sequence_sub <- substring(combined.dat.ranked.filter$initiation_sequence, 6, 9)

combined.dat.ranked.filter$perfect_kozak <- NA
combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$perfect_kozak <- 1
kozak_x_lables <-combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$initiation_sequence_sub

combined.dat.ranked.filter.sequence.group <- combined.dat.ranked.filter %>% group_by(initiation_sequence_sub, perfect_kozak, upstream_kozak_strength, stage) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), median_IR = median(Initiation_rate), mean_IR_RFP = mean(Initiation_rate_RPF), mean_IR = mean(Initiation_rate), sequence_count = n(), median_TE=median(TE), IR_SD=sd(Initiation_rate_RPF))

combined.dat.ranked.filter.sequence.group <- filter (combined.dat.ranked.filter.sequence.group, sequence_count >= 5)
dim(combined.dat.ranked.filter.sequence.group)
# Mean grouping
groups <- as.data.table(combined.dat.ranked.filter.sequence.group)
groups <- groups[, .(mean_IR = mean(mean_IR), median_IR = median(median_IR), sequence_count = sum(sequence_count)), by = .(initiation_sequence_sub)]
LSU_ranking_new <- groups[groups$initiation_sequence_sub %in% c("AAGC", "AAAC", "TGGA"),][c(2,1,3)];LSU_ranking_new
groups$initiation_sequence_sub <- factor(groups$initiation_sequence_sub, levels = groups$initiation_sequence_sub[order(groups$mean_IR)])
mean_IR_20 <- ggplot(data=groups, aes(x=initiation_sequence_sub, y=1, fill=log2(mean_IR))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mean_IR_20
# Median grouping
groups$initiation_sequence_sub <- factor(groups$initiation_sequence_sub, levels = groups$initiation_sequence_sub[order(groups$median_IR)])
median_IR_20 <- ggplot(data=groups, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
median_IR_20

groups <- as.data.table(combined.dat.ranked.filter.sequence.group)
groups <- groups[, .(mean_IR = mean(mean_IR), median_IR = median(median_IR), upstream_kozak_strength = mean(upstream_kozak_strength)), by = .(initiation_sequence_sub)]
groups$initiation_sequence_sub <- factor(groups$initiation_sequence_sub, levels = groups$initiation_sequence_sub[order(groups$mean_IR)])
tile_KS_20 <- ggplot(data=groups, aes(x=initiation_sequence_sub, y=1, fill=upstream_kozak_strength)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
tile_KS_20

ggsave(p(plotFolder, "LSU_PERSTAGE_grouped_bar_mean_IR_20.pdf"), mean_IR_20, dpi = 300, width = 30, heigh = 5)
ggsave(p(plotFolder, "LSU_PERSTAGE_grouped_bar_median_IR_20.pdf"), median_IR_20, dpi = 300, width = 30, heigh = 5)
ggsave(p(plotFolder, "LSU_PERSTAGE_grouped_KS_20.pdf"), tile_KS_20, dpi = 300, width = 30, heigh = 5)

# Kozak / IR correlation
mat <- groups[, c("upstream_kozak_strength", "median_IR")]
corPlot <- ggplot(data = mat, aes(x = upstream_kozak_strength, y = log2(median_IR))) +
  geom_point() +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_classic()
corPlot
ggsave(p(plotFolder, "corPlot_kozak_IR_Ranking_meanofmeans.pdf"), corPlot)


################################### Kozak heatmap median
plotFolder <- "/export/valenfs/projects/Hakon/RCP_SEQ/plots/new_IR_plots/"
rnaFilt <- 10
rnaFiltNew <- 5 # Used filters: 6.14, 1
combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM >= rnaFilt, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate_RPF > 0, !is.na(Initiation_rate_RPF), is.finite(Initiation_rate_RPF), leader_length >= 100)

seqs <- combined.dat.ranked.filter$initiation_sequence
rate <- combined.dat.ranked.filter$Initiation_rate_RPF
start <- 6; stop <- 14;  center <- ceiling((stop - start + 1)/2)
min.observations <- ">q1"; skip.startCodon = T; type = "IR"
oldHM <- ORFik:::kozakHeatmap(seqs, rate,
                              start, stop, center,
                              min.observations, skip.startCodon, type = type)

# new version
combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM_new >=rnaFiltNew, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate_RPF > 0, !is.na(Initiation_rate_RPF), is.finite(Initiation_rate_RPF), leader_length >= 100)
seqs <- combined.dat.ranked.filter$initiation_sequence
rate <- combined.dat.ranked.filter$Initiation_rate_RPF
newHM <- ORFik:::kozakHeatmap(seqs, rate,
                              start, stop, center,
                              min.observations, skip.startCodon, type = type)
oldHM
newHM
ggsave(p(plotFolder, "old_hm_byIR.pdf"), oldHM, dpi = 300)
ggsave(p(plotFolder, "new_hm_byIR.pdf"), newHM, dpi = 300)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# IR ranking median split by kozak
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM >=rnaFilt, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate_RPF > 0, !is.na(Initiation_rate_RPF), is.finite(Initiation_rate_RPF), leader_length >= 100)
combined.dat.ranked.filter$initiation_sequence_sub <- substring(combined.dat.ranked.filter$initiation_sequence, 6, 9) # was 12

combined.dat.ranked.filter$perfect_kozak <- NA
combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$perfect_kozak <- 1
kozak_x_lables <-combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$initiation_sequence_sub

combined.dat.ranked.filter.sequence.group <- combined.dat.ranked.filter %>% group_by(initiation_sequence_sub, perfect_kozak, upstream_kozak_strength) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), mean_IR_RFP = mean(Initiation_rate_RPF), sequence_count = n(), median_TE=median(TE), IR_SD=sd(Initiation_rate_RPF))

combined.dat.ranked.filter.sequence.group <- filter (combined.dat.ranked.filter.sequence.group, sequence_count >= 20)
dim(combined.dat.ranked.filter.sequence.group)

combined.dat.ranked.filter.sequence.group$initiation_sequence_sub <- factor(combined.dat.ranked.filter.sequence.group$initiation_sequence_sub, levels = combined.dat.ranked.filter.sequence.group$initiation_sequence_sub[order(combined.dat.ranked.filter.sequence.group$median_IR_RFP)])
old <- combined.dat.ranked.filter.sequence.group

median_IR_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR_RFP))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
tile_KS_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=upstream_kozak_strength)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
median_IR_20
tile_KS_20
# New
combined.dat.ranked.filter <- adams %>% filter(complete_CDS_totalRNA_FPKM_new >=rnaFiltNew, complete_leader_SSU_FPKM >=0, complete_CDS_RFP_FPKM >=0, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, histone==FALSE, Initiation_rate_RPF > 0, !is.na(Initiation_rate_RPF), is.finite(Initiation_rate_RPF), leader_length >= 100)
combined.dat.ranked.filter$initiation_sequence_sub <- substring(combined.dat.ranked.filter$initiation_sequence, 6, 9)

combined.dat.ranked.filter$perfect_kozak <- NA
combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$perfect_kozak <- 1
kozak_x_lables <-combined.dat.ranked.filter[combined.dat.ranked.filter$upstream_kozak_strength > -7.5,]$initiation_sequence_sub

combined.dat.ranked.filter.sequence.group <- combined.dat.ranked.filter %>% group_by(initiation_sequence_sub, perfect_kozak, upstream_kozak_strength) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), mean_IR_RFP = mean(Initiation_rate_RPF), sequence_count = n(), median_TE=median(TE), IR_SD=sd(Initiation_rate_RPF))
a <- (combined.dat.ranked.filter %>% group_by(initiation_sequence_sub, stage) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), median_IR = median(Initiation_rate), mean_IR = mean(Initiation_rate), sequence_count = n())); a[a$initiation_sequence_sub %in% c("AAGC", "AAAC", "TGGA"),]
combined.dat.ranked.filter.sequence.group <- filter (combined.dat.ranked.filter.sequence.group, sequence_count >= 20)
dim(combined.dat.ranked.filter.sequence.group)
combined.dat.ranked.filter.sequence.group$initiation_sequence_sub <- factor(combined.dat.ranked.filter.sequence.group$initiation_sequence_sub, levels = combined.dat.ranked.filter.sequence.group$initiation_sequence_sub[order(combined.dat.ranked.filter.sequence.group$median_IR_RFP)])


median_IR_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR_RFP))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
tile_KS_20 <- ggplot(data=combined.dat.ranked.filter.sequence.group, aes(x=initiation_sequence_sub, y=1, fill=upstream_kozak_strength)) +
  geom_tile() +
  # scale_fill_gradientn(colors = c("#146DA8", "white", "#A87225")) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
ggsave(p(plotFolder, "bar_median_IR_20.pdf"), median_IR_20, dpi = 300, width = 20, heigh = 5)
ggsave(p(plotFolder, "bar_tile_KS_20.pdf"), tile_KS_20, dpi = 300, width = 20, heigh = 5)

# Check new ordering
b <- as.character(old[order(old$median_IR_RFP),]$initiation_sequence_sub)
d <- as.character(combined.dat.ranked.filter.sequence.group[order(combined.dat.ranked.filter.sequence.group$median_IR_RFP),]$initiation_sequence_sub)
length(d);length(b)
b <- b[b %in% d]
all(d %in% b)
sum(b == d)
all(b == d)
# The 3 experimentally verified ones
"AAGC" %in% b
"AAAC" %in% b
"TGGA" %in% b
print("low: AAGC");paste("RFP old:",which(b == "AAGC"));paste("RFP:",which(d == "AAGC"))
print("middle: AAAC");paste("RFP old:",which(b == "AAAC"));paste("RFP:",which(d == "AAAC"))
print("high TGGA");paste("RFP old:",which(b == "TGGA"));paste("RFP:",which(d == "TGGA"))
