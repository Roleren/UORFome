############################ UORF BARS #######################################################################
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(cowplot)
library(ggpubr)
library(viridis)

raw.dat.64cell<- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/matrix_files/64cell_WT_reps_1_2_peaks_removed_translating_filter.csv", header = TRUE, sep = ",", data.table = F)
raw.dat.sphere<-fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/matrix_files/sphere_WT_reps_1_2_3_peaks_removed_translating_filter.csv", header = TRUE, sep = ",", data.table = F)
raw.dat.shield<-fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/matrix_files/shield_WT_reps_1_2_3_peaks_removed_translating_filter.csv", header = TRUE, sep = ",", data.table = F)

#uORF counts
uorf.ATG.64cell<-fread("/export/valenfs/projects/adam/TCP_seq/Bams_for_eivind_16_09_18/uORF_count/64cell_most_highly_expressed_ATG_uorf_counts.csv", header = TRUE, sep = ",", data.table = F)
uorf.ATG.sphere<-fread("/export/valenfs/projects/adam/TCP_seq/Bams_for_eivind_16_09_18/uORF_count/Sphere_most_highly_expressed_ATG_uorf_counts.csv", header = TRUE, sep = ",", data.table = F)
uorf.ATG.shield<-fread("/export/valenfs/projects/adam/TCP_seq/Bams_for_eivind_16_09_18/uORF_count/Shield_most_highly_expressed_ATG_uorf_counts.csv", header = TRUE, sep = ",", data.table = F)

raw.dat.64cell$ATG_count <- 0
raw.dat.sphere$ATG_count <- 0
raw.dat.shield$ATG_count <- 0

raw.dat.64cell$ATG_count <- uorf.ATG.64cell$ATG_uORF_count[ match(raw.dat.64cell$`#gene`, uorf.ATG.64cell$`#gene`)]
raw.dat.sphere$ATG_count <- uorf.ATG.sphere$ATG_uORF_count[ match(raw.dat.sphere$`#gene`, uorf.ATG.sphere$`#gene`)]
raw.dat.shield$ATG_count <- uorf.ATG.shield$ATG_uORF_count[ match(raw.dat.shield$`#gene`, uorf.ATG.shield$`#gene`)]

raw.dat.64cell$ATG_count[is.na(raw.dat.64cell$ATG_count)] <- 0
raw.dat.sphere$ATG_count[is.na(raw.dat.sphere$ATG_count)] <- 0
raw.dat.shield$ATG_count[is.na(raw.dat.shield$ATG_count)] <- 0

#combine
raw.dat.64cell$stage <- "64cell"
raw.dat.sphere$stage <- "sphere"
raw.dat.shield$stage <- "shield"

combined.dat <- rbind(raw.dat.64cell, raw.dat.sphere, raw.dat.shield)
combined.dat$stage <- factor(combined.dat$stage, levels = c("64cell", "sphere", "shield"))

#add dervided columns
combined.dat <- combined.dat %>% mutate(TE                                  = complete_CDS_RFP_FPKM    / complete_CDS_totalRNA_FPKM)
combined.dat <- combined.dat %>% mutate(Scanning_efficiency                 = complete_leader_SSU_FPKM / complete_CDS_totalRNA_FPKM)
combined.dat <- combined.dat %>% mutate(Scanning_efficiency_exluding_5      = leader_SSU_FPKM          / complete_CDS_totalRNA_FPKM)
combined.dat <- combined.dat %>% mutate(TIS_recognition_rate_exluding_5     = Start_codon_SSU_FPKM     / leader_SSU_FPKM)
combined.dat <- combined.dat %>% mutate(TIS_recognition_rate                = Start_codon_SSU_FPKM     / (leader_SSU_FPKM+TSS_beginning_SSU_FPKM))
combined.dat <- combined.dat %>% mutate(Leader_LSU_translating_efficiency   = complete_leader_LSU_FPKM / complete_CDS_totalRNA_FPKM)
combined.dat <- combined.dat %>% mutate(Leader_RFP_translating_efficiency   = complete_leader_RFP_FPKM / complete_CDS_totalRNA_FPKM)
combined.dat <- combined.dat %>% mutate(Initiation_rate                     = complete_CDS_LSU_FPKM    / complete_leader_SSU_FPKM)
combined.dat <- combined.dat %>% mutate(Initiation_rate_RPF                 = complete_CDS_RFP_FPKM    / complete_leader_SSU_FPKM)
#combined.dat <- combined.dat %>% mutate(Initiation_rate                     = complete_CDS_LSU_FPKM    / initiation_coverage_SSU) #50nt upstream of TIS
#combined.dat <- combined.dat %>% mutate(Initiation_rate_RPF                 = complete_CDS_RFP_FPKM    / initiation_coverage_SSU) #50nt upstream of TIS
combined.dat <- combined.dat %>% mutate(SSU_TSS_TIS                         = Start_codon_SSU_FPKM     / TSS_beginning_SSU_FPKM)
combined.dat <- combined.dat %>% mutate(CDS_SSU_RFP                         = complete_CDS_RFP_FPKM    / complete_CDS_SSU_FPKM)

#Kozak
combined.dat.ranked <- combined.dat %>% arrange(desc(kozak_strength))
combined.dat.ranked$type <- "Other" #set catagories before filtering
combined.dat.ranked[1:500,]$type <- "Strong kozak"
combined.dat.ranked[(nrow(combined.dat.ranked) - 500):nrow(combined.dat.ranked),]$type <- "Weak kozak"



###################################################################################################
# Barplots of uORF counts
###################################################################################################
combined.dat.ranked$ATG_count[is.na(combined.dat.ranked$ATG_count)] <- 0
combined.dat.ranked_new <-  readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
combined.dat.ranked$complete_CDS_totalRNA_FPKM_new <- combined.dat.ranked_new$complete_CDS_totalRNA_FPKM_new
combined.dat.ranked.filtered <- filter(combined.dat.ranked,
                                       complete_CDS_totalRNA_FPKM_new >= 1,
                                       ATG_count <= 4
)

pallet_5  <- c("#146DA8", "#4EB3F5", "#FFD291", "#F5B04E", "#A87225")


uORF_IR <- ggplot(combined.dat.ranked.filtered, aes(x=ATG_count, y=log2(leader_SSU_FPKM), fill=as.factor(ATG_count))) +
  geom_boxplot(width=0.5, outlier.colour=NA) +
  scale_fill_manual(values=pallet_5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  ylab("5' UTR SSU fpkm")

uORF_leader_TE <- ggplot(combined.dat.ranked.filtered, aes(x=ATG_count, y=log2(Leader_RFP_translating_efficiency), fill=as.factor(ATG_count))) +
  geom_boxplot(width=0.5, outlier.colour=NA) +
  scale_fill_manual(values=pallet_5) +
  coord_cartesian(ylim = c(-9, 5)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  ylab("5' UTR translational efficiency")

uORF_CDS_TE <- ggplot(combined.dat.ranked.filtered, aes(x=ATG_count, y=log2(TE), fill=as.factor(ATG_count))) +
  geom_boxplot(width=0.5, outlier.colour=NA) +
  scale_fill_manual(values=pallet_5) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  ylab("Translational efficiency")

uORF_bars <- grid.arrange(uORF_IR, uORF_leader_TE, uORF_CDS_TE, nrow=1)

ggsave("/export/valenfs/projects/Hakon/RCP_SEQ/plots/new_IR_plots/uORF_bars_selected.pdf", uORF_bars, height=150, width=200, units = 'mm', limitsize = F)
ggsave("/export/valenfs/projects/Hakon/RCP_SEQ/plots/new_IR_plots/uORF_bars_selected.png", uORF_bars, height=150, width=200, units = 'mm', limitsize = F)


#################################### UORF STOP BAR PLOTS #####################################################

library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

#directionality in total-rna library updated

raw.dat.64cell <- raw.dat.64cell %>% mutate(TE = complete_CDS_RFP_FPKM / complete_CDS_totalRNA_FPKM)
raw.dat.sphere <- raw.dat.sphere %>% mutate(TE = complete_CDS_RFP_FPKM / complete_CDS_totalRNA_FPKM)
raw.dat.shield <- raw.dat.shield %>% mutate(TE = complete_CDS_RFP_FPKM / complete_CDS_totalRNA_FPKM)

raw.dat.64cell <- raw.dat.64cell %>% mutate(SE = complete_leader_SSU_FPKM / complete_CDS_totalRNA_FPKM)
raw.dat.sphere <- raw.dat.sphere %>% mutate(SE = complete_leader_SSU_FPKM / complete_CDS_totalRNA_FPKM)
raw.dat.shield <- raw.dat.shield %>% mutate(SE = complete_leader_SSU_FPKM / complete_CDS_totalRNA_FPKM)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#uORF up dpown regions, with down regions from uORF stop to TIS

uORF_stats_64cell <- read.csv("~/analysis/TCP/uORF/64cell_uorf_near_canonical_up_down_counts_from_uORF_stop_codons.csv", header = T)
uORF_stats_sphere <- read.csv("~/analysis/TCP/uORF/sphere_uorf_near_canonical_up_down_counts_from_uORF_stop_codons.csv", header = T)
uORF_stats_shield <- read.csv("~/analysis/TCP/uORF/shield_uorf_near_canonical_up_down_counts_from_uORF_stop_codons.csv", header = T)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Add in TE / RNA
uORF_stats_64cell$CDS_TE <- 0;
uORF_stats_64cell$CDS_TE <- raw.dat.64cell$TE[match(uORF_stats_64cell$gene_id, raw.dat.64cell$X.gene_id)]
uORF_stats_64cell$CDS_SE <- 0;
uORF_stats_64cell$CDS_SE <- raw.dat.64cell$SE[match(uORF_stats_64cell$gene_id, raw.dat.64cell$X.gene_id)]
uORF_stats_64cell$CDS_RFP <- 0;
uORF_stats_64cell$CDS_RFP <- raw.dat.64cell$complete_CDS_RFP_FPKM[match(uORF_stats_64cell$gene_id, raw.dat.64cell$X.gene_id)]
uORF_stats_64cell$CDS_RNA <- 0;
uORF_stats_64cell$CDS_RNA <- raw.dat.64cell$complete_CDS_totalRNA_FPKM[match(uORF_stats_64cell$gene_id, raw.dat.64cell$X.gene_id)]
uORF_stats_64cell$leader_SSU <- 0;
uORF_stats_64cell$leader_SSU <- raw.dat.64cell$complete_leader_SSU_FPKM[match(uORF_stats_64cell$gene_id, raw.dat.64cell$X.gene_id)]
uORF_stats_64cell$leader_length <- 0;
uORF_stats_64cell$leader_length <- raw.dat.64cell$leader_length[match(uORF_stats_64cell$gene_id, raw.dat.64cell$X.gene_id)]

uORF_stats_sphere$CDS_TE <- 0;
uORF_stats_sphere$CDS_TE <- raw.dat.sphere$TE[match(uORF_stats_sphere$gene_id, raw.dat.sphere$X.gene_id)]
uORF_stats_sphere$CDS_SE <- 0;
uORF_stats_sphere$CDS_SE <- raw.dat.sphere$SE[match(uORF_stats_sphere$gene_id, raw.dat.sphere$X.gene_id)]
uORF_stats_sphere$CDS_RFP <- 0;
uORF_stats_sphere$CDS_RFP <- raw.dat.sphere$complete_CDS_RFP_FPKM[match(uORF_stats_sphere$gene_id, raw.dat.sphere$X.gene_id)]
uORF_stats_sphere$CDS_RNA <- 0;
uORF_stats_sphere$CDS_RNA <- raw.dat.sphere$complete_CDS_totalRNA_FPKM[match(uORF_stats_sphere$gene_id, raw.dat.sphere$X.gene_id)]
uORF_stats_sphere$leader_SSU <- 0;
uORF_stats_sphere$leader_SSU <- raw.dat.sphere$complete_leader_SSU_FPKM[match(uORF_stats_sphere$gene_id, raw.dat.sphere$X.gene_id)]
uORF_stats_sphere$leader_length <- 0;
uORF_stats_sphere$leader_length <- raw.dat.sphere$leader_length[match(uORF_stats_sphere$gene_id, raw.dat.sphere$X.gene_id)]

uORF_stats_shield$CDS_TE <- 0;
uORF_stats_shield$CDS_TE <- raw.dat.shield$TE[match(uORF_stats_shield$gene_id, raw.dat.shield$X.gene_id)]
uORF_stats_shield$CDS_SE <- 0;
uORF_stats_shield$CDS_SE <- raw.dat.shield$SE[match(uORF_stats_shield$gene_id, raw.dat.shield$X.gene_id)]
uORF_stats_shield$CDS_RFP <- 0;
uORF_stats_shield$CDS_RFP <- raw.dat.shield$complete_CDS_RFP_FPKM[match(uORF_stats_shield$gene_id, raw.dat.shield$X.gene_id)]
uORF_stats_shield$CDS_RNA <- 0;
uORF_stats_shield$CDS_RNA <- raw.dat.shield$complete_CDS_totalRNA_FPKM[match(uORF_stats_shield$gene_id, raw.dat.shield$X.gene_id)]
uORF_stats_shield$leader_SSU <- 0;
uORF_stats_shield$leader_SSU <- raw.dat.shield$complete_leader_SSU_FPKM[match(uORF_stats_shield$gene_id, raw.dat.shield$X.gene_id)]
uORF_stats_shield$leader_length <- 0;
uORF_stats_shield$leader_length <- raw.dat.shield$leader_length[match(uORF_stats_shield$gene_id, raw.dat.shield$X.gene_id)]

uORF_stats_64cell$stage <- "64cell"
uORF_stats_sphere$stage <- "Sphere"
uORF_stats_shield$stage <- "Shield"

uorf_stats_combined <- rbind(uORF_stats_64cell, uORF_stats_sphere, uORF_stats_shield)

# Test load location

uorf_stats_combined <- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/combined/uORF_up_down_stop_codon/uorfs.csv", header = T, data.table = F)

uorf_stats_combined <- uorf_stats_combined %>% mutate(upstream_SSU   =( leader_SSU_count_upstream / uorf_distance_to_TSS) )
uorf_stats_combined <- uorf_stats_combined %>% mutate(upstream_LSU   =( leader_LSU_count_upstream / uorf_distance_to_TSS) )
uorf_stats_combined <- uorf_stats_combined %>% mutate(downstream_SSU =( leader_SSU_count_downstream / uorf_distance_to_TIS) )
uorf_stats_combined <- uorf_stats_combined %>% mutate(downstream_LSU =( leader_LSU_count_downstream / uorf_distance_to_TIS) )
uorf_stats_combined <- uorf_stats_combined %>% mutate(up_down_ratio_SSU =( upstream_SSU / downstream_SSU ) )
uorf_stats_combined <- uorf_stats_combined %>% mutate(up_down_ratio_LSU =( upstream_LSU / downstream_LSU ) )

#normalise downstream values by CDS RNA
uorf_stats_combined <- uorf_stats_combined %>% mutate(downstream_SSU_by_RNA =( (leader_SSU_count_downstream / CDS_RNA )/ uorf_distance_to_TIS) )
uorf_stats_combined <- uorf_stats_combined %>% mutate(downstream_LSU_by_RNA =( (leader_LSU_count_downstream / CDS_RNA )/ uorf_distance_to_TIS) )

uorf_stats_combined <- uorf_stats_combined %>% mutate(SSU_LSU_stop_ratio = (leader_SSU_stop_codon/leader_LSU_stop_codon))

#use ntile to bin mfe into 4 equal groups
uorf_stats_combined <- uorf_stats_combined %>% mutate(kozak_quartile = ntile(uorf_stats_combined$uorf_kozak, 4))

#combine the low + high quantiles to middle.
uorf_stats_combined[uorf_stats_combined$kozak_quartile == 1, ]$kozak_quartile <- "Low"
uorf_stats_combined[uorf_stats_combined$kozak_quartile == 2, ]$kozak_quartile <- "Mid"
uorf_stats_combined[uorf_stats_combined$kozak_quartile == 3, ]$kozak_quartile <- "Mid"
uorf_stats_combined[uorf_stats_combined$kozak_quartile == 4, ]$kozak_quartile <- "High"
uorf_stats_combined$kozak_quartile <- factor(uorf_stats_combined$kozak_quartile, levels=c("Low", "Mid", "High"))


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#A87225 Dark Orange  (C)
#F5B04E Orange  (T)
#FFC26B Light Orange
#4EB3F5 Light Blue (G)
#146DA8 Dark Blue (A)

pallet_5  <- c("#A87225", "#F5B04E", "#FFD291", "#4EB3F5", "#146DA8")
pallet_4n <- c("#A87225", "#146DA8", "#F5B04E", "#4EB3F5")  #mucleotides
pallet_4  <- c("#A87225", "#F5B04E", "#4EB3F5", "#146DA8")
pallet_3o <- c("#A87225", "#F5B04E", "#FFD291") #oranges
pallet_3  <- c("#A87225", "#FFD291", "#146DA8") #brown, orange, blue
pallet_2  <- c("#A87225", "#146DA8")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#filtering for all plots

start_codons=c("ATG", "CTG", "GTG")

#uorf_stats_combined.filtered <- filter(uorf_stats_combined, CDS_RNA >= 10, uorf_distance_to_TSS >= 50, uorf_stop_distance_to_TIS >= 50, gene_overlaps_another_gene == 0, leader_overlap_an_upstream_gene == 0, gene_overlaps_ncRNA == 0, uorf_start_codon %in% start_codons)
uorf_stats_combined.filtered <- filter(uorf_stats_combined, CDS_RNA >= 1, uorf_distance_to_TSS >= 50, uorf_stop_distance_to_TIS >= 50, gene_overlaps_another_gene == 0, leader_overlap_an_upstream_gene == 0, gene_overlaps_ncRNA == 0, uorf_start_codon %in% start_codons)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#VS TE

ATG_TE_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(CDS_TE), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  coord_cartesian(ylim = c(-5.5, 4)) +
  scale_fill_manual(values = rev(pallet_3o)) +
  ylab("log2(Translational efficiency") +
  theme_bw() +
  theme(legend.position="none") +
  NULL

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#VS Downstrem LSU

ATG_LSU_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(downstream_LSU), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  coord_cartesian(ylim = c(-7, 9)) +
  scale_fill_manual(values = rev(pallet_3o)) +
  ylab("log2(Downstream 80S)") +
  xlab("uORF stop codon") +
  theme_bw() +
  theme(legend.position="none") +
  NULL

out.supp <- grid.arrange(ATG_LSU_vs_stop_codon,ATG_TE_vs_stop_codon, ncol=2)
#ggsave("/Home/ii/adamg/analysis/RCP/uORF_plots/uORF_ATG_CTG_GTC_TE_and_Downstream.png", out.supp, height=150, width=100, units = 'mm')

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#VS Downstrem LSU by RNA

ATG_LSU_by_rna_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(downstream_LSU_by_RNA), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  coord_cartesian(ylim = c(-11, 5)) +
  scale_fill_manual(values = rev(pallet_3o)) +
  ylab("log2(Downstream 80S)") +
  xlab("uORF stop codon") +
  theme_bw() +
  theme(legend.position="none") +
  NULL

out.supp2 <- grid.arrange(ATG_LSU_by_rna_vs_stop_codon, ATG_TE_vs_stop_codon, ncol=2)
ggsave("/Users/adamgiess/Desktop/RCP/uORF_plots/june_files/uORF_ATG_CTG_GTC_TE_and_Downstream_by_RNA.png", out.supp2, height=150, width=100, units = 'mm')

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# ratio of 40S/80 stop codon coverage

stop_SSU_LSU_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(SSU_LSU_stop_ratio), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  scale_fill_manual(values = rev(pallet_3o)) +
  coord_cartesian(ylim = c(-9, 7)) +
  ylab("log2(40S / 80S stop ratio)") +
  xlab("uORF stop codon") +
  theme_bw() +
  theme(legend.position="none") +
  NULL

ggsave("/Users/adamgiess/Desktop/RCP/uORF_plots/june_files/uORF_ATG_CTG_GTC_stop_vs_SSU_by_LSU_over_stop.png", stop_SSU_LSU_vs_stop_codon, height=150, width=50, units = 'mm')
ggsave("/Users/adamgiess/Desktop/RCP/uORF_plots/june_files/uORF_ATG_CTG_GTC_stop_vs_SSU_by_LSU_over_stop.pdf", stop_SSU_LSU_vs_stop_codon, height=150, width=50, units = 'mm')


out.supp3 <- grid.arrange(ATG_LSU_by_rna_vs_stop_codon, stop_SSU_LSU_vs_stop_codon, ATG_TE_vs_stop_codon, ncol=3)
ggsave("/Users/adamgiess/Desktop/RCP/uORF_plots/june_files/uORF_ATG_CTG_GTC_stop_TE_and_Downstream_by_RNA.png", out.supp3, height=150, width=150, units = 'mm')
ggsave("/Users/adamgiess/Desktop/RCP/uORF_plots/june_files/uORF_ATG_CTG_GTC_stop_TE_and_Downstream_by_RNA.pdf", out.supp3, height=150, width=150, units = 'mm')



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Stats 1FPKM
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

#VS Downstream LSU by RNA

median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$CDS_TE)  #*0.5586155*  #*-0.8400726*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAG")$CDS_TE)  #*0.5598959*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$CDS_TE)  #*0.5394704*  #*-0.8903843*

log2(median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$CDS_TE)) #*-0.8400726*
log2(median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$CDS_TE)) #*-0.8903843*

median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA", CDS_TE > 0)$CDS_TE, na.rm = T)  # *0.5965621*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAG", CDS_TE > 0)$CDS_TE, na.rm = T)  # *0.5976088*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA", CDS_TE > 0)$CDS_TE, na.rm = T)  # *0.5677746*

median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$downstream_LSU_by_RNA)  #*0.06463817* == -4.001
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAG")$downstream_LSU_by_RNA)  #*0.06595991*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$downstream_LSU_by_RNA)  #*0.07483407*

log2(median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$downstream_LSU_by_RNA)) # *-3.95147*
log2(median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$downstream_LSU_by_RNA)) # *-3.740161*

median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA", downstream_LSU_by_RNA > 0)$downstream_LSU_by_RNA, na.rm = T)  #*0.1336634*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAG", downstream_LSU_by_RNA > 0)$downstream_LSU_by_RNA, na.rm = T)  #*0.127952*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA", downstream_LSU_by_RNA > 0)$downstream_LSU_by_RNA, na.rm = T)  #*0.1473359*

median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$SSU_LSU_stop_ratio, na.rm = T)  #*0.22*
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAG")$SSU_LSU_stop_ratio, na.rm = T)  #0.1111111
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$SSU_LSU_stop_ratio, na.rm = T)  #*0.05263158*

log2(median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$SSU_LSU_stop_ratio, na.rm = T)) # *-2.184425*
log2(median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$SSU_LSU_stop_ratio, na.rm = T)) # *-4.247928*

median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA", SSU_LSU_stop_ratio > 0)$SSU_LSU_stop_ratio, na.rm = T)  #2.11
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAG", SSU_LSU_stop_ratio > 0)$SSU_LSU_stop_ratio, na.rm = T)  #1.25
median(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA", SSU_LSU_stop_ratio > 0)$SSU_LSU_stop_ratio, na.rm = T)  #1



#Wilcoxon rank sum test with continuity correction
wilcox.test(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$downstream_LSU_by_RNA, filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$downstream_LSU_by_RNA, alternative="two.sided") #Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
#W = *1190600000*, p-value < 2.2e-16

wilcox.test(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$CDS_TE, filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$CDS_TE, alternative="two.sided") #Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
#W = 1298200000, p-value = 6.128e-11

wilcox.test(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$SSU_LSU_stop_ratio, filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$SSU_LSU_stop_ratio, alternative="two.sided") #Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
#W = 250820000, p-value < 2.2e-16
