#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO (this is not used in article anymore)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Used to make uORF plots for article

library(ORFikPipeline)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

# Will load data if you have run this before
makeCombinedUORFs <- function(path = "/export/valenfs/projects/Hakon/RCP_SEQ/uorf_stats_combined.csv") {
  if (file.exists(path))
    return(fread(path))

  raw.dat.64cell<-fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/matrix_files/64cell_WT_reps_1_2_peaks_removed_translating_filter.csv", header = TRUE, sep = ",")
  raw.dat.sphere<-fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/matrix_files/sphere_WT_reps_1_2_3_peaks_removed_translating_filter.csv", header = TRUE, sep = ",")
  raw.dat.shield<-fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/matrix_files/shield_WT_reps_1_2_3_peaks_removed_translating_filter.csv", header = TRUE, sep = ",")

  raw.dat.64cell <- raw.dat.64cell %>% mutate(TE = complete_CDS_RFP_FPKM / complete_CDS_totalRNA_FPKM)
  raw.dat.sphere <- raw.dat.sphere %>% mutate(TE = complete_CDS_RFP_FPKM / complete_CDS_totalRNA_FPKM)
  raw.dat.shield <- raw.dat.shield %>% mutate(TE = complete_CDS_RFP_FPKM / complete_CDS_totalRNA_FPKM)

  raw.dat.64cell <- raw.dat.64cell %>% mutate(SE = complete_leader_SSU_FPKM / complete_CDS_totalRNA_FPKM)
  raw.dat.sphere <- raw.dat.sphere %>% mutate(SE = complete_leader_SSU_FPKM / complete_CDS_totalRNA_FPKM)
  raw.dat.shield <- raw.dat.shield %>% mutate(SE = complete_leader_SSU_FPKM / complete_CDS_totalRNA_FPKM)

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  #uORF up dpown regions, with down regions from uORF stop to TIS

  #uORF_stats_64cell <- read.csv("~/analysis/TCP/uORF/64cell_uorf_near_canonical_up_down_counts_from_uORF_stop_codons.csv", header = T)
  #uORF_stats_sphere <- read.csv("~/analysis/TCP/uORF/sphere_uorf_near_canonical_up_down_counts_from_uORF_stop_codons.csv", header = T)
  #uORF_stats_shield <- read.csv("~/analysis/TCP/uORF/shield_uorf_near_canonical_up_down_counts_from_uORF_stop_codons.csv", header = T)

  uORF_stats_64cell <- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/uORF_matrix/uorfs_64cell_up_down_over_stop.csv.gz", header = T)
  uORF_stats_sphere <- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/uORF_matrix/uorfs_sphere_up_down_over_stop.csv.gz", header = T)
  uORF_stats_shield <- fread("/export/valenfs/projects/adam/TCP_seq/RCP_files/uORF_matrix/uorfs_shield_up_down_over_stop.csv.gz", header = T)

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  #Add in TE / RNA
  uORF_stats_64cell$CDS_TE <- 0;
  uORF_stats_64cell$CDS_TE <- raw.dat.64cell$TE[match(uORF_stats_64cell$gene_id, raw.dat.64cell$`#gene_id`)]
  uORF_stats_64cell$CDS_SE <- 0;
  uORF_stats_64cell$CDS_SE <- raw.dat.64cell$SE[match(uORF_stats_64cell$gene_id, raw.dat.64cell$`#gene_id`)]
  uORF_stats_64cell$CDS_RFP <- 0;
  uORF_stats_64cell$CDS_RFP <- raw.dat.64cell$complete_CDS_RFP_FPKM[match(uORF_stats_64cell$gene_id, raw.dat.64cell$`#gene_id`)]
  uORF_stats_64cell$CDS_RNA <- 0;
  uORF_stats_64cell$CDS_RNA <- raw.dat.64cell$complete_CDS_totalRNA_FPKM[match(uORF_stats_64cell$gene_id, raw.dat.64cell$`#gene_id`)]
  uORF_stats_64cell$leader_SSU <- 0;
  uORF_stats_64cell$leader_SSU <- raw.dat.64cell$complete_leader_SSU_FPKM[match(uORF_stats_64cell$gene_id, raw.dat.64cell$`#gene_id`)]
  uORF_stats_64cell$leader_length <- 0;
  uORF_stats_64cell$leader_length <- raw.dat.64cell$leader_length[match(uORF_stats_64cell$gene_id, raw.dat.64cell$`#gene_id`)]
  uORF_stats_64cell$leader_RFP_FPKM <- raw.dat.64cell$leader_RFP_FPKM[match(uORF_stats_64cell$gene_id, raw.dat.64cell$`#gene_id`)]

  uORF_stats_sphere$CDS_TE <- 0;
  uORF_stats_sphere$CDS_TE <- raw.dat.sphere$TE[match(uORF_stats_sphere$gene_id, raw.dat.sphere$`#gene_id`)]
  uORF_stats_sphere$CDS_SE <- 0;
  uORF_stats_sphere$CDS_SE <- raw.dat.sphere$SE[match(uORF_stats_sphere$gene_id, raw.dat.sphere$`#gene_id`)]
  uORF_stats_sphere$CDS_RFP <- 0;
  uORF_stats_sphere$CDS_RFP <- raw.dat.sphere$complete_CDS_RFP_FPKM[match(uORF_stats_sphere$gene_id, raw.dat.sphere$`#gene_id`)]
  uORF_stats_sphere$CDS_RNA <- 0;
  uORF_stats_sphere$CDS_RNA <- raw.dat.sphere$complete_CDS_totalRNA_FPKM[match(uORF_stats_sphere$gene_id, raw.dat.sphere$`#gene_id`)]
  uORF_stats_sphere$leader_SSU <- 0;
  uORF_stats_sphere$leader_SSU <- raw.dat.sphere$complete_leader_SSU_FPKM[match(uORF_stats_sphere$gene_id, raw.dat.sphere$`#gene_id`)]
  uORF_stats_sphere$leader_length <- 0;
  uORF_stats_sphere$leader_length <- raw.dat.sphere$leader_length[match(uORF_stats_sphere$gene_id, raw.dat.sphere$`#gene_id`)]
  uORF_stats_sphere$leader_RFP_FPKM <- raw.dat.sphere$leader_RFP_FPKM[match(uORF_stats_sphere$gene_id, raw.dat.sphere$`#gene_id`)]

  uORF_stats_shield$CDS_TE <- 0;
  uORF_stats_shield$CDS_TE <- raw.dat.shield$TE[match(uORF_stats_shield$gene_id, raw.dat.shield$`#gene_id`)]
  uORF_stats_shield$CDS_SE <- 0;
  uORF_stats_shield$CDS_SE <- raw.dat.shield$SE[match(uORF_stats_shield$gene_id, raw.dat.shield$`#gene_id`)]
  uORF_stats_shield$CDS_RFP <- 0;
  uORF_stats_shield$CDS_RFP <- raw.dat.shield$complete_CDS_RFP_FPKM[match(uORF_stats_shield$gene_id, raw.dat.shield$`#gene_id`)]
  uORF_stats_shield$CDS_RNA <- 0;
  uORF_stats_shield$CDS_RNA <- raw.dat.shield$complete_CDS_totalRNA_FPKM[match(uORF_stats_shield$gene_id, raw.dat.shield$`#gene_id`)]
  uORF_stats_shield$leader_SSU <- 0;
  uORF_stats_shield$leader_SSU <- raw.dat.shield$complete_leader_SSU_FPKM[match(uORF_stats_shield$gene_id, raw.dat.shield$`#gene_id`)]
  uORF_stats_shield$leader_length <- 0;
  uORF_stats_shield$leader_length <- raw.dat.shield$leader_length[match(uORF_stats_shield$gene_id, raw.dat.shield$`#gene_id`)]
  uORF_stats_shield$leader_RFP_FPKM <- raw.dat.shield$leader_RFP_FPKM[match(uORF_stats_shield$gene_id, raw.dat.shield$`#gene_id`)]

  uORF_stats_64cell$stage <- "64cell"
  uORF_stats_sphere$stage <- "Sphere"
  uORF_stats_shield$stage <- "Shield"

  uorf_stats_combined <- rbind(uORF_stats_64cell, uORF_stats_sphere, uORF_stats_shield)

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
  fwrite(uorf_stats_combined, file = path)
  return(uorf_stats_combined)
}

uorf_stats_combined <- makeCombinedUORFs()

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
#uorf_stats_combined.filtered <- filter(uorf_stats_combined, CDS_RNA >= 1, uorf_distance_to_TSS >= 50, uorf_stop_distance_to_TIS >= 50, gene_overlaps_another_gene == 0, leader_overlap_an_upstream_gene == 0, gene_overlaps_ncRNA == 0, uorf_start_codon %in% start_codons, leader_LSU_count_downstream > 10)
# New filter, close to TIS
# Used filter values, 4 and 10
uorf_stats_combined.filtered <- filter(uorf_stats_combined, CDS_RNA >= 1, uorf_stop_distance_to_TIS <= 50,
                                       gene_overlaps_another_gene == 0, leader_overlap_an_upstream_gene == 0,
                                       gene_overlaps_ncRNA == 0, uorf_start_codon %in% start_codons,
                                       leader_RFP_FPKM > -1)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#VS TE

ATG_TE_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(CDS_TE), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  coord_cartesian(ylim = c(-5.5, 4)) +
  scale_fill_manual(values = rev(pallet_3o)) +
  ylab("log2(CDS translational efficiency") +
  theme_bw() +
  theme(legend.position="none") +
  NULL

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#VS Downstrem LSU

ATG_LSU_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(downstream_LSU), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  coord_cartesian(ylim = c(-7, 9)) +
  scale_fill_manual(values = rev(pallet_3o)) +
  ylab("log2(downstream 80S)") +
  theme_bw() +
  theme(legend.position="none") +
  NULL

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#VS Downstrem LSU by RNA

ATG_LSU_by_rna_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(downstream_LSU_by_RNA), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  coord_cartesian(ylim = c(-11, 5)) +
  scale_fill_manual(values = rev(pallet_3o)) +
  ylab("log2(downstream 80S)") +
  theme_bw() +
  theme(legend.position="none") +
  NULL

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# ratio of 40S/80 stop codon coverage

stop_SSU_LSU_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(SSU_LSU_stop_ratio), fill=as.factor(uorf_stop_codon))) +
  geom_boxplot(width=0.6, outlier.colour=NA) +
  scale_fill_manual(values = rev(pallet_3o)) +
  coord_cartesian(ylim = c(-9, 7)) +
  theme_bw() +
  theme(legend.position="none") +
  NULL

out.supp3 <- grid.arrange(ATG_LSU_by_rna_vs_stop_codon, stop_SSU_LSU_vs_stop_codon, ATG_TE_vs_stop_codon, ncol=3)



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Stats 1FPKM
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#Wilcoxon rank sum test with continuity correction
wilcox.test(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$downstream_LSU_by_RNA, filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$downstream_LSU_by_RNA, alternative="two.sided") #Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
#W = 1228700000, p-value < 2.2e-16

wilcox.test(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$CDS_TE, filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$CDS_TE, alternative="two.sided") #Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
#W = 1298200000, p-value = 6.128e-11

wilcox.test(filter (uorf_stats_combined.filtered, uorf_stop_codon == "TAA")$SSU_LSU_stop_ratio, filter (uorf_stats_combined.filtered, uorf_stop_codon == "TGA")$SSU_LSU_stop_ratio, alternative="two.sided") #Wilcoxon rank sum test (equivalent to the Mann-Whitney test)
#W = 250820000, p-value < 2.2e-16

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# NEW VERSION WITH COVERAGE OF START REGION
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
df <- read.experimentl("Chew13")
df <- df[c(3,7,10),]
# New rna seq
dfr <- read.experimentl("Val19")
dfr <- dfr[c(2,3,6:11),]
fa <- ORFik:::findFa(df)
txdb <- loadTxdb(df, seqlevelsStyle(fa))
loadRegions(txdb, c("leaders", "cds"))
seqlevels(leaders)[26] <- "MT"

leaders <- reassignTSSbyCage(leaders, "/export/valenfs/data/processed_data/CAGE/nepal_2013_zebrafish/final_results/aligned_GRCz10/09_shield.1.bam", preCleanup = FALSE)
uorfs <- findUORFs(leaders, fa, startCodon = c("ATG|CTG|GTG"), stopCodon = c("TAA|TAG|TGA"), longestORF = FALSE)
seqlevels(uorfs)[26] <- "MT"
print(paste("Adam has: ", nrow(filter(uorf_stats_combined, uorf_start_codon %in% start_codons)), "uorfs"))
print(paste("I get: ", length(uorfs), "uorfs"))

outputLibs(list(df, dfr))

stopcodons <- as.character(ORFik:::txSeqsFromFa(stopCodons(uorfs, TRUE), fa))
dt <- data.table(txNames = txNames(uorfs), orfName = names(uorfs),uorfID = ORFik:::orfID(uorfs),
                 distToCds = distToCds(uorfs, leaders), distToTSS = distToTSS(uorfs, leaders),
                 uorf_stop_codon = stopcodons)
#cdsRNA <- countTable(dfr, "cds", "fpkm")
cdsRNA <- countTable(dfr, "cds", "fpkm", collapse = TRUE)[, c("RNA_64cell", "RNA_Sphere", "RNA_Shield")]
dt$cdsRNA <- rowSums(cdsRNA)[chmatch(dt$txNames, names(cds))]

cdsRFP <- makeSummarizedExperimentFromBam(df,region = cds, saveName = "/export/valenfs/projects/Hakon/AdamVienna/QC_STATS/countTable_cdsRFP.rds")
sums <- rowSums(assay(cdsRFP))
dt$cdsRFP <- ORFik:::fpkm_calc(sums, lengthSize = widthPerGroup(cds), librarySize = sum(sums))[chmatch(dt$txNames, names(cds))]

uorfRFP <- makeSummarizedExperimentFromBam(df, region = uorfs, saveName = "/export/valenfs/projects/Hakon/AdamVienna/QC_STATS/countTable_uorfs.rds")
sums <- rowSums(assay(uorfRFP))
dt$uorfRFP <- ORFik:::fpkm_calc(sums, lengthSize = widthPerGroup(uorfs), librarySize = sum(sums))
uorfStart <- startRegion(uorfs, leaders, upstream = 5, downstream = 5)
uorfstartRFP <- makeSummarizedExperimentFromBam(df, region = uorfStart, saveName = "/export/valenfs/projects/Hakon/AdamVienna/QC_STATS/countTable_uorfstartregion.rds")
sums <- rowSums(assay(uorfstartRFP))
dt$uorfstartRFP <- sums

dt$CDS_TE <- (dt$cdsRFP) / (dt$cdsRNA + 0.000001)
saveRDS(dt, "/export/valenfs/projects/Hakon/AdamVienna/new_TE_uorf_table.rds")
############################# STARTING POINT WHEN RUN ALREADY! ####################################################
dt <- readRDS("/export/valenfs/projects/Hakon/AdamVienna/new_TE_uorf_table.rds")
uorf_stats_combined <- dt
uorf_stats_combined$CDS_RNA <- uorf_stats_combined$cdsRNA
uorf_stats_combined$uorf_stop_distance_to_TIS <- uorf_stats_combined$distToCds

# New stats and plots
summary(dt)
############################################# NEW PLOTS ##################################################################
# With < 50 no other filter
uorf_stats_combined.filtered <- filter(uorf_stats_combined,  uorf_stop_distance_to_TIS <= 50)
dim(uorf_stats_combined.filtered)

dist <- c(50, 100)
rnafilter <- c(-1, 1)
uorfRFPfilt <- c(-1, 1)
uorfstartfilt <- c(-1, 2)
d <- list()
setting <- c()
nGenes <- c()
for (i in dist) {
  for (j in rnafilter) {
    for (k in uorfRFPfilt) {
      for (l in uorfstartfilt) {
        if (i == 50) {
          tab <- filter(uorf_stats_combined,  uorf_stop_distance_to_TIS <= 50, cdsRNA > j, uorfRFP > k, uorfstartRFP > l)
          nGenes <- c(nGenes, nrow(tab))
          d <- c(d, list(tab))
          setting <- c(setting, paste("stop_distance_TIS <= 50, cdsRNA >", j, ", uorfRFP >", k,", uorfstartRFP > ",l,")"))
        } else {
          tab <- filter(uorf_stats_combined,  uorf_stop_distance_to_TIS >= 100, cdsRNA > j, uorfRFP > k, uorfstartRFP > l)
          nGenes <- c(nGenes, nrow(tab))
          d <- c(d, list(filter(uorf_stats_combined,  uorf_stop_distance_to_TIS >= 100, cdsRNA > j, uorfRFP > k, uorfstartRFP > l)))
          setting <- c(setting, paste("stop_distance_TIS >= 100, cdsRNA >", j, ", uorfRFP >", k,", uorfstartRFP > ",l,")"))
        }
      }
    }
  }
}

#uorf_stats_combined.filtered <- filter(uorf_stats_combined,  uorf_stop_distance_to_TIS <= 50,
#                                       cdsRNA > 1, uorfRFP > 1, uorfstartRFP > 2)
gg <- list()
x = 1
for (uorf_stats_combined.filtered in d) {
  ATG_TE_vs_stop_codon <- ggplot(data=uorf_stats_combined.filtered, aes( x=uorf_stop_codon, y=log2(CDS_TE), fill=as.factor(uorf_stop_codon))) +
    geom_boxplot(width=0.6, outlier.colour=NA) +
    coord_cartesian(ylim = c(-5.5, 4)) +
    scale_fill_manual(values = rev(pallet_3o)) +
    ylab("log2(CDS translational efficiency") +
    theme_bw() +
    labs(title = setting[x], subtitle = paste("nGenes =", nGenes[x])) +
    theme(legend.position="none", plot.title = element_text(size = 9))
  gg <- c(gg, list(ATG_TE_vs_stop_codon))
  x = x + 1
}
do.call("grid.arrange", c(gg, ncol=4))

list50 <- setDT(d[[8]])
list100 <- setDT(d[[16]])

a <- list50[, .(mean = mean(CDS_TE), median = median(CDS_TE)), by = uorf_stop_codon]
a[,2:3] <- a[,2:3] - rbind(a[3,2:3], a[3,2:3], a[3,2:3])
print("< 50")
a
b <- list100[, .(mean = mean(CDS_TE), median = median(CDS_TE)), by = uorf_stop_codon]
b[,2:3] <- b[,2:3] - rbind(b[1,2:3], b[1,2:3], b[1,2:3])
print("> 100")
b
