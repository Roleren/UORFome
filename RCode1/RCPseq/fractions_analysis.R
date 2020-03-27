#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Extended Data Figure S14: Coverage metaplots of all RCP-seq fractions.


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Original version (coverage plot) (USED IN ARTICLE)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(ggridges)   #geom_density_ridges(aes(fill= as.factor(Treatment))) +

#make coverage plots: Leader, CDS, Traier plots scaled to 100nt
coverage_1   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction9_S5_R1_001_gene_matrix_coverage.csv", header=T))
coverage_2   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction10_S7_R1_001_gene_matrix_coverage.csv", header=T))
coverage_3   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction12_S9_R1_001_gene_matrix_coverage.csv", header=T))
coverage_4   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction13_S11_R1_001_gene_matrix_coverage.csv", header=T))
coverage_5   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction14_S1_R1_001_gene_matrix_coverage.csv", header=T))
coverage_6   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction15_S2_R1_001_gene_matrix_coverage.csv", header=T))
coverage_7   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction16_S3_R1_001_gene_matrix_coverage.csv", header=T))
coverage_8   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction17_S4_R1_001_gene_matrix_coverage.csv", header=T))
coverage_9   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction18_S6_R1_001_gene_matrix_coverage.csv", header=T))
coverage_10   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction19_S8_R1_001_gene_matrix_coverage.csv", header=T))
coverage_11   <- (fread("/export/valenfs/projects/adam/TCP_seq/valen_5/scaled_plots_per_fraction_tRNA_multi_marking_overlaps/Fraction20_S10_R1_001_gene_matrix_coverage.csv", header=T))

coverage_1$Fraction   <- "09"
coverage_2$Fraction   <- "10"
coverage_3$Fraction   <- "12"
coverage_4$Fraction   <- "13"
coverage_5$Fraction   <- "14"
coverage_6$Fraction   <- "15"
coverage_7$Fraction   <- "16"
coverage_8$Fraction   <- "17"
coverage_9$Fraction   <- "18"
coverage_10$Fraction  <- "19"
coverage_11$Fraction  <- "20"

coverage <- rbindlist(list(coverage_1, coverage_2, coverage_3, coverage_4, coverage_5, coverage_6, coverage_7, coverage_8, coverage_9, coverage_10,  coverage_11))
coverage$Feature  <- factor(coverage$Feature, levels=c("leader", "cds", "trailer"), labels=c("leader", "CDS", "trailer"))

a <- copy(coverage)
coverage <- a
#remove genes with problematic repetative peaks
transcripts_to_exclude <- c("ENSDARG00000077330","ENSDARG00000102873","ENSDARG00000089382")
coverage <- subset(coverage, !Gene %in% transcripts_to_exclude)

coverage[, median := median(Count), by = .(Fraction, Feature, Gene)]
hits <- coverage$Count > coverage$median*200 + 25
genes <- unique(coverage$Gene[hits])

transcripts_to_exclude <- genes
coverage <- subset(coverage, !Gene %in% transcripts_to_exclude)

coverage.sum.subset <- coverage %>% group_by(Gene, Fraction) %>% filter(sum(Count) > 0) #remove genes without counts

# Raw #############################################################################################

#sum count per position for each gene
coverage.sum.subset.raw.sum <- coverage.sum.subset               %>% ungroup() %>% group_by(Fraction, Positon, Feature) %>% summarise(count_sum=sum(Count))

# Variable max y
plot_02 <- ggplot(data=as.data.frame(coverage.sum.subset.raw.sum ), aes(x=Positon, ymax = count_sum , ymin = 0, y=count_sum, colour = as.factor(Fraction))) +
  geom_line(aes(colour = as.factor(Fraction ))) +
  geom_ribbon(stat="identity", position = "identity", aes(fill= as.factor(Fraction), alpha=0.2)) +
  theme_bw() +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Position in scaled transcript region") + ylab("Mean of standardised counts") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(paste0("0-10 Genes n=", length(unique(coverage.sum.subset$Gene))) ) +
  scale_colour_viridis(option = "cividis", discrete = T) +
  scale_fill_viridis(option = "cividis", discrete = T) +
  facet_grid(Fraction ~ Feature, scales = "free")

ggsave(paste0(outFolder,"run5_transcript_coverage_metaplots_frations_raw_filtered_freeScale.png"), plot_02, width=300, height=500, unit="mm", dpi=300)
ggsave(paste0(outFolder,"run5_transcript_coverage_metaplots_frations_raw_filtered_freeScale.pdf"), plot_02, width=300, height=500, unit="mm", dpi=300)

# Equal max y
plot_02 <- ggplot(data=as.data.frame(coverage.sum.subset.raw.sum ), aes(x=Positon, ymax=count_sum, ymin=0, y=0, colour = as.factor(Fraction))) +
  geom_line(aes(colour = as.factor(Fraction ))) +
  geom_ribbon(stat="identity", position = "identity", aes(fill= as.factor(Fraction), alpha=0.2, ymin = 0)) +
  theme_bw() +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Position in scaled transcript region") + ylab("Mean of standardised counts") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle(paste0("0-10 Genes n=", length(unique(coverage.sum.subset$Gene))) ) +
  scale_colour_viridis(option = "cividis", discrete = T) +
  scale_fill_viridis(option = "cividis", discrete = T) +
  facet_grid(Fraction ~ Feature)

ggsave(paste0(outFolder,"run5_transcript_coverage_metaplots_frations_raw_filtered.png"), plot_02, width=300, height=500, unit="mm", dpi=100)



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# DATA (my version, only for sanity test)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load Experiment
library(ORFikPipeline)
dfR <- read.experimentl("Valen5_2018") # Load fractions library
dfR@expInVarName <- FALSE

# Annotation
txdb <- loadTxdb(dfR@txdb)
txNames <- filterTranscripts(txdb, 100, 100, 100, T) # Filter
loadRegions(txdb, names.keep = txNames)

# Output data
outFolder <- "/export/valenfs/projects/Hakon/RCP_SEQ/plots/"
outputLibs(dfR)

a <- transcriptWindow(leaders, cds, trailers, df = dfR, outdir = outFolder, allTogether = TRUE,
                      colors = "blue", returnPlot = T)
