# TODO: Option to remove 3' UTR tail, something is strange

#SSU = "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/combined_time_points/all_three_SSU_peaks_removed_200_removed_translating_lengths_25_35.bam"
#LSU = "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/combined_time_points/all_three_LSU_peaks_removed_200_selected_translating_lengths_25_35.bam"


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(gridExtra)

Palette1 <- c('skyblue4', 'orange')

outdir <- p(mainFolder, "/")
sample_name <- "tcp_plot_old"

#coverage plots: Leader, CDS, Traier plots scaled to 100nt
#coverage <- read.csv(args[3], header=T)
coverage <- read.csv("/export/valenfs/projects/adam/TCP_seq/valen_8/processed_28_09_18_most_expressed_transcript/scaled_feature_plots_translating_filter/SSU_peaks_removed_200_removed_translating_lengths_25_35_gene_matrix_coverage.csv", header=T)
coverage <- as.data.table(coverage)
#c <- copy(coverage)
coverage$Feature  <- factor(coverage$Feature, levels=c("leader", "cds", "trailer"), labels=c("leader", "CDS", "trailer"))
coverage$Fraction <- factor(coverage$Fraction, levels=c("SSU", "LSU"), labels = c("SSU", "LSU"))

# Keep only genes with coverage > 1 on SSU
coverage[, `:=` (gene_sum = sum(Count)), by = list(Gene, Fraction)] 
valid_genes <- coverage[Fraction == "SSU" & gene_sum > 0 ,]
unique_valid_genes <- unique(valid_genes$Gene)
coverage_sum_subset <- coverage[Gene %in% unique_valid_genes,]

coverage_sum_subset[, `:=` (windowSD=sd(Count), windowMean=mean(Count)),
                    by = list(Gene, Fraction)] 
coverage_sum_subset[,  zscore := (Count-windowMean)/windowSD]
coverage_zscore <- coverage_sum_subset[, .(zscore_sum=sum(zscore, na.rm = T), zscore_mean=mean(zscore, na.rm = T)),
                                       by = list(Fraction, Positon, Feature)]
coverage_zscore[, `:=` (fraction_min=min(zscore_mean)), by = Fraction]

plot_01 <- ggplot(data=as.data.frame(coverage_zscore), aes(x=Positon, ymax=zscore_mean, ymin=fraction_min, y=zscore_mean, colour = as.factor(Fraction))) +
  geom_ribbon(stat="identity", position = "identity", aes(fill= as.factor(Fraction), alpha=0.5)) +
  geom_line() +
  theme_bw() +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=Palette1) +
  scale_color_manual(values=Palette1) +
  ggtitle(paste0("0-10 Genes n=", length(unique(coverage_sum_subset$Gene))) ) +
  xlab("Scaled position in transcript") + ylab("Zscore mean over transcript") +
  theme(legend.position="none") +
  facet_grid(Fraction ~ Feature, scales = "free")

ggsave(paste(outdir, sample_name, "_ranked.pdf", sep=""), plot_01, width=200, height=150, unit="mm", dpi=100, limitsize = FALSE)
