source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")
exp_dir = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/"
stats_folder <- paste0(exp_dir, "QC_STATS/")

df <- read.experimentl(exper_name); df@expInVarName <- FALSE

# Heatmaps
txNames <- filterTranscripts(df@txdb, 100, 100, 100, longestPerGene = FALSE)
txNames <- txNames[-which(startSites(leaders[txNames], F) < 50)] # Remove transcripts too close to scaffold starts
# TIS
hmFolder <- p(stats_folder, "heatmaps")
dir.create(hmFolder)
heatMapL(cds[txNames], mrna[txNames], df, outdir = hmFolder, upstream = 50, downstream = 49, shifting = "5prime", acLen = 20:70)
# TSS
heatMapL(leaders[txNames], extendLeaders(mrna[txNames], 50), df, outdir = p(exp_dir, "hm_"), upstream = 50, downstream = 49, shifting = c("5prime", "3prime"), acLen = 20:70, location = "TSS")

# Cor plot between experiments

dfV <- df
dfL <- read.experimentl("Lee13") 
dfl <- list(dfV, dfL); names(dfl) <- c("Val16", "Lee13")

data_for_pairs <- makeSummarizedExperimentFromBam(dfV, region = mrna, geneOrTxNames = "tx", longestPerGene = FALSE,
                                                  type = "fpkm")
data_for_pairs_Lee <- makeSummarizedExperimentFromBam(dfL, region = mrna, geneOrTxNames = "tx", longestPerGene = FALSE,
                                                  type = "fpkm")
data_for_pairs_both <- cbind(data_for_pairs, data_for_pairs_Lee)
paired_plot <- ggpairs(as.data.frame(data_for_pairs_both), columns = 1:ncol(data_for_pairs_both))
ggsave(paste0(stats_folder, "cor_plot_comparison_Lee2013.png"), paired_plot, height=500, width=500, units = 'mm', dpi=300)

# Log2
paired_plot_log2 <- ggpairs(as.data.frame(log2(data_for_pairs_both)), columns = 1:ncol(data_for_pairs_both))
ggsave(paste0(stats_folder, "cor_plot_comparison_Lee2013_log2.png"), paired_plot_log2, height=500, width=500, units = 'mm', dpi=300)
