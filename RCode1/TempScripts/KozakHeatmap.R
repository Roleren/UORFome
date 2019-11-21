# Kozak heatmap
###################################################################################################
# Initiation rate by nucleotide and position
###################################################################################################

source("/export/valenfs/projects/uORFome/RCode1/loadUorfome.R")
library(ggpubr)
library(ggthemes)
plotFolder <- "/export/valenfs/projects/Håkon/AdamVienna/plots/new_plots/"
setwd("/export/valenfs/projects/Håkon/AdamVienna/")



ready <- readRDS("expression_both_withGorilla_filteredRFPSSU.rds")
# ER group
ready_ER <- ready[fraction == "GO_ER",] 
dim(ready_ER)
ready_ER_filtered <- ready_ER[complete_CDS_totalRNA_FPKM >=10& complete_leader_SSU_FPKM >=0& complete_CDS_RFP_FPKM >=0& overlapping_gene==FALSE& leader_potentially_overlaps_upstream_gene==FALSE& gene_potentially_overlaps_downstream_leader==FALSE& gene_overlaps_non_coding_transcript==FALSE& histone==FALSE& Initiation_rate_RPF > 0& !is.na(Initiation_rate_RPF)& is.finite(Initiation_rate_RPF)& leader_length >= 100, ]


#' Make region heatmap relative to scoring
#'
#' Given sequences, DNA or RNA.
#' And some score, ribo-seq fpkm, TE etc.
#' Create a heatmap divided per letter in seqs, by how strong the score is.
#'
#' @param seqs the sequences (character vector, DNAStringSet)
#' @param rate a scoring vector (equal size to seqs)
#' @param start position in seqs to start at (first is 1)
#' @param stop position in seqs to stop at (first is 1)
#' @param center position in seqs to center at (first is 1), center will
#' be +1 in heatmap
#' @param min.observations How many observations per position per letter to accept?
#' numeric or quantile, default (">q1", bigger than quartile 1 (25 percentile)).
#' You can do (10), to get all with more than 10 observations.
#' @param skip.startCodon startCodon is defined as after centering (position 1, 2 and 3).
#' Should they be skipped ? default (FALSE).
#' Not relevant if you are not doing Translation initiation sites (TIS).
#' @param xlab Region you are checking, default (TIS)
#' @param type What type is the rate scoring ? default (ribo-seq)
#' @return a ggplot of the heatmap
#' @importFrom data.table melt
#' @export
#' @examples

#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
#'   txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                           package = "GenomicFeatures")
#'   #Extract sequences of Coding sequences.
#'   cds <- loadRegion(txdbFile, "cds")
#'   tx <- loadRegion(txdbFile, "mrna")
#'   #Either you import fasta file of ranges, or you have some BSgenome.
#'
#'   bam_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#'   rfp <- readBam(bam_file)
#'   fpkm <- fpkm(cds, rfp)
#'
#'   # Get region to check
#'   kozakRegions <- startRegionString(cds, tx, BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#'                                     , upstream = 4, 5)
#'   kozakHeatmap(kozakRegions, fpkm, 1, 9)
#'
#'
kozakHeatmap <- function(seqs, rate, start, stop, center = ceiling((stop - start + 1)/2),
                         min.observations = ">q1", skip.startCodon = FALSE,
                         xlab = "TIS", type = "ribo-seq") {
  if (length(seqs) != length(rate)) stop("Length of rate and seq must be equal!")
  if (length(seqs) == 0 | length(rate) == 0) stop("Length of rate and seq must be > 0!")
  if (is.null(names(seqs))) names(seqs) <- as.character(seq(1, length(seqs)))
  
  dt <- data.table(X.gene_id = names(seqs))
  vars <- c()
  for (i in seq(start, stop)) {
    dt[,paste0("seq", i)] <- substring(seqs, i, i)
    vars <- c(vars, paste0("seq", i))
  }
  dt$rate <- rate
  dt.melt <- melt(dt, id.vars = c("X.gene_id", "rate"))
  
  # codon.table <- dt.melt %>% group_by(variable, value) %>% summarise(
  #   median_IR=median(rate, na.rm=T),
  #   count_seq_pos_with_count= n()
  # )
  codon.table <- dt.melt[, .(median_score = median(rate, na.rm=T),
                             count_seq_pos_with_count = .N),
                         by = .(variable, value)]
  
  uniques <- rev(unique(dt.melt$value))
  codon.table$value <- factor(codon.table$value, uniques)
  xPos <- seq(start, stop + 1) - start + 1 - center
  xPos <- xPos[-center]
  xPos[xPos > 0] <- paste0("+", xPos[xPos > 0])
  codon.table$variable <- factor(codon.table$variable, levels=vars,
                                 labels=xPos)
  
  if (skip.startCodon) {
    codon.table[codon.table$variable %in% c("+1", "+2", "+3"), ]$median_score <- NA
  }
  
  #output position where we have more than 1st quartile observations etc
  print("Distribution of observations per position per letter")
  print(summary(codon.table$count_seq_pos_with_count))
  if (min.observations == ">q1") {
    quart <- summary(codon.table$count_seq_pos_with_count)[2]
  } else if (is.numeric(min.observations)) {
    quart <- min.observations
  } else stop("min.observations must be >q1 or a numeric!")
  
  print(paste0("Picking: >", quart, " observations"))
  codon.table.filtered <- codon.table[count_seq_pos_with_count > quart,]
  if (nrow(codon.table.filtered) == 0) stop("No rows passed min.observations!")
  
  plot_matrix2_log <- ggplot(data=codon.table.filtered,
                             aes(x=variable, y=value, fill=log2(median_score))) +
    theme(panel.background=element_rect(fill="lightgrey", colour="lightgrey")) +
    geom_tile(color = "lightgrey") +
    scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = 'lightgrey',
                         name = paste0("log2(median_", type,")")) +
    xlab(paste0("Position realitive to ", xlab)) +
    ylab("Nucleotide")
  
  rows <- seq(from = length(unique(codon.table.filtered$value)) - 0.5, to = 0.5)
  names(rows) <- unique(codon.table.filtered$value)
  xmin <- -0.5
  for(col in unique(codon.table.filtered$variable)){
    mat <- codon.table.filtered[codon.table.filtered$variable == col,]
    highest <- rows[names(rows) == mat$value[which.max(mat$median_score)]]
    xmin = xmin + 1;xmax = xmin + 1; ymin = highest;ymax = ymin + 1
    if (length(highest) == 0) next
    input <- paste0("geom_rect(aes(xmin=",xmin,",xmax=",xmax,",ymin=",
                    ymin,",ymax=",ymax,"), color='black', size=0.5, fill=NA)")
    plot_matrix2_log <- plot_matrix2_log +
      eval(parse(text = input))
  }
  return(plot_matrix2_log)
}

seqs <- ready_ER_filtered$initiation_sequence
rate <- ready_ER_filtered$Initiation_rate_RPF
start <- 6; stop <- 14;  center <- ceiling((stop - start + 1)/2)
min.observations <- ">q1"; skip.startCodon = T; type = "IR"
plot_matrix2_log <- kozakHeatmap(seqs, rate, 
                                 start, stop, center, 
                                 min.observations, skip.startCodon, type = type)
plot_matrix2_log
ggsave(paste0(plotFolder, "IR_by_nucleotide_10FPKM_100nt_leaders_count1000_median.png"), plot_matrix2_log, height=100, width=250, units = 'mm', limitsize = F)


# Rerun from /export/valenfs/projects/adam/TCP_seq/RCP_files/RCP_plots/initiation_plotsinitiation_sequences_4_upstream_context.R
ready_adam <- setDT(readRDS("kozakTxWithGo.rds"))
dim(ready_adam)
ready_filtered <- ready_adam[complete_CDS_totalRNA_FPKM >=10& complete_leader_SSU_FPKM >=0& complete_CDS_RFP_FPKM >=0& overlapping_gene==FALSE& leader_potentially_overlaps_upstream_gene==FALSE& gene_potentially_overlaps_downstream_leader==FALSE& gene_overlaps_non_coding_transcript==FALSE& histone==FALSE& Initiation_rate_RPF > 0& !is.na(Initiation_rate_RPF)& is.finite(Initiation_rate_RPF)& leader_length >= 100, ]
plot_matrix2_log <- kozakHeatmap(ready_filtered$initiation_sequence, 
                                 ready_filtered$Initiation_rate_RPF,
                                 start, stop, center, 
                                 min.observations, skip.startCodon, type = type)
plot_matrix2_log
ggsave(paste0(plotFolder, "IR_original_by_nucleotide_10FPKM_100nt_leaders_count1000_median.png"), plot_matrix2_log, height=100, width=250, units = 'mm', limitsize = F)

# Boxplot difference in groups
# General
ready_filtered <-  cbind(ready_filtered, fraction = "Other")
ready_filtered[ready_filtered$transcript_id %in% ready_ER_filtered$transcript_id,]$fraction <- "ER"
ready_filtered <- ready_filtered[complete_leader_SSU_FPKM >0 & complete_CDS_RFP_FPKM > 0, ]
ggplot(data = ready_filtered, aes(x = fraction, y = Initiation_rate_RPF)) + 
  geom_boxplot() + 
  ylim(0, 20)

# Relative to kozak median
dt <- ready_filtered

dtt <- dt %>% group_by(initiation_sequence_sub, perfect_kozak, upstream_kozak_strength) %>% summarise(median_IR_RFP = median(Initiation_rate_RPF), mean_IR_RFP = mean(Initiation_rate_RPF), sequence_count = n(), median_TE=median(TE), IR_SD=sd(Initiation_rate_RPF))
tile_IE_median_20 <- ggplot(data=dtt, aes(x=initiation_sequence_sub, y=1, fill=log2(median_IR_RFP))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  #scale_fill_viridis(discrete=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NULL
tile_IE_median_20

dttt <- merge(dt, dtt, by = "initiation_sequence_sub")

dttt$IR_dif_ER_Other <- dttt$Initiation_rate_RPF - dttt$median_IR_RFP
# box
ggplot(data = dttt, aes(x = fraction, y = IR_dif_ER_Other)) + 
  geom_boxplot() + 
  ylim(-5, 5) + 
  ylab("IR: singles - median ")
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

# statistical test
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
