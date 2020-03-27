s#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Here are some of the TE AND SE box plot figures in RCP-seq article
# Also some uORF stuff
library(ORFikPipeline)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(cowplot)

###################################################################################################
# Colours
###################################################################################################
new_pallet_4 <- c("#E495A5","#ABB065","#39BEB1","#ACA4E2")
new_pallet_3 <- c("#39BEB1", "#ACA4E2", "#E495A5")
tl <- theme(legend.position = c(0.2, 0.67), legend.background=element_blank())
tlb <- theme(legend.position = c(0.8, 0.3), legend.background=element_blank(), axis.ticks.y=element_blank(), axis.text.y = element_blank())
tit <- labs(color = "First nucleotide")
titb <- labs(color = "Motif")

# Figure 2 D & E
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# SE
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# New RNA seq
combined.dat.ranked <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq2020feb.rds")
combined.dat.ranked$complete_CDS_totalRNA_FPKM_new <- combined.dat.ranked$complete_CDS_totalRNA_FPKM_new2

# Add top
combined.dat.ranked$TOP_plusC <- "Not TOP other"
combined.dat.ranked$TOP_plusC[ combined.dat.ranked$TSS_start_sequence == "C" ] <- "No TOP C"
combined.dat.ranked$TOP_plusC[ combined.dat.ranked$TOP == TRUE ] <- "TOP"
# Filter
combined.dat.ranked.filter <- combined.dat.ranked %>% filter(complete_CDS_totalRNA_FPKM_new >=10, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, TSS_start_sequence!="N", histone==FALSE)
table(combined.dat.ranked.filter$TSS_start_sequence)
table(combined.dat.ranked.filter$TOP_plusC)


data.subset.by.rna <-
  combined.dat.ranked.filter %>%
  select (
    X.gene_id,
    stage,
    TSS_start_sequence,
    Scanning_efficiency
  )
data.subset.by.rna.melt<- melt(data.subset.by.rna, id.vars = c("X.gene_id", "stage", "TSS_start_sequence"))
se1 <- ggplot(data=data.subset.by.rna.melt, aes((value), colour = TSS_start_sequence)) +
  stat_ecdf() +
  scale_x_log10() +
  coord_cartesian(xlim=c(0.1,10)) +
  scale_color_manual(values=new_pallet_4) +
  ylab("") +
  xlab("") +
  theme_bw() +
  tl + tit
ste1 <- copy(data.subset.by.rna.melt)
#------#
# TOP Motif
#------#


data.subset.by.rna<-
  combined.dat.ranked.filter %>%
  select (
    X.gene_id,
    stage,
    TOP_plusC,
    Scanning_efficiency
  )
data.subset.by.rna.melt<- melt(data.subset.by.rna, id.vars = c("X.gene_id", "stage", "TOP_plusC"))

#se2 <- ggplot(data=data.subset.by.rna.melt, aes((value), colour = TOP)) +
se2 <- ggplot(data=data.subset.by.rna.melt, aes((value), colour = TOP_plusC)) +
  stat_ecdf() +
  scale_x_log10() +
  coord_cartesian(xlim=c(0.1,10)) +
  scale_color_manual(values=new_pallet_3) +
  ylab("") +
  xlab("") +
  theme_bw() +
  tlb + titb

se1 <- se1 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
se2 <- se2 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Statistics
ste2 <- copy(data.subset.by.rna.melt)
notopc <- ste2[ste2$TOP_plusC == "No TOP C",]
notopother <- ste2[ste2$TOP_plusC == "Not TOP other",]

wilcox.test(notopc$value, notopother$value)
###################################################################################################
# TE (LSU)
###################################################################################################
combined.dat.ranked.filter <- combined.dat.ranked %>% filter(complete_CDS_totalRNA_FPKM_new >=10, overlapping_gene==FALSE, leader_potentially_overlaps_upstream_gene==FALSE, gene_potentially_overlaps_downstream_leader==FALSE, gene_overlaps_non_coding_transcript==FALSE, TSS_start_sequence!="N", histone==FALSE)
combined.dat.ranked.filter$TE_LSU <- combined.dat.ranked.filter$CDS_LSU_FPKM / combined.dat.ranked.filter$complete_CDS_totalRNA_FPKM_new

data.subset.by.rna<-
  combined.dat.ranked.filter %>%
  select (
    X.gene_id,
    stage,
    TSS_start_sequence,
    TE_LSU
  )

data.subset.by.rna.melt<- melt(data.subset.by.rna, id.vars = c("X.gene_id", "stage", "TSS_start_sequence"))

te1 <- ggplot(data=data.subset.by.rna.melt, aes((value), colour = TSS_start_sequence)) +
  stat_ecdf() +
  scale_x_log10() +
  coord_cartesian(xlim=c(0.1,10)) +
  scale_color_manual(values=new_pallet_4) +
  ylab("") +
  xlab("") +
  theme_bw() +
  tl + tit

st1 <- copy(data.subset.by.rna.melt)

#------#
# TOP Motif
#------#

data.subset.by.rna<-
  combined.dat.ranked.filter %>%
  select (
    X.gene_id,
    stage,
    TOP_plusC,
    TE_LSU
  )

data.subset.by.rna.melt<- melt(data.subset.by.rna, id.vars = c("X.gene_id", "stage", "TOP_plusC"))

te2 <- ggplot(data=data.subset.by.rna.melt, aes((value), colour = TOP_plusC)) +
  stat_ecdf() +
  scale_x_log10() +
  coord_cartesian(xlim=c(0.1,10)) +
  scale_color_manual(values=new_pallet_3) +
  ylab("") +
  xlab("") +
  theme_bw() +
  tlb + titb

te1 <- te1 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
te2 <- te2 +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))


combse <- grid.arrange(plot_grid(se1, se2, nrow = 1, align = 'hv', axis  = 'tb'), bottom = "log10(Scanning efficiency)")
combte <- grid.arrange(plot_grid(te1, te2, nrow = 1, align = 'hv', axis  = 'tb'), bottom = "log10(Translational efficiency)")
comb <- grid.arrange(combse, combte, nrow = 2, left = "Cumulative fraction", widths = 20)

ggslackR(plot = comb, width = 200, height = 200)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# STATISTICS
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#################### No TOP C vs Not TOP other
st <- copy(data.subset.by.rna.melt)
notopc <- st[st$TOP_plusC == "No TOP C",]
notopother <- st[st$TOP_plusC == "Not TOP other",]
wilcox.test(notopc$value, notopother$value)

st <- copy(data.subset.by.rna.melt)
top <- st[st$TOP_plusC == "TOP",]
notopother <- st[st$TOP_plusC == "Not TOP other",]
wilcox.test(top$value, notopother$value)

##################### C vs others
# SE
C <- ste1[ste1$TSS_start_sequence == "C",]
Cother <- ste1[ste1$TSS_start_sequence %in% c("A", "G"),]
wilcox.test(C$value, Cother$value)
# TE
C <- st1[st1$TSS_start_sequence == "C",]
Cother <- st1[st1$TSS_start_sequence %in% c("A", "G"),]
wilcox.test(C$value, Cother$value)
# C & T vs other
# SE
Tt <- ste1[ste1$TSS_start_sequence %in% c("T"),]
Tother <- ste1[ste1$TSS_start_sequence %in% c("A", "G"),]
wilcox.test(Tt$value, Tother$value)
# TE
Tt <- st1[st1$TSS_start_sequence %in% c("T"),]
Tother <- st1[st1$TSS_start_sequence %in% c("A", "G"),]
wilcox.test(Tt$value, Tother$value)
