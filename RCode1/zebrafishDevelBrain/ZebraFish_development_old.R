# DATA:
#brain_day_night
#brain_development
#development
#pre_fertilization

#TODO:
# 1. FPKM per library, average replicates
# 2. Table of increased time fpkms (dev, dev brain)
# 3. Day/night, find periodicity (switch on off day/night) ()
# 4. Intersection with the Eivind peptide list (dev brain, daynight brain to find neural peptides) (go from peptide name to tx in fasta)
# 5. Heatmaps /filter on fpkm, transcriptNormalized columns(time),
rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./TempScripts/tcp_pipeline.R")
source("./DataBaseSetup.R")
dbDisconnect(uorfDB)
rm(uorfDB)
library(ORFik)
#library(Rsubread)
library(DESeq2)
library(SummarizedExperiment)
library(pheatmap)
library(BiocParallel)
register(MulticoreParam(20))
setwd("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep")

gtfAnno <- "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81_chr.gtf"
gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
getFasta("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa")
txdb <- loadTxdb(gtfPath, seqlevelsStyle(fa));

tx <- loadRegion(txdb, "mrna")
names(tx) <- ORFik:::txNamesToGeneNames(names(tx), txdb)

# Protein
genes <- getPeptides()
write.table(genes, file = "GenesUsed/PeptideGenes.csv")


############################ brain_development ############################
libPath <- "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/brain_development/aligned"
libs <- allBamFilesInFolder(libPath)
df <- data.frame(RNA = libs, stage = c("10dfp", "10dfp", "10dfp", "21dfp", "21dfp", "21dfp", "6dfp", "6dfp", "6dfp"),
                 rep = rep(seq.int(3), 3), stringsAsFactors = F)
df <- df[c(7,8,9,1,2,3,4,5,6),]
df$stage <- factor(df$stage, levels = c("6dfp", "10dfp", "21dfp"))
df$rep <- factor(df$rep)

saveName <- "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/tables/experiments_brain_devel.rds"
#final <- makeSummarizedExperimentFromBam(df, txdb, saveName)
final <- readRDS(saveName)
dds <- scoreSummarizedExperiment(final)
fpkmCollapsed <- scoreSummarizedExperiment(final, "fpkm")



# New heatmap with only neuro peptides
isNeuroPep <- (rownames(dds) %in% genes) & (rowMax(fpkmCollapsed) > 1)
dff <- as.data.frame(colData(dds)[,c("SAMPLE","replicate")])
p <- pheatmap(assay(dds)[isNeuroPep,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=dff)
p
ggsave("brain_devel_heatmap.pdf", p, height = 30, limitsize = FALSE)
write.table(isNeuroPep, "GenesUsed/neuroPeptideGenes_brain_devel.csv")

############################ DAY / NIGHT ############################################
libPath2 <- "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/brain_day_night/aligned"
libs2 <- allBamFilesInFolder(libPath2)
dfT <- data.frame(RNA = libs2, stage = c(rep("6am", 8), rep("6pm", 8), rep("midnight", 8), rep("noon", 8)),
                 rep = rep(seq.int(8), 4))
dfT$RNA <- as.character(dfT$RNA)
dfT$stage <- factor(dfT$stage, levels = c("6am", "noon", "6pm", "midnight"))
saveName <- "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/tables/experiments_brain_day_night.rds"
#finalTime <- makeSummarizedExperimentFromBam(dfT, txdb, saveName)
finalTime <- readRDS(saveName)

dds <- scoreSummarizedExperiment(finalTime)
fpkmCollapsed <- scoreSummarizedExperiment(finalTime, "fpkm")
isNeuroPep <- (rownames(dds) %in% genes) & (rowMax(fpkmCollapsed) > 1)
dff <- as.data.frame(colData(dds)[,c("SAMPLE","replicate")])
p <- pheatmap(assay(dds)[isNeuroPep,], cluster_rows=TRUE, show_rownames=FALSE,
              cluster_cols=FALSE, annotation_col=dff)
p
ggsave("brain_day_night_heatmap.pdf", p, height = 30, limitsize = FALSE)
write.table(isNeuroPep, "GenesUsed/neuroPeptideGenes_day_night.csv")
de <- DESeqDataSet(finalTime, design = ~ SAMPLE)
dds <- DESeq(de)
res <- results(dds)

############################ Devel ############################
libPath3 <- "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/development/aligned"
libs3 <- allBamFilesInFolder(libPath3)
df <- data.frame(RNA = libs3, stage = c(rep("1KCell", 2), rep("2-4Cell", 2), rep("256Cell", 1),
                                        rep("28hpf", 2), rep("2dpf", 2), rep("5dpf", 2),
                                        rep("Bud", 2), rep("Bud_MZDicer", 1), rep("Bud_MZDicer_miR430", 1),
                                        rep("Dome", 2), rep("Dome_MZOep", 1), rep("Dome_Squint_1", 1),
                                        rep("Oblong", 1), rep("Shield", 3), rep("Shield_MZOep_1", 1),
                                        rep("Shield_Squint", 1)),
                 rep = c(seq.int(2), seq.int(2), seq.int(1), seq.int(2), seq.int(2), seq.int(2),
                         seq.int(2), seq.int(1), seq.int(1), seq.int(2), seq.int(1), seq.int(1),
                         seq.int(1),  seq.int(3), seq.int(1), seq.int(1)),
                 stringsAsFactors = FALSE)
df <- df[c(3,4,5,1,2,20, seq(16,19), seq(21,25), seq(12,15), seq(6,11)),]
df$stage <- factor(df$stage, levels = unique(df$stage))
df$rep <- factor(df$rep)
df <- df[-c(3,6,9,10),] # remove bad  TODO USE!
df <- df[-c(10, 14, 15),] # Remove mzdizer etc
df <- df[-c(10),] # remove squint

saveName <- "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/tables/experiments_development.rds"
#finalDevel <- makeSummarizedExperimentFromBam(df, txdb, saveName)
finalDevel <- readRDS(saveName)

dds <- scoreSummarizedExperiment(finalDevel)
fpkmCollapsed <- scoreSummarizedExperiment(finalDevel, "fpkm")

isNeuroPep <- (rownames(dds) %in% genes) & (rowMax(fpkmCollapsed) > 1)
dff <- as.data.frame(colData(dds)[,c("SAMPLE","replicate")])
p <- pheatmap(assay(dds)[isNeuroPep,], cluster_rows=TRUE, show_rownames=FALSE,
              cluster_cols=FALSE, annotation_col=dff)
p
ggsave("development_heatmap.pdf", p, height = 30, limitsize = FALSE)
write.table(isNeuroPep, "GenesUsed/neuroPeptideGenes_development.csv")

############################ Merged BRAIN ############################

merged <- cbind(assay(final), assay(finalTime))
colnames(merged) <- NULL

s <- c(as.character(df$stage), as.character(dfT$stage))
s <- factor(s, levels = unique(s))
colData <- DataFrame(SAMPLE=s,
                     replicate=c(df$rep, dfT$rep),
                     row.names=c(paste(df$stage, df$rep, sep = "-"), paste(dfT$stage, dfT$rep, sep = "-")))

res <- SummarizedExperiment(assays=list(counts=merged), rowRanges=tx, colData=colData)

fpkmCollapsed <- scoreSummarizedExperiment(res, "fpkm")

dds <- collapseReplicates(res, groupby = rep.int(1, ncol(assay(res))))
saveRDS(dds, "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/tables/experiments_merged.rds")
isNeuroPep <- (rownames(dds) %in% genes) & (rowMax(fpkmCollapsed) > 1)
#assay(dds) <- log10(assay(dds))
p <- pheatmap(assay(dds)[isNeuroPep,], cluster_rows=TRUE, show_rownames=FALSE,
              cluster_cols=FALSE);
p
ggsave("merged_heatmap.pdf", p, height = 30, limitsize = FALSE)
write.table(isNeuroPep, "GenesUsed/neuroPeptideGenes_merged.csv")

############################ tmhmm ############################

genesHmm <- getPeptides("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/tmhmm/output_hmm.txt")
finalGenes <- genes[genes %in% genesHmm]
write.table(finalGenes, file = "GenesUsed/PeptideGenes_AfterTmhmm.csv")

############################# POST SCRIPTING ###########################
all_genes <- read.table("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/allGenes.csv")$x
day_night <- read.table("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/neuroPeptideGenes_day_night.csv")$x
devel <- read.table("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/neuroPeptideGenes_development.csv")$x
brain_devel <- read.table("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/neuroPeptideGenes_brain_devel.csv")$x
merged <- read.table("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/neuroPeptideGenes_merged.csv")$x
tmhmm <- read.table("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/PeptideGenes_AfterTmhmm.csv")$x
# tests
c(length(all_genes), length(day_night),length(devel),length(brain_devel), length(merged)) == length(all_genes)
all(tmhmm %in% all_genes)
sum(day_night);sum(devel);sum(brain_devel);sum(merged)




day_night[!(all_genes[day_night] %in% tmhmm)] <- FALSE
devel[!(all_genes[devel] %in% tmhmm)] <- FALSE
brain_devel[!(all_genes[brain_devel] %in% tmhmm)] <- FALSE
merged[!(all_genes[merged] %in% tmhmm)] <- FALSE

sum(day_night);sum(devel);sum(brain_devel);sum(merged)

dir.create("/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/", showWarnings = FALSE)
write.table(all_genes, "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/allGenes.csv")
write.table(day_night, "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/neuroPeptideGenes_day_night.csv")
write.table(devel, "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/neuroPeptideGenes_development.csv")
write.table(brain_devel, "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/neuroPeptideGenes_brain_devel.csv")
write.table(merged, "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/neuroPeptideGenes_merged.csv")
write.table(tmhmm, "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/PeptideGenes_AfterTmhmm.csv")
write.table(all_genes[merged], "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/GenesUsed/new/final_genes.csv")
