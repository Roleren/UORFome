#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# DATA
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load Packages 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
rm(list=ls())
library(ORFikPipeline)
source("/export/valenfs/projects/uORFome/RCode1/zebrafishDevelBrain/Peptide_analysis.R")
library(DESeq2)
library(SummarizedExperiment)
library(pheatmap)
library(BiocParallel)
register(MulticoreParam(20))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
setwd("/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep")
genesDir <- "./GenesUsed_fpkm0.1"

dfr <- read.experimentl("zf_brain_devel")
fa <- ORFik:::findFa(dfr)
txdb <- loadTxdb(dfr, seqlevelsStyle(fa));
tx <- loadRegion(txdb, "mrna")
names(tx) <- ORFik:::txNamesToGeneNames(names(tx), txdb)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# FILTERS (protein and fpkm)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Protein
#genes <- getPeptides()
write.table(genes, file = p(genesDir, "/PeptideGenes.csv"))
genes <- read.table(p(genesDir, "/PeptideGenes.csv"))$x
fpkm_filt <- 0.1

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Count tables (For the 3 experiments)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

############################ brain_development #######################################
libPath <- "/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/brain_development/aligned"

saveName <- "/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/tables/experiments_brain_devel.rds"
#final <- makeSummarizedExperimentFromBam(dfr, txdb, saveName)
final <- readRDS(saveName)
dds <- scoreSummarizedExperiment(final)
fpkmCollapsed <- scoreSummarizedExperiment(final, "fpkm")

# New heatmap with only neuro peptides
isNeuroPep <- (rownames(dds) %in% genes) & (rowMax(fpkmCollapsed) > fpkm_filt)
dff <- as.data.frame(colData(dds)[,c("SAMPLE","replicate")])
p <- pheatmap(assay(dds)[isNeuroPep,], cluster_rows=TRUE, show_rownames=FALSE,
              cluster_cols=FALSE, annotation_col=dff)
p
ggsave("brain_devel_heatmap.pdf", p, height = 30, limitsize = FALSE)
write.table(isNeuroPep, p(genesDir, "/neuroPeptideGenes_brain_devel.csv"))

############################ DAY / NIGHT ############################################
libPath2 <- "/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/brain_day_night/aligned"
#dfT <- create.experimentl(libPath2, exper = "zf_brain_dn")
dfT <- read.experimentl("zf_brain_dn")



saveName <- "/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/tables/experiments_brain_day_night.rds"
#finalTime <- makeSummarizedExperimentFromBam(dfT, txdb, saveName)
finalTime <- readRDS(saveName)

dds <- scoreSummarizedExperiment(finalTime)
fpkmCollapsed <- scoreSummarizedExperiment(finalTime, "fpkm")
isNeuroPep <- (rownames(dds) %in% genes) & (rowMax(fpkmCollapsed) > fpkm_filt)
dff <- as.data.frame(colData(dds)[,c("SAMPLE","replicate")])
p <- pheatmap(assay(dds)[isNeuroPep,], cluster_rows=TRUE, show_rownames=FALSE,
              cluster_cols=FALSE, annotation_col=dff)
p
ggsave("brain_day_night_heatmap.pdf", p, height = 30, limitsize = FALSE)
write.table(isNeuroPep, p(genesDir, "/neuroPeptideGenes_day_night.csv"))
de <- DESeqDataSet(finalTime, design = ~ SAMPLE)
dds <- DESeq(de)
res <- results(dds)

############################ Merged BRAIN ############################

merged <- cbind(assay(final), assay(finalTime))
colnames(merged) <- NULL

s <- c(as.character(colData(final)$SAMPLE), as.character(colData(finalTime)$SAMPLE))
s <- factor(s, levels = unique(s))
colData <- DataFrame(SAMPLE=s,
                     replicate=c(colData(final)$replicate, colData(finalTime)$replicate),
                     row.names=c(paste(colData(final)$SAMPLE, colData(final)$replicate, sep = "-"),
                                 paste(colData(finalTime)$SAMPLE, colData(finalTime)$replicate, sep = "-")))

res <- SummarizedExperiment(assays=list(counts=merged), rowRanges=rowRanges(final), colData=colData)

fpkmCollapsed <- scoreSummarizedExperiment(res, "fpkm")

dds <- collapseReplicates(res, groupby = rep.int(1, ncol(assay(res))))
saveRDS(dds, "/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/tables/experiments_merged.rds")
dds <- readRDS("/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/tables/experiments_merged.rds")
isNeuroPep <- (rownames(dds) %in% genes) & (rowMax(fpkmCollapsed) > fpkm_filt)
#assay(dds) <- log10(assay(dds))
p <- pheatmap(assay(dds)[isNeuroPep,], cluster_rows=TRUE, show_rownames=FALSE,
              cluster_cols=FALSE);
p
ggsave("merged_heatmap.pdf", p, height = 30, limitsize = FALSE)
write.table(isNeuroPep, p(genesDir, "/neuroPeptideGenes_merged.csv"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# tmhmm
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

genesHmm <- getPeptides("/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/tmhmm/output_hmm.txt")
finalGenes <- genes[genes %in% genesHmm]
write.table(finalGenes, file = p(genesDir, "/PeptideGenes_AfterTmhmm.csv"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# POST SCRIPTING (update tables to new gene folder)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
############################# POST SCRIPTING ###########################
all_genes <- read.table(p(genesDir, "/allGenes.csv"))$x
day_night <- read.table(p(genesDir, "/neuroPeptideGenes_day_night.csv"))$x
brain_devel <- read.table(p(genesDir, "/neuroPeptideGenes_brain_devel.csv"))$x
merged <- read.table(p(genesDir, "/neuroPeptideGenes_merged.csv"))$x
tmhmm <- read.table(p(genesDir, "/PeptideGenes_AfterTmhmm.csv"))$x
# tests
c(length(all_genes), length(day_night),length(brain_devel), length(merged)) == length(all_genes)
all(tmhmm %in% all_genes)
sum(day_night);sum(brain_devel);sum(merged)


day_night[!(all_genes[day_night] %in% tmhmm)] <- FALSE
brain_devel[!(all_genes[brain_devel] %in% tmhmm)] <- FALSE
merged[!(all_genes[merged] %in% tmhmm)] <- FALSE

sum(day_night);sum(brain_devel);sum(merged)


final_genes <- as.character(all_genes[merged])
dir.create(p(genesDir, "/with_tmm_filter/"), showWarnings = FALSE)
write.table(all_genes, p(genesDir, "/with_tmm_filter/allGenes.csv"))
write.table(day_night, p(genesDir, "/with_tmm_filter/neuroPeptideGenes_day_night.csv"))
write.table(brain_devel, p(genesDir, "/with_tmm_filter/neuroPeptideGenes_brain_devel.csv"))
write.table(merged, p(genesDir, "/with_tmm_filter/neuroPeptideGenes_merged.csv"))
write.table(tmhmm, p(genesDir, "/with_tmm_filter/PeptideGenes_AfterTmhmm.csv"))
write.table(all_genes[merged], p(genesDir, "/with_tmm_filter/final_genes.csv"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# ZFIN analysis
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# neuropeptides
#updated list of known neuropeptides (171) with zfin ID (columns: zfin ID, symbol, full name, species, class (r = receptor, p = peptide, s = peptide processing)
knownNeuroPep <- fread("20191127_Neuropeptides_Known_171.csv", header = F, col.names = c("ZFIN_ID", "gene", "description", "species", "type"))
ensembl <- useMart("ensembl", dataset="drerio_gene_ensembl")
ens_IDs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
      filters = "external_gene_name", 
      values = knownNeuroPep$gene, 
      mart = ensembl)
ens <- ens_IDs$ensembl_gene_id

sum(ens %in% all_genes)
sum(ens %in% final_genes)
final_genes_old <- read.table("/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/GenesUsed/new/final_genes.csv")$x
sum(ens %in% final_genes_old)
