#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Motif search: GAG-(X)3-ACA/G/C (not used in article)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# Test for RCP-seq paper
library(ORFikPipeline)
library(ggpubr)
library(Biostrings)
plotFolder <- "/export/valenfs/projects/Hakon/AdamVienna/plots/new_plots/"

# Anotation
gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath)
getFasta("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa")
seqlevelsStyle(txdb)  <- seqlevelsStyle(fa)
loadRegions(txdb, "leaders")
getSequencesFromFasta(leaders)

# Pairwise alignment
# GAG-(X)3-ACA/G/C
motif <- "GAGNNNACN"
res <- pairwiseAlignment(seqs, motif, type = "local-global")
leadersER <- leaders[!(substr(res@pattern, 9, 9) == "T")]
res <- res[!(substr(res@pattern, 9, 9) == "T")]
leadersER <- leadersER[res@score >= max(res@score)]
res <- res[res@score >= max(res@score)]

# Check with ER tx list
ready <- readRDS("expression_both_withGorilla_filteredRFPSSU.rds")
length(leadersER)
dim(ready)
dim(ready[fraction == "GO_ER",])

sum(names(leadersER) %in% ready$transcript_id)
sum(names(leadersER) %in% ready[fraction == "GO_ER",]$transcript_id)

