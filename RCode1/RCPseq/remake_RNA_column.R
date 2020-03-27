# Using ORFik fpkm calc for fpkm
library(ORFikPipeline)
dt <-  readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")

############## load rna exp
dfr <- read.experimentl("Val19", expInVarName = TRUE);
loadRegions(dfr, c("cds"))
dfrs <- dfr[c(7, 10, 2),]; dfrs # Subset experiment
outputLibs(dfrs, cds)

fpkm64 <- fpkm(cds, Val19_RNA_64cell, librarySize = "overlapping")
fpkmSphere <- fpkm(cds, Val19_RNA_Sphere, librarySize = "overlapping")
fpkmShield <- fpkm(cds, Val19_RNA_Shield, librarySize = "overlapping")
fpkm <- melt(data.table("64cell" = fpkm64, sphere = fpkmSphere, shield = fpkmShield, transcript_id = names(cds)),
             id.vars = "transcript_id",
             variable.name = "stage", value.name = "complete_CDS_totalRNA_FPKM_new2")
fpkmMerge <- data.table::merge.data.table(dt, fpkm, by = c("transcript_id", "stage"), all.x = TRUE, all.y = FALSE)
fpkmMerge[is.na(complete_CDS_totalRNA_FPKM_new2), complete_CDS_totalRNA_FPKM_new2 := 0]

saveRDS(fpkmMerge, "/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq2020feb.rds")


# Statistics
summary(fpkmMerge$complete_CDS_totalRNA_FPKM_new) # DESeq
summary(fpkmMerge$complete_CDS_totalRNA_FPKM_new2) # corrected
cor.test(fpkmMerge$complete_CDS_totalRNA_FPKM_new, fpkmMerge$complete_CDS_totalRNA_FPKM_new2)
