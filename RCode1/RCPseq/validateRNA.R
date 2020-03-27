library(ORFikPipeline)
dfr <- read.experimentl("Val19", expInVarName = FALSE)
loadRegions(dfr, c("leaders", "cds"))

# RNA seq
rna <- makeSummarizedExperimentFromBam(saveName = "/export/valenfs/data/processed_data/RNA-seq/Valen_2019_zebrafish_1/aligned/QC_STATS/countTable_cds.rds")
test <- makeSummarizedExperimentFromBam(dfr[1,], region = "cds", geneOrTxNames = "tx", longestPerGene = FALSE)

# From test
ds <- DESeq2::DESeqDataSet(test, design = ~1)
head(DESeq2::fpkm(ds, robust = TRUE))
head(DESeq2::fpkm(ds, robust = FALSE))
fin <- DESeq2::fpkm(ds, robust = FALSE)

# From rna
dsrna <- DESeq2::DESeqDataSet(rna, design = ~SAMPLE)
head(DESeq2::fpkm(dsrna, robust = TRUE))
head(DESeq2::fpkm(dsrna, robust = FALSE))
fin2 <- DESeq2::fpkm(dsrna, robust = FALSE)[,1]
# Conclusion: RNA is equal when robust == FALSE



############### See if you can remake deseq2 fpkm
adams <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/adam_matrix_combined_nonfiltered_new_RNAseq.rds")
adams <- adams[stage == "64cell",]

fin2 <- DESeq2::fpkm(dsrna, robust = FALSE)[,1]
rnaSample <- assay(rna)[,1]
res <- ORFik:::fpkm_calc(rnaSample, widthPerGroup(rowRanges(rna), FALSE), librarySize = sum(rnaSample))
resORI <- ORFik:::fpkm(cds, RNA, librarySize = "overlapping"); names(resORI) <- names(cds)
resORI[testNames]
#  Conclusion, equal when librarySize = sum(rnaSample)

############### In adam table new with test:

testNames <- head(adams$transcript_id)
fin[chmatch(head(adams$transcript_id),  rownames(fin), nomatch = 0),] # Pure deseq
hed <- head(adams$complete_CDS_totalRNA_FPKM_new);names(hed) <-  head(adams$transcript_id); hed # From matrix

# ORFik fpkm
# Ori
stdORFOri <- ORFik:::fpkm(cds, RNA); names(stdORFOri) <- names(cds)
stdORFOri[testNames]
# full
stdORF <- ORFik:::fpkm_calc(countOverlaps(cds, RNA), widthPerGroup(cds, TRUE), librarySize = length(RNA)); names(stdORF) <- names(cds)
stdORF[testNames]
# Only overlapping
stdORF1 <- ORFik:::fpkm_calc(countOverlaps(cds, RNA), widthPerGroup(cds, TRUE), librarySize = sum(countOverlaps(RNA, cds) > 0)); names(stdORF1) <- names(cds)
stdORF1[testNames]
