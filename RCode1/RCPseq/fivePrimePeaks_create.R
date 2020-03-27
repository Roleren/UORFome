library(ORFikPipeline)
Palette1 <- c('skyblue4', 'orange')


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Use TISU definition to find SE difference in 4Ei vs WT
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create new fivePrimePeak plots

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
dfr <- read.experimentl("Val19", expInVarName = TRUE)
loadRegions(dfr, c("leaders", "cds"), names.keep = filterTranscripts(dfr, 1, 1, 0))
leadersCage <- loadRegion("/export/valenfs/projects/uORFome/Annotations/Zebrafish/shield_transcripts.gff3.db",
                          "leaders")
leadersCage <- leadersCage[names(leadersCage) %in% names(leaders)]
leadersTemp <- leaders;leadersTemp[names(leadersCage)] <- leadersCage;leadersCage <- leadersTemp
length(leadersCage) == length(leaders); identical(names(leadersCage),names(leaders))

# TISU transcripts
genes <- unique(readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/TISU-SE.rds")$gene_id)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load libraries
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# SSU
SSUWT64     <- fimport("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_SSU_reps_1_2_peaks_removed_translating_filter.bam", leaders)
SSUWTSHIELD <- readBam("/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_SSU_reps_1_2_3_peaks_removed_translating_filter.bam", leaders)
# 4Ei 10 uM (64cell)
df4Ei    <- read.experimentl("val_4Ei64", expInVarName = TRUE)[2,]; outputLibs(df4Ei, chrStyle = leaders)
# 4Ei 10 uM (shield)
df4EiShi <- read.experimentl("val_4EiShi.csv", expInVarName = TRUE)[2,]; outputLibs(df4EiShi, chrStyle = leaders)
# 4Ei 0.1 uM (64 cell)
df4Ei0.1 <- read.experimentl("val_10", expInVarName = TRUE)[3, ]; outputLibs(df4Ei0.1, chrStyle = leaders)

# RNA controls
dfrs <- dfr[c(7, 6, 2, 3, 11, 1, 4),]
outputLibs(dfrs, leaders)
merged64  <- c(Val19_RNA_64cell_r3, Val19_RNA_64cell_r1)
mergedShi <- c(Val19_RNA_Shield_r1, Val19_RNA_Shield_r2, Val19_RNA_Shield_r3)
rna64fpkm        <- fpkm(cds, merged64, librarySize = "overlapping")
rnaShieldfpkm    <- fpkm(cds, mergedShi, librarySize = "overlapping")
rna4Ei64fpkm     <- fpkm(cds, Val19_RNA_4Ei_64cell_rNA, librarySize = "overlapping")
rna4EiShieldfpkm <- fpkm(cds, Val19_RNA_4Ei_Shield_rNA, librarySize = "overlapping")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# FPKM tables
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# SSU
dt.ssu <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leaders), dfr), txNames = names(leaders),
                     leaders64_fpkm = fpkm(leaders, SSUWT64, librarySize = "overlapping"),
                     leadersShield_fpkm = fpkm(leaders, SSUWTSHIELD, librarySize = "overlapping"),
                     leaders4Ei64_fpkm = fpkm(leaders, val_4Ei64_SSU, librarySize = "overlapping"),
                     leaders4EiShi_fpkm = fpkm(leaders, val_4EiShi_SSU, librarySize = "overlapping"))
# SSU TSS PEAKS
tss <- startSites(leadersCage, asGR = T)
SSUWT64TSS        <- SSUWT64[countOverlaps(SSUWT64, tss) > 0]
SSUWTSHIELDTSS    <- SSUWTSHIELD[countOverlaps(SSUWTSHIELD, tss) > 0]
val_4Ei64_SSUTSS  <- val_4Ei64_SSU[countOverlaps(val_4Ei64_SSU, tss) > 0]
val_4EiShi_SSUTSS <- val_4EiShi_SSU[countOverlaps(val_4EiShi_SSU, tss) > 0]
# 4Ei 0.1uM only for TSS
val_10_SSUTSS <- val_10_SSU[countOverlaps(val_10_SSU, tss) > 0]

dt.ssu.tss <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leadersCage), dfr), txNames = names(leadersCage),
                         leaders64_fpkm = fpkm(leadersCage, SSUWT64TSS, librarySize = sum(countOverlaps(SSUWT64TSS, leadersCage) > 0)),
                         leadersShield_fpkm = fpkm(leadersCage, SSUWTSHIELDTSS, librarySize = sum(countOverlaps(SSUWTSHIELDTSS, leadersCage) > 0)),
                         leaders4Ei64_fpkm = fpkm(leadersCage, val_4Ei64_SSUTSS, librarySize = sum(countOverlaps(val_4Ei64_SSUTSS, leadersCage) > 0)),
                         leaders4EiShi_fpkm = fpkm(leadersCage, val_4EiShi_SSUTSS, librarySize = sum(countOverlaps(val_4EiShi_SSUTSS, leadersCage) > 0)))

leadersCage100 <- leadersCage[widthPerGroup(leadersCage, FALSE) >= 100]
tss100 <- startSites(leadersCage100, asGR = T)
nottss100 <- startRegion(leadersCage100, upstream =  -1, downstream = 99)
dt.ssu.tssCounts <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leadersCage100), dfr), txNames = names(leadersCage100),
                               leaders64_counts = countOverlaps(tss100, convertToOneBasedRanges(SSUWT64TSS)),
                               leaders4Ei64_counts = countOverlaps(tss100, convertToOneBasedRanges(val_4Ei64_SSUTSS)),
                               leaders4Ei_0.1_64_counts = countOverlaps(tss100, convertToOneBasedRanges(val_10_SSUTSS)),
                               leadersShi_counts = countOverlaps(tss100, convertToOneBasedRanges(SSUWTSHIELDTSS)),
                               leaders4EiShi_counts = countOverlaps(tss100, convertToOneBasedRanges(val_4EiShi_SSUTSS)))

dt.ssu.nottssCounts <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(leadersCage100), dfr), txNames = names(leadersCage100),
                                 leaders64_counts = countOverlaps(nottss100, convertToOneBasedRanges(SSUWT64)),
                                 leaders4Ei64_counts = countOverlaps(nottss100, convertToOneBasedRanges(val_4Ei64_SSU)),
                                 leaders4Ei_0.1_64_counts = countOverlaps(nottss100, convertToOneBasedRanges(val_10_SSU)),
                                 leadersShi_counts = countOverlaps(nottss100, convertToOneBasedRanges(SSUWTSHIELD)),
                                 leaders4EiShi_counts = countOverlaps(nottss100, convertToOneBasedRanges(val_4EiShi_SSU)))


# RNA
dt.rna <- data.table(X.gene_id = ORFik:::txNamesToGeneNames(names(cds), dfrs), txNames = names(cds),
                     rna64fpkm, rnaShieldfpkm, rna4Ei64fpkm, rna4EiShieldfpkm)
dt.rna <- dt.rna[txNames %in% dt.ssu$txNames,]
#dt.rna[,c(3:6)] <- dt.rna[,c(3:6)] + 1
identical(dt.rna$txNames, dt.ssu$txNames)
################ FPKM filters
# Filters for RNA-seq
rnaFilt <- 10
validInAll <- rowSums(dt.rna[,3:6] > rnaFilt) > 4
validIn64WT <- rowSums(dt.rna[,c(3) ] > rnaFilt) == 1
validIn64 <- rowSums(dt.rna[,c(3, 5) ] > rnaFilt) == 2
validInShiWT <- rowSums(dt.rna[,c(4) ] > rnaFilt) == 1
validInShi <- rowSums(dt.rna[,c(4, 6) ] > rnaFilt) == 2
# Filters for SSU libs
ssuFilt <- 1
validIn64WTSSU <- rowSums(dt.ssu[,c(3) ] > ssuFilt) == 1
validInShiWTSSU <- rowSums(dt.ssu[,c(4) ] > ssuFilt) == 1
# Filter for SSU TSS
validIn64WTSSUTSS <- rowSums(dt.ssu.tss[,c(3) ] > ssuFilt) == 1
validInShiWTSSUTSS <- rowSums(dt.ssu.tss[,c(4) ] > ssuFilt) == 1
# Filter for SSU TSS (tss vs downstream to 100)
validInTSSCounts <- rowSums(dt.ssu.tssCounts[,c(3) ] > 4) == 1

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# MAKE SE
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

cleanAndSave <- function(dt, saveName) {
  res <- melt(dt, id.vars = c("X.gene_id", "txNames", "type"), value.name = "SE")
  res$stage <- "64cell"
  res[grep("shield", variable), stage := "shield"]
  res$condition <- "WT"
  res[grep("4Ei", variable), condition := "4Ei"]
  saveRDS(res, saveName)
}
##################### MAKE SE without TSS filter
dt <- data.table(dt.ssu[,1:2], (dt.ssu[, 3:6] + 0.1) / dt.rna[,3:6] )
dt$type <- "Background"
dt[X.gene_id %in% genes, type := "TISU"]
colnames(dt) <- c("X.gene_id", "txNames", "64_SE" , "shield_SE", "64_4Ei_SE", "shield_4Ei_SE", "type")
# 64
cleanAndSave(dt[validIn64 & validIn64WTSSU,],
             "/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_64cell.rds")
# shield
cleanAndSave(dt[validInShi & validInShiWTSSU,],
             "/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_shield.rds")

##################### MAKE SE with TSS filter
### TSS (filtering out SSU reads not overlapping TSS)
dt <- data.table(dt.ssu.tss[,1:2], (dt.ssu.tss[, 3:6] + 0.1) / dt.rna[,3:6])
dt$type <- "Background"
dt[X.gene_id %in% genes, type := "TISU"]
colnames(dt) <- c("X.gene_id", "txNames", "64_SE" , "shield_SE", "64_4Ei_SE", "shield_4Ei_SE", "type")
# 64
cleanAndSave(dt[validIn64 & validIn64WTSSUTSS,],
             "/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_tss_64cell.rds")
# Shield
cleanAndSave(dt[validInShi & validInShiWTSSUTSS,],
             "/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_tss_shield.rds")

####################### MAKE SSU RATIO with TSS filter (vs downstream to 100)
# With pseudo counts
dt <- data.table(dt.ssu.tssCounts[,1:2], (dt.ssu.tssCounts[, 3:7] + 0.1) / (dt.ssu.nottssCounts[,3:7] + 0.1))
dt$type <- "Background"
dt[X.gene_id %in% genes, type := "TISU"]
colnames(dt) <- c("X.gene_id", "txNames", "64_WT_counts" , "64_4Ei_10_counts", "64_4Ei_0.1_counts",
                  "Shi_WT_counts", "Shi_4Ei_10_counts", "type")
res <- melt(dt[validInTSSCounts,], id.vars = c("X.gene_id", "txNames", "type"), value.name = "counts")

# pseud count Non
dt.ssu.tssCounts.m <- copy(dt.ssu.tssCounts)
dt.ssu.tssCounts.m$leaders4Ei_10_counts <- dt.ssu.tssCounts.m$leaders4Ei64_counts + dt.ssu.tssCounts.m$leaders4EiShi_counts
dt.ssu.tssCounts.m$leaders4Ei64_counts <- NULL
dt.ssu.tssCounts.m$leaders4EiShi_counts <- NULL
dt.ssu.tssCounts.m$leadersShi_counts <- NULL

dt.ssu.nottssCounts.m <- copy(dt.ssu.nottssCounts)
dt.ssu.nottssCounts.m$leaders4Ei_10_counts <- dt.ssu.nottssCounts.m$leaders4Ei64_counts + dt.ssu.nottssCounts.m$leaders4EiShi_counts
dt.ssu.nottssCounts.m$leaders4Ei64_counts <- NULL
dt.ssu.nottssCounts.m$leaders4EiShi_counts <- NULL
dt.ssu.nottssCounts.m$leadersShi_counts <- NULL

dt <- data.table(dt.ssu.tssCounts.m[,1:2], (dt.ssu.tssCounts.m[, 3:5]) / (dt.ssu.nottssCounts.m[,3:5]))
dt$type <- "Background"
dt[X.gene_id %in% genes, type := "TISU"]
colnames(dt) <- c("X.gene_id", "txNames", "64_WT_counts" , "64_4Ei_0.1_counts",
                  "4Ei_10_counts", "type")
res <- melt(dt[validInTSSCounts,], id.vars = c("X.gene_id", "txNames", "type"), value.name = "counts")
saveRDS(res, "/export/valenfs/projects/Hakon/RCP_SEQ/matrix_new2020feb_tss_countsratio.rds")