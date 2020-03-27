#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis of test data of C11 made by Yamilla, MARCH 2020

############################## CREATE DATA ###################################################
library(ORFikPipeline)
df <- read.experimentl("Val20HumMer")
outputLibs(df, type = "bedo")

# heatmaps
heatMapRegion(df, c("TIS","TSS", "TTS"))
# With cage
heatMapRegion(df, c("TIS","TSS"), paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/TSS_cage/"),
              cage = "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz")


# Merged Fraction 11-13 SSU
library(ORFikPipeline)
df <- read.experimentl("Val20HumMer")
outputLibs(df, type = "bedo")
heatMapRegion(df, c("TSS"), paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/TSS_cage/"),
              cage = "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz")


s#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Metacoverage, IR & TOP motif
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis march 2020

# 1 Metacoverage
library(ORFikPipeline)
df <- read.experimentl("Val20HumMer", expInVarName = TRUE)
outputLibs(df, type = "bedo")
dfr <- read.experimentl("Val20Hum")
tables <- countTable(dfr, region = "mrna", type = "fpkm")
tables <- tables[, "RNA_"]; summary(tables)
mrna <- loadRegion(df, "mrna")
rnaFilt <- (tables$RNA_ > 1)
filter <- rnaFilt  & (names(mrna) %in% filterTranscripts(df, 100, 100, 100))
mrna <- mrna[filter]

t <- removeBadTxByRegion(mrna, bamVarName(df)[4])
loadRegions(df, names.keep = names(t))
transcriptWindow(leaders, get("cds", mode = "S4"),
                 trailers, df = df, outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/"),
                 allTogether = TRUE, colors = c(rep("orange", 2),rep("skyblue4", 2)),
                 scores = c("sum", "zscore"))

# Use BAM files
df <- read.experimentl("Val20HumMer")
outputLibs(df, chrStyle = t)
transcriptWindow(leaders, get("cds", mode = "S4"),
                 trailers, df = df, outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/"),
                 allTogether = TRUE, colors = c(rep("orange", 2),rep("skyblue4", 2)),
                 scores = c("sum", "zscore"))
slackrUpload(filename = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/","val20humMer_cp_all_zscore_.png"), channels = "visualizations")


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# IR
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library(ORFikPipeline)
df <- read.experimentl("Val20HumMer", expInVarName = TRUE)
outputLibs(df, type = "bedo")
dfr <- read.experimentl("Val20Hum")
tables <- countTable(dfr, region = "mrna", type = "fpkm")
tables <- tables[, "RNA_"]; summary(tables)
mrnaOri <- loadRegion(df, "mrna")
rnaFilt <- (tables$RNA_ > 1)
mrna <- mrnaOri[rnaFilt]

loadRegions(df, parts = c("leaders", "cds"))

if (!all(names(leaders) %in% names(cds))) stop("Not all leaders present!")


dt <- data.table()
tablesLeaders <- countTable(df, "leaders", type = "fpkm")
tablesCDS <- countTable(df, "cds", type = "fpkm")
tablesCDS <- tablesCDS[(names(cds) %in% names(leaders)) & (names(cds) %in% names(mrna)),]
tablesLeaders <- tablesLeaders[names(leaders) %in% names(mrna),]
if (nrow(tablesLeaders) == nrow(tablesCDS)) stop("not matching rows!")

dt$txNames <- names(leaders)[names(leaders) %in% names(mrna)]
dt$LSU_CDS_FPKM <- tablesCDS$LSU_rp + tablesCDS$LSU_WT
dt$SSU_LEADERS_FPKM <- tablesLeaders$LSU_rp + tablesLeaders$LSU_WT
dt$RNA_MRNA_FPKM <- tables$RNA_[names(mrnaOri) %in% dt$txNames]
dtOri <- copy(dt)
dt <- dt[SSU_LEADERS_FPKM > 0,]
dt$IR <- dt$LSU_CDS_FPKM / dt$SSU_LEADERS_FPKM
dt$SE <- dt$SSU_LEADERS_FPKM / dt$RNA_MRNA_FPKM
dtOri2 <- copy(dt)

seqs <- startRegionString(cds[dt$txNames], tx = mrna, faFile = df, upstream = 6, downstream = 2)
hits <- nchar(seqs) == max(nchar(seqs)); summary(hits)
seqs <- seqs[hits]
kozakHeat <- kozakHeatmap(seqs = seqs, rate = dt$IR[hits],
                          start = 1, stop = max(nchar(seqs)), center = 7, type = "IR");kozakHeat
ggslackR()

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TOP motif (Effect on SE FPKM)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
cage <- "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz"

leadersCage <- reassignTSSbyCage(leaders[dt$txNames], cage)

seqs <- startRegionString(leadersCage, NULL, df, 0, 4)
rate <- dt$SE

comb <- TOP.Motif.ecdf(seqs, rate)
ggslackR(plot = comb)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# SSU processivity
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

