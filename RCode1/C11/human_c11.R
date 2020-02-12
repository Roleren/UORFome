############################## CREATE DATA ###################################################
source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")
df <- read.experimentl("Val20HumMer")
outputLibs(df, type = "bedo")

# heatmaps
heatMapRegion(df, c("TIS","TSS", "TTS"))
# With cage
heatMapRegion(df, c("TIS","TSS"), paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/TSS_cage/"),
              cage = "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz")


# Fraction 11-13 SSU
source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")
df <- read.experimentl("Val20HumMer")
outputLibs(df, type = "bedo")
heatMapRegion(df, c("TSS"), paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/TSS_cage/"), 
              cage = "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz")