### DATA ###
# Load Experiments
dfA <- getArcher19()
dfA <- dfA[dfA$libtype == "SSU",]
dfA@expInVarName <- FALSE
dfa <- getArcher16()
dfa <- dfa[dfa$libtype == "SSU",]
dfl <- list(dfA, dfa)
names(dfl) <- c("19", "16")

# Annotation
outFolder <- "/export/valenfs/projects/Håkon/RCP_SEQ/plots/"
txdb <- loadTxdb(dfA@txdb)
loadRegions(txdb, extension = "")
txNames <- filterTranscripts(txdb, 60, 100, 60, F)
leaders <- leaders[txNames]; cds <- cds[txNames];trailers <- trailers[txNames];mrna <- mrna[txNames]
# Output data
outputLibs(dfA)
outputLibs(dfa)

### CAGE ###
# CAGE DATA test
CAGE <- fread.bed("/export/valenfs/data/processed_data/CAGE/McManus_2017_yeast/aligned/merged_CAGE_5prime.bed")
CAGE_old <- fread.bed("/export/valenfs/data/processed_data/CAGE/wery_2015_S_cerevisiae/final_results/aligned_R64_1_1/merged_CAGE_5prime.bed")
leadersCage <- reassignTSSbyCage(leaders, CAGE, removeUnused = F, cageMcol = F)
leadersCage <- leadersCage[widthPerGroup(leadersCage) >= 60]
leadersCage_old <- reassignTSSbyCage(leaders, CAGE_old, removeUnused = F, cageMcol = F)
leadersCage_old <- leadersCage[widthPerGroup(leadersCage) >= 60]

### PLOTS ###
#1 Metacoverage
summary(readWidths(SSU_CS_r1))
transcriptWindow(leaders, cds, trailers, df = dfl, outdir = outFolder, allTogether = TRUE)
#2 TSS (5' and 3')
## Archer 16 and 19 data
#non CAGE
interval <- 18:43
heatMapL(leaders, extendLeaders(mrna, 76), dfl, upstream = 75, downstream = 59, 
            outdir = outFolder, acLen = interval,
            scores = "zscore", shifting = "5prime", location = "TSS", skip.last = T)
heatMapL(leaders, extendLeaders(mrna, 76), dfl, upstream = 30, downstream = 74, 
         outdir = outFolder, acLen = 18:43,
         scores = "zscore", shifting = "3prime", location = "TSS", skip.last = T)
# CAGE
heatMapL(leadersCage, extendLeaders(leadersCage, 76), dfl, upstream = 75, downstream = 59, 
         outdir = p(outFolder), acLen = interval,
         scores = "zscore", shifting = "5prime", location = "TSS_CAGE_", skip.last = T)
heatMapL(leadersCage_old, extendLeaders(leadersCage_old, 76), dfl, upstream = 75, downstream = 59, 
         outdir = p(outFolder), acLen = interval,
         scores = "zscore", shifting = "5prime", location = "TSS_CAGE_old_", skip.last = T)
# TEST
# new CAGE
region <- ORFik:::startRegion(leadersCage, leadersCage, upstream = -2, downstream = 14)
counts <- countOverlaps(region, SSU_CS_r1);summary(counts); sum(counts > 10)
yeastLeadersCage <- leadersCage[(counts > 10)]
counts <- fpkm(yeastLeadersCage, SSU_CS_r1);summary(counts)
yeastLeadersCage <- yeastLeadersCage[counts > 93]
tcpHeatMap_single(yeastLeadersCage, extendLeaders(yeastLeadersCage, extension = 50), SSU_CS_r2, upstream = 30, downstream = 45, outdir = paste0(heatMapsFolder, "test_yeast_SSU_5prime_TSS_trans.pdf"), scores = "transcriptNormalized", acLen = 16:45)
# OLD CAGE
region <- ORFik:::startRegion(leadersCage_old, leadersCage_old, upstream = -2, downstream = 14)
counts <- countOverlaps(region, SSU_CS_r1);summary(counts); sum(counts > 200)
yeastLeadersCage_old <- leadersCage_old[!(counts > 200)]
counts <- fpkm(yeastLeadersCage_old, SSU_CS_r1);summary(counts)
yeastLeadersCage_old <- yeastLeadersCage_old[counts > 5]
tcpHeatMap_single(leadersCage_old, extendLeaders(leadersCage_old, extension = 50), SSU_CS_r2, upstream = 30, downstream = 45, outdir = paste0(heatMapsFolder, "test_yeast_SSU_5prime_TSS_trans.pdf"), scores = "transcriptNormalized", acLen = 16:70)
#3 TIS (5' and 3')
interval <- 18:45
heatMapL(cds, mrna, dfA, upstream = 59, downstream = 74, 
            outdir = outFolder, acLen = interval,
            scores = "zscore", shifting = "5prime")

# Make countOverlaps with score
heatMapL(yeastLeadersCage, extendLeaders(yeastLeadersCage, 76), dfl, upstream = 75, downstream = 59, 
         outdir = p(outFolder), acLen = interval,
         scores = "zscore", shifting = "5prime", location = "TEST_TSS_CAGE_", skip.last = T)
heatMapL(yeastLeadersCage_old, extendLeaders(yeastLeadersCage_old, 76), dfl, upstream = 75, downstream = 59, 
         outdir = p(outFolder), acLen = interval,
         scores = "zscore", shifting = "5prime", location = "TEST_TSS_CAGE_old_", skip.last = T)
