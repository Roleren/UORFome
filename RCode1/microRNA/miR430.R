# Create coverage plots from 3' UTR miR430 targets
source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")
source("./microRNA/MZdicer.R")
hm <- "/export/valenfs/projects/Hakon/mir430/heatmaps/"

# Pick which dataset
dfB <- getBazzini12() # bazzini
dfR <- getRCPMZ() # RCP

# Get annotation
txdbB <- loadTxdb(dfB, "NCBI")
txdbR <- loadTxdb(dfR, "NCBI")
loadRegions(txdbB, extension = "B")
loadRegions(txdbR, extension = "R")

# miR430 target transcripts, filter on score -0.20
txNamesB <- getMir430Tx(txdb = txdbB, tx = mrnaB)
txNamesR <- getMir430Tx(txdb = txdbR, tx = mrnaR)

####################################### FPKM PlOTS ######################################
# Bazzini replicate (tx)
dt <- makeSummarizedExperimentFromBam(dfB, txdbB, saveName = "/export/valenfs/projects/Hakon/mir430/count_lists/countList.rds",
                                       longestPerGene = TRUE, geneOrTxNames = "gene", type = "fpkm")
dif <- SEdif(dt, dfB)
#dB <- SESplit(dif, names(mrnaB)[names(mrnaB) %in% txNamesB] %in% names(leadersB),  dtS = dt[, c(1,2, 7, 8)], df = dfB)
dB <- SESplit(dif, names(mrnaB) %in% txNamesB,  dtS = dt[, c(1,2, 7, 8)], df = dfB)
difPlot(dB, lim = 6)

# RCP version (mrna)
dtR <- makeSummarizedExperimentFromBam(dfR, txdbR, saveName = "/export/valenfs/projects/Hakon/mir430/count_lists/countList_RCP.rds",
                                       longestPerGene = TRUE, geneOrTxNames = "gene", type = "fpkm")
dif <- SEdif(dtR, dfR)
#dRT <- SESplit(dif, names(mrnaR)[names(mrnaR) %in% txNamesR] %in% names(leadersR), dtS = dtR[, c(1,2,3, 7, 8, 9)], df = dfR, score = 2)
dRT <- SESplit(dif, names(mrnaR) %in% txNamesR, dtS = dtR[, c(1,2,3, 7, 8, 9)], df = dfR, score = 7)
difPlot(dRT, filename = c("/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_LSU.png", 
                        "/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_SSU.png"),
        type = "RCP", lim = 6)


# Use mrna RNA-seq here
# RCP version (leader)
dt <- makeSummarizedExperimentFromBam(dfR, txdbR, saveName = "/export/valenfs/projects/Hakon/mir430/count_lists/countList_RCP_leaders.rds",
                                        longestPerGene = TRUE, geneOrTxNames = "gene", region = "leaders", type = "fpkm")
rnaInds <- (1:ncol(dt))[1:ncol(dt) %% 3 == 0]
dt[, rnaInds] <- dtR[names(mrnaR) %in% names(leadersR), rnaInds, with = F]

dif <- SEdif(dt, dfR)
dR <- SESplit(dif, names(leadersR) %in% txNamesR, dtS = dt[, c(1,2,3, 7, 8, 9)],
              df = dfR, score = 2)
difPlot(dR, filename = c("/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_LSU_leaders.png", 
                         "/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_SSU_leaders.png"),
        type = "RCP", lim  = 30)

# RCP version (CDS)
dt <- makeSummarizedExperimentFromBam(dfR, txdbR, saveName = "/export/valenfs/projects/Hakon/mir430/count_lists/countList_RCP_cds.rds",
                                        longestPerGene = TRUE, geneOrTxNames = "gene", region = "cds", type = "fpkm")
dt[, rnaInds] <- dtR[names(mrnaR) %in% names(cdsR), rnaInds, with = F]
dif <- SEdif(dt, dfR)
dR <- SESplit(dif, names(cdsR) %in% txNamesR, dtS = dtR[names(mrnaR) %in% names(cdsR), c(1,2,3, 7, 8, 9)],
              df = dfR, score = 2)
difPlot(dR, filename = c("/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_LSU_cds.png", 
                         "/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_SSU_cds.png"),
        type = "RCP")

# RCP version (trailer)
dt <- makeSummarizedExperimentFromBam(dfR, txdbR, saveName = "/export/valenfs/projects/Hakon/mir430/count_lists/countList_RCP_trailers.rds",
                                        longestPerGene = TRUE, 
                                        geneOrTxNames = "gene", region = "trailers", type = "fpkm")
dt[, rnaInds] <- dtR[names(mrnaR) %in% names(trailersR), rnaInds, with = F]
dif <- SEdif(dt, dfR)
dR <- SESplit(dif, names(trailersR) %in% txNamesR, dtS = dtR[names(mrnaR) %in% names(trailersR), c(1,2,3, 7, 8, 9)],
              df = dfR, score = 2)
difPlot(dR, filename = c("/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_LSU_trailers.png", 
                         "/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_RCP_SSU_trailers.png"),
        type = "RCP")

############################### COVERAGE PLOTS #####################################################
dfB <- dfB[dfB$libtype != "RNA",] # We don't need RNA samples loaded
dfR <- dfR[dfR$libtype != "RNA",]
outputLibs(dfB, mrnaB)
outputLibs(dfR, mrnaR)

#RNA normalizing
dfr <- getRCPMZ()
dfr <- dfr[dfr$libtype == "RNA",]
outputLibs(dfr, mrnaR)

# Split into target and non-target miR430 annotation
filterTranscriptsSplit(splitList = list(txNamesR, txNamesB), extensions = c("R", "B"), 
                       names = c("validMir", "validNames"))
splitRegions(splitList = list(validMirR, validNamesR, validMirB, validNamesB), 
             extensions = c("R", "B"), splitExt = c("T", "N"))

# Bazzini
# Create coverage plots of miR430 targets
transcriptWindow(leadersBT, cdsBT, trailersBT, df = dfB, outdir = p(hm,"cp_targets_baz_"), 
                 allTogether = TRUE, colors = c(rep("orange", 3), rep("skyblue4", 3)))

# Create coverage of non-miR430 targets
transcriptWindow(leadersBN, cdsBN, trailersBN, df = dfB, outdir = p(hm,"cp_nontargets_baz_"), 
                 allTogether = TRUE, colors = c(rep("orange", 3), rep("skyblue4", 3)))

# RCP-seq
# Create coverage plots of miR430 targets
transcriptWindow(leadersRT, cdsRT, trailersRT, df = dfR, outdir = p(hm,"cp_targets_RCP_"), 
                 allTogether = TRUE)

# Create coverage of non-miR430 targets
transcriptWindow(leadersRN, cdsRN, trailersRN, df = dfR, outdir = p(hm,"cp_nontargets_RCP_"), 
                 allTogether = TRUE)

# Region windows TSS
hits <- names(mrnaR[which(startSites(mrnaR) > 80)])

# Region windows TSS
regionWindow(leadersRT[names(leadersRT) %in% hits], extendLeaders(mrnaR[names(mrnaR) %in% hits]), dfR, outdir = hm, title = "TSS metacoverage targets")
regionWindow(leadersRN[names(leadersRN) %in% hits], extendLeaders(mrnaR[names(mrnaR) %in% hits]), dfR, outdir = hm, title = "TSS metacoverage Non-targets")

regionWindow(leadersRT[names(leadersRT) %in% hits], extendLeaders(mrnaR[names(mrnaR) %in% hits]), dfR, outdir = hm, title = "TSS metacoverage targets", dfr = dfr)
regionWindow(leadersRN[names(leadersRN) %in% hits], extendLeaders(mrnaR[names(mrnaR) %in% hits]), dfR, outdir = hm, title = "TSS metacoverage Non-targets", dfr = dfr)
# Interest Region STOP site Cover plots 
# TODO


############################### HEATMAP plots #####################################################
# TIS
interval <- 25:36

# Bazzini targets
heatMapL(cdsBT, mrnaBT, dfB, upstream = 75, downstream = 74, 
            outdir = p(hm, "hm_targets_"), 
            acLen = interval, scores = "sum", shifting = "5prime")
# Bazzini non-targets
heatMapL(cdsBN, mrnaBN, dfB, upstream = 75, downstream = 74, 
            outdir = p(hm, "hm_nontargets_"), 
            acLen = interval, scores = "sum", shifting = "5prime")

# RCP targets
heatMapL(cdsRT, mrnaRT, dfR, upstream = 75, downstream = 74, 
            outdir = p(hm, "hm_targets_"), 
            scores = "sum", shifting = "5prime", skip.last = T)
# RCP non-targets
heatMapL(cdsRN, mrnaRN, dfR, upstream = 75, downstream = 74, 
            outdir = p(hm, "hm_nontargets_"), 
            scores = "sum", shifting = "5prime", skip.last = T)

# TSS (Leader)
dfR$LSU <- NULL
# RCP targets
heatMapL(leadersRT, extendLeaders(mrnaRT, 76), dfR, upstream = 75, downstream = 74, 
            outdir = p(hm, "hm_targets_"), 
            scores = "sum", shifting = "5prime", location = "TSS", skip.last = T)
# RCP non-targets
heatMapL(leadersRN[!(names(leadersRN) %in% "ENSDART00000167860")], extendLeaders(mrnaRN, 76), 
            dfR, upstream = 75, downstream = 74, 
            outdir = p(hm, "hm_nontargets_"), 
            scores = "sum", shifting = "5prime", location = "TSS", skip.last = T)

# Region windows TIS
regionWindow(cdsRT, mrnaR, dfR, outdir = hm, title = "cds metacoverage targets")
regionWindow(cdsRN, mrnaR, dfR, outdir = hm, title = "cds metacoverage Non-targets")
