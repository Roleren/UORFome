rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/") 
library(ORFik)
source("./pipelineSetup.R")
dbDisconnect(uorfDB)
rm(uorfDB)

source("./TempScripts/tcp_pipeline.R")

setwd(p(mainFolder, "/AdamVienna/"))
plotFolder <- "/export/valenfs/projects/uORFome/AdamVienna/plots/new_plots/"

gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath);
getFasta("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa")
seqlevelsStyle(txdb)  <- seqlevelsStyle(fa)
leaders <- fiveUTRsByTranscript(txdb, use.names = T)
cds <- cdsBy(txdb, use.names = TRUE)
tx <- exonsBy(txdb, use.names = TRUE)

txnames <- names(cds[widthPerGroup(cds) >= 50])
txnames <- txnames[txnames %in% names(leaders[widthPerGroup(leaders) >= 100])]
leaders <- leaders[txnames]
cds <- cds[txnames]
tx <- tx[txnames]
########################## heatmaps ###############################
df <- getTCPdfAll()
df$RFP <- NULL
df$LSU <- NULL
df$RNA <- NULL
gtfPath <- "/export/valenfs/projects/adam/TCP_seq/transcript_GFF3/sphere_transcripts.gff3"
txdb <- GenomicFeatures::makeTxDbFromGFF(gtfPath)
# SSU
tcpHeatMap(txdb, df = df, outdir = p(plotFolder, "heatmaps/heatMap_"), shifting = "5prime")
tcpHeatMap(txdb, df = df, outdir = p(plotFolder, "heatmaps/heatMap_"), shifting = "3prime")

# RFP
df <- getTCPdfAll()
df$SSU <- NULL
df$LSU <- NULL
df$RNA <- NULL
tcpHeatMap(txdb, df = df[3,], outdir = p(plotFolder, "heatmaps/heatMap_RFP_"), shifting = "5prime")
tcpHeatMap(txdb, df = df[3,], outdir = p(plotFolder, "heatmaps/heatMap_RFP_"), shifting = "3prime")

# LSU
df <- getTCPdfAll()
df$SSU <- NULL
df$RFP <- NULL
df$RNA <- NULL
tcpHeatMap_int(cds, tx, df = df[3,], outdir = p(plotFolder, "heatmaps/heatMap_LSU_"), shifting = "5prime",
               upstream = 100, downstream = 49)
tcpHeatMap_int(cds, tx, df = df[3,], outdir = p(plotFolder, "heatmaps/heatMap_LSU_"), shifting = "3prime",
               upstream = 100, downstream = 49)

# periodicity
periodicityChecker(txdb, df = df[3,], shifting = "5prime")