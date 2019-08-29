# Create coverage plots from 3' UTR miR430 targets
setwd("/export/valenfs/projects/uORFome/RCode1/") 
source("./pipelineSetup.R")
source("./TempScripts/tcp_pipeline.R")
source("./TempScripts/MZdicer.R")
source("./SummarizedExperimentHelpers.R")

# Get annotation
txdb <- loadTxdb("/export/valenfs/projects/uORFome/Annotations/Zebrafish/Danio_rerio.Zv9.79.gtf.db")
seqlevelsStyle(txdb) <- "NCBI"
tx <- loadRegion(txdb, "mrna")
df <- getmzDicerDf()

# miR430 target genes
genes <- getMir430()

# Get miR430 target transcripts
len <- GenomicFeatures::transcriptLengths(txdb)
match <- len[len$gene_id %in% genes,]
match <- match[match$tx_name %in% names(tx),]
txNames <- match$tx_name
length(txNames)

####################################### FPKM PlOTS ######################################
res <- makeSummarizedExperimentFromBam(df, txdb, saveName = "/export/valenfs/projects/Håkon/mir430/countList.rds",
                                       longestPerGene = TRUE, 
                                       geneOrTxNames = "gene")
dt <- scoreSummarizedExperiment(res, score = "fpkm")
dt <- as.data.table(dt)
dif <- SEdif(dt)
d <- SESplit(dif, names(tx) %in% txNames)

lim <- 100
p <- ggplot() +
  geom_point(data = d[ismiRNA == F,], aes(x = RNA, y = RFP, alpha = 0.1, color = "gray")) + 
  geom_point(data = d[ismiRNA == T,], aes(x = RNA, y = RFP, alpha = 0.1, color = "red")) + 
  xlim(-lim, lim) +
  ylim(-lim, lim) + 
  xlab(expression("WT vs MZ dicer "*Delta*" mRNA")) + 
  ylab(expression("WT vs MZ dicer "*Delta*" RFP")) + 
  scale_color_manual(values= c("gray", "red")) +
  facet_grid( ~ stage, scales = "free")
p  
ggsave(p, filename = "/export/valenfs/projects/Håkon/mir430/2dFigure.png", width = 10, height = 4)

############################### COVERAGE PLOTS #####################################################
longEnough <- filterTranscripts(txdb, 100, 100, 100)
validMir <- txNames[txNames %in% longEnough]

# Create coverage plots of miR430 targets
leaders = fiveUTRsByTranscript(txdb,use.names = T)[validMir]
cds <- cdsBy(txdb,"tx", use.names = TRUE)[validMir]
trailers = threeUTRsByTranscript(txdb, use.names = TRUE)[validMir]
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[validMir]
transcriptWindow(leaders, cds, trailers, df = df, outdir = p(mainFolder, "/tcp_plots/mir430/targets_"), 
                 allTogether = TRUE)

# Create coverage of non-miR430 targets
validNames <- longEnough
validNames <- validNames[!(validNames %in% txNames)]
leaders = fiveUTRsByTranscript(txdb,use.names = T)[validNames]
cds <- cdsBy(txdb,"tx", use.names = TRUE)[validNames]
trailers = threeUTRsByTranscript(txdb, use.names = TRUE)[validNames]
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[validNames]
transcriptWindow(leaders, cds, trailers)

# Interest Region TIS Cover plots 
cds <- cdsBy(txdb,"tx", use.names = TRUE)[txNames] # targets
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[txNames]
regionWindow(cds, tx)

cds <- cdsBy(txdb,"tx", use.names = TRUE)[validNames] # non-targets
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[validNames]
regionWindow(cds, tx, outdir = p(mainFolder, "/tcp_plots/TIS_region/normal_"))

# Interest Region STOP site Cover plots 
trailers = threeUTRsByTranscript(txdb, use.names = TRUE)[txNames] # targets
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[txNames]
regionWindow(trailers, tx, outdir = p(mainFolder, "/tcp_plots/STOP_region/target_"))

trailers = threeUTRsByTranscript(txdb, use.names = TRUE)[validNames] # non-targets
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[validNames]
regionWindow(trailers, tx, outdir = p(mainFolder, "/tcp_plots/STOP_region/normal_"))


