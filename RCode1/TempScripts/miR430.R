# Create coverage plots from 3' UTR miR430 targets
setwd("/export/valenfs/projects/uORFome/RCode1/") 
source("./pipelineSetup.R")
source("./TempScripts/tcp_pipeline.R")

# miR430 target genes
path <- "/export/valenfs/projects/adam/TCP_seq/data_files/targetscanFish_180829_tidy.csv"
targets <- read.csv(path)
targets <- targets[!(is.na(targets$miR430_score) | is.nan(targets$miR430_score)),]
targets <- targets[targets$miR430_score <= -0.20,]
genes <- targets$gene_id

# Get annotation
# gtfPath <- "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81.gtf"
gtfPath <- p(mainFolder, "/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath); seqlevelsStyle(txdb) <- "NCBI"


# txdb <- makeTxDbFromGFF(file = gtfPath)
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)

# Get miR430 target transcripts
len <- GenomicFeatures::transcriptLengths(txdb)
match <- len[len$gene_id %in% genes,]
match <- match[match$tx_name %in% names(tx),]
txNames <- match$tx_name
longEnough <- filterTranscripts(txdb, 100, 100, 100)
txNames <- txNames[txNames %in% longEnough]

# Create coverage plots of miR430 targets
leaders = fiveUTRsByTranscript(txdb,use.names = T)[txNames]
cds <- cdsBy(txdb,"tx", use.names = TRUE)[txNames]
trailers = threeUTRsByTranscript(txdb, use.names = TRUE)[txNames]
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[txNames]
transcriptWindow(leaders, cds, trailers, outdir = p(mainFolder, "/tcp_plots/mir430/targets_"))

# Create coverage of non-miR430 targets
validNames <- longEnough
validNames <- validNames[!(validNames %in% txNames)]
leaders = fiveUTRsByTranscript(txdb,use.names = T)[validNames]
cds <- cdsBy(txdb,"tx", use.names = TRUE)[validNames]
trailers = threeUTRsByTranscript(txdb, use.names = TRUE)[validNames]
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[validNames]
transcriptWindow(leaders, cds, trailers)

# Create counts per library per transcipt
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[txNames] # targets
countsPerLibraryOverTranscript(tx)
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[validNames] # non-targets
countsPerLibraryOverTranscript(tx, output = p(mainFolder, "/tcp_plots/countsPerTranscript_Normal.csv"))

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

