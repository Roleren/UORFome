setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path as codeFolder
source("./DataBaseSetup.R") # <- load uORFome pipeline
nmdFolder <- p(mainFolder,"/nmdProject/")
# Data
faiName <- "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa"
fa <- FaFile(file = faiName)
txdb <- GenomicFeatures::makeTxDbFromGFF("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81_chr.gtf",                                         format = "gtf")
cds <- cdsBy(txdb,"tx", use.names = TRUE)
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)

# load their nmd data
library(xlsx)
df2 <- xlsx::read.xlsx2(p(nmdFolder, "TCT_SuppTable.xlsx"), sheetName = "Table S3", startRow = 5,
                        colClasses = c("character", "character", "double", "double", "double",
                                       "double", "double", "character"))
big <- df2[(df2$YC.1 > df2$YC),]
bigger <- df2[(df2$YC.1 > df2$YC) & df2$p.value < 0.06, ]
theirHigh <- length(unique(bigger$EnsemblID))
theirAll <- length(unique(big$EnsemblID))

# ribo-seq folder
rfpFolder <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/shoelaces_GRCz10"
rfpFiles <- list.files(rfpFolder)
rfpFiles <- rfpFiles[grep(pattern = ".wig.gz", rfpFiles)]
rfpFiles <- paste0(rfpFolder,"/", rfpFiles)

# load ribo-seq
for(i in seq.int(1, length(rfpFiles), by = 2)){
  rfpForward <- rfpFiles[i]
  stage <- sub("_trimmed",x = sub(pattern = "\\..*",x = sub(".*GRCz10/",x = rfpForward,""), ""), "")
  print(p("stage: ", stage))
  rfp <- import.wig(rfpForward) # forward direction
  strand(rfp) <- "+"
  rfp2 <- import.wig(rfpFiles[i+1]) # backward direction
  strand(rfp2) <- "-"
  rfp <- sort(c(rfp, rfp2))
  rfp <- rep(rfp, rfp$score)
  rfp$score <- NULL
  
  df <- predictNMD(cds, rfp, fa)
  # View and save
  # View(goodHits)
  save(df, file = paste0(nmdFolder, "potentialOutOfFrameStopCodons_",stage,".Rdata"))
  
  # Statistics
  goodHits <- df[df$outOfFrame == T & df$hasStopCodon == T & df$distanceCDSEnd > 4, ]
  hits <- goodHits[goodHits$transcript %in% bigger$EnsemblID,]
  hitsBig <- goodHits[goodHits$transcript %in% big$EnsemblID,]
  overlapHigh <- length(unique(hits$transcript))
  overlapAll <- length(unique(hitsBig$transcript))
  
  print(paste("found", nrow(hits), "high confidence stop codon sites"))
  print(paste(round(overlapHigh/theirHigh, 2)*100, "% of high confidence set"))
  print(paste(round(overlapAll/theirAll, 2)*100, "% of whole set"))
  print(paste("found", nrow(goodHits), "candidate out of frame stop sites,",
              length(unique(goodHits$transcript)), "unique transcript candidates"))
  
  #lower <- df2[(df2$YC.1 < df2$YC) & df2$p.value < 0.06, ]
  overlapAllNext <- tail(overlapAll, 1) + 1
  stop = c(1:overlapAll,overlapAllNext:(overlapAllNext+length(unique(goodHits$transcript))-overlapAll-1))
  stopNext <- tail(stop, 1) + 1
  nmd = c(1:overlapAll,stopNext:(stopNext+theirAll-overlapAll-1))
  overlaps = c(1:overlapHigh, stopNext:(stopNext+theirHigh-overlapHigh-1))
  
  plotNMD(stop, nmd, overlaps, stage)
}


predictNMD <- function(cds, rfp, fa){
  # algorithm
  cov <- ORFik:::coverageByWindow(rfp, cds, is.sorted = T)
  cov <- cov[sum(cov) > 1]
  print(paste(round(length(cov) / length(cds), 2), "% of CDS had reads accepted by filter"))
  print(summary(sum(cov)))
  
  hits <- runValue(cov) > pmax(median(sum(cov)), 1) # median is 37 on average
  cov <- cov[any(hits)]
  hits <- hits[any(hits)]
  
  # now we know which are hits, now lets find their positions.
  len <- BiocGenerics::lengths(cov)
  sums <- cumsum(runLength(cov))
  locations <- sums[hits]
  names <-  names(unlist(locations, use.names = T))
  grl <- ORFik:::pmapFromTranscriptF(IRanges(start = unlist(locations, use.names = F), width = 1), cds[names],
                                     indices = seq.int(length(names)))
  grl <- ORFik:::removeMetaCols(grl)
  gr <- unlistGrl(grl) # <- here are the positions in genomic coordinates
  
  # Now lets make a window
  window <- ORFik:::windowPerGroup(gr = gr, tx = tx, downstream = 3, upstream = 3) # window of 3 each direction from point
  seqs <- ORFik:::txSeqsFromFa(grl = window, faFile = fa, is.sorted = T) # window per sequence
  
  taaHits <- grep(x = seqs, pattern = "TAA")
  tagHits <- grep(x = seqs, pattern = "TAG")
  tgaHits <- grep(x = seqs, pattern = "TGA")
  
  df <- data.frame(transcript = names, gr)
  df$end <- NULL
  df$distanceCDS <- unlist(locations, use.names = F)
  df$frame <- ORFik:::isInFrame(df$distanceCDS)
  df$outOfFrame <- df$frame != 0
  df$distanceCDSEnd <- widthPerGroup(cds[names]) - df$distanceCDS
  df$seqs <- seqs
  df$hasStopCodon <- seq.int(length(gr)) %in% unique(c(taaHits, tagHits, tgaHits))
  return(df)
}

plotNMD <- function(stop, nmd, overlaps, stage){
  library(VennDiagram)
  pdf(p(nmdFolder, paste0("nmdVenn_", stage,".pdf")));
  venn.plot <- venn.diagram(
    x = list (
      stop, nmd, overlaps
    ),
    filename = NULL,
    main = "NMD vs stop codon peaks",
    sub = "For nmd increasing with ribo seq",
    category = c("stop codon peaks", "NMD", "significant"),
    euler.d = TRUE
  );
  grid.draw(venn.plot);
  dev.off();
  
}