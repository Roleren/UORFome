#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Test of timings of uORFome pipeline, latest update March 2020.

# ORFik timings:
#   
#   old uORFs in database:
#   system.time(startRegion <- windowPerGroup(startSites(grl, T, T, T), tx, -3, 9))
# user  system elapsed 
# 113.047  18.654 131.648 

if (requireNamespace("uORFomePipe")) {
  library(uORFomePipe)
} else stop("You do not have uORFomePipe installed")

setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path as codeFolder
source("./uORFome/Init_Variables.R") # Make this into experiments instead

dataBaseFolder <- p(mainFolder,"/dataBase")
grl = getUorfsInDb()
#gg <- grl
#grl <- startSites(grl, asGR = T, keep.names = T)
getAll()
startRegion <- startRegion(grl, tx, T, -3, 9)
rfps <- list.files("/export/valenfs/data/processed_data/Ribo-seq/All_human_bedo/", full.names = TRUE)[1:2]
RFPPath <- rfps[1]
RFP <- fimport(RFPPath)
Rprof()
out <- ORFik:::allFeaturesHelper(grl, RFP, RNA = NULL, tx, fiveUTRs, cds , threeUTRs,
                          faFile = NULL, riboStart = 26, riboStop = 34,
                          sequenceFeatures = FALSE, grl.is.sorted = TRUE,
                          weight.RFP = "score", weight.RNA = 1L,
                          st = startRegion)
Rprof(NULL)
system.time(out <- ORFik:::allFeaturesHelper(grl, RFP, RNA = NULL, tx, fiveUTRs, cds , threeUTRs,
                                             faFile = NULL, riboStart = 26, riboStop = 34,
                                             sequenceFeatures = FALSE, grl.is.sorted = TRUE,
                                             weight.RFP = "score", weight.RNA = 1L,
                                             st = startRegion))
# TImings:
system.time(ORFik:::optimizeReads(tx, RFP))
# user  system elapsed 
# 0.563   0.096   0.680 
system.time(floss(grl, RFP, cds, 26, 34, weight = "score"))
# user  system elapsed 
# 8.108   1.595   1.250 
system.time(orfScore(grl, RFP, weight = "score", is.sorted = T)) 
# user  system elapsed 
# 38.409   7.496   8.549 
system.time(entropy(grl, RFP, weight = "score",is.sorted = T))
# user  system elapsed 
# 100.187   4.725   9.667 
system.time(ds <- disengagementScore(grl, RFP, tx, weight = "score", RFP.sorted = T))
# user  system elapsed 
# 30.497   2.416  32.911
system.time(insideOutsideORF(grl, RFP, tx, weight = "score", RFP.sorted = T, ds = ds))
# user  system elapsed 
# 24.546   0.770  21.061 
system.time(ORFik:::startRegionCoverage(grl, RFP, tx, weight = "score"))
# user  system elapsed 
# 22.446   3.050  25.538 

# total: 70
# test

bed122 <- function(grl, bedName, fixChromoNaming = FALSE){
  # TODO, check seqlevels style should be forced to UCSC ?
  if(!ORFik:::is.grl(class(grl))) stop("grl, must be of class GRangesList")
  if (fixChromoNaming) print(seqlevels(grl))
  grl <- sortPerGroup(grl,ignore.strand = T) # <- sort bed way!
  
  dt.grl <- data.table(seqnames = ORFik:::seqnamesPerGroup(grl, F))
  dt.grl$start <- as.integer(ORFik:::firstStartPerGroup(grl,keep.names = F) -1)
  dt.grl$end <- ORFik:::lastExonEndPerGroup(grl, keep.names = F) #non inclusive end
  dt.grl$name <- names(grl)
  dt.grl$score <- widthPerGroup(grl, keep.names = F)
  dt.grl$strand <- strandPerGroup(grl, F)
  dt.grl$thickStart <- dt.grl$start
  dt.grl$thickEnd <- dt.grl$end
  dt.grl$rgb <- rep(0, length(grl))
  dt.grl$blockCount <- ORFik:::numExonsPerGroup(grl)
  blockSizes <- paste(width(grl), collapse = ",")
  names(blockSizes) <- NULL
  dt.grl$blockSizes <- blockSizes
  relativeStarts <- (start(grl) -1) - dt.grl$start
  blockStarts <- paste(relativeStarts, collapse = ",")
  names(blockStarts) <- NULL
  dt.grl$blockStarts <- blockStarts
  
  #chromStart + chromStarts[last] + blockSizes[last])
  #must equal chromEnd.
  data.table::fwrite(x = dt.grl, file = bedName,
                     sep = "\t", col.names = F, row.names = F, quote = F)
  return(NULL)
}

#
bplapply(c("1", "2"), function(x) {x})


grl <- gg
a <- pmapToTranscripts(grl[1:1000], grl[1:1000]);a
b <- pmapToTranscriptF(grl[1:1000], gg[1:1000]);b

identical(a, b)

inter <- 1:1000; print(1000)
system.time(pmapToTranscripts(grl[inter], gg[inter]))
system.time(pmapToTranscriptF(grl[inter], gg[inter]))      

inter <- 1:10000; print(10000)
system.time(pmapToTranscripts(grl[inter], grl[inter]))
system.time(pmapToTranscriptF(grl[inter], grl[inter]))      

inter <- 1:100000; print(100000)
system.time(pmapToTranscripts(grl[inter], grl[inter]))
system.time(pmapToTranscriptF(grl[inter], grl[inter]))      

inter <- 1:1000000; print(1000000)
system.time(pmapToTranscripts(grl[inter], grl[inter]))
system.time(pmapToTranscriptF(grl[inter], grl[inter]))

print("all")
Rprof()
pmapToTranscripts(grl, grl)
pmapToTranscriptF(grl, grl)
Rprof(NULL)
system.time(pmapToTranscripts(grl, grl))
system.time(pmapToTranscriptF(grl, grl)) 

inter <- 1:1000000; print(1000000)
system.time(a <- pmapToTranscriptF(startSites(grl[inter], asGR = T), grl[inter]))
system.time(b <- pmapToTranscripts(startSites(grl[inter], asGR = T), grl[inter]))
identical(a,b)


df.rfp <- read.experimentl("uORFome_rfp")
RFP <- fimport("/export/valenfs/data/processed_data/Ribo-seq/All_human_bedo/Andreev_DE_2015.Human.HEK293.RPF.GRCh38.SRR1173909.reads_merged.bedo")

startRegion(grl, tx, T, -3, 9)
