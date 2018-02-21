setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./uorfomeGeneratorHelperFunctions.R")
library(doParallel)

maxCores = 10
cl <- makeCluster(maxCores)
registerDoParallel(cl)

foreach(i=1:10) %dopar% {
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./DataBaseCreator.R")
  rfpDir <- "/export/valenfs/data/processed_data/Ribo-seq/fantom_human_bed/"
  rfps <- grep(x = list.files(rfpDir),
               pattern = "merged", value = T)
  
  if (i > 5) {
    rfps <- grep(x = rfps, pattern = "HEK293", value = T)[i-5]
    standardCage = "./../DATA/CAGE/human/kidney%2c%20adult%2c%20pool1.CNhs10622.10017-101C8.hg38.nobarcode.ctss.bed.gz"
    rnas <- list.files("/export/valenfs/data/processed_data/RNA-seq/andreev_DE_2015_human/final_results/aligned_GRCh38/")[i -5]
    rnas <- paste0("/export/valenfs/data/processed_data/RNA-seq/andreev_DE_2015_human/final_results/aligned_GRCh38/", rnas)
  } else {
    rfps <- grep(x = rfps, pattern = "brain", value = T)[i]
    rnas <- list.files("/export/valenfs/data/processed_data/RNA-seq/gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/")[i]
    rnas <- paste0("/export/valenfs/data/processed_data/RNA-seq/gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/", rnas)
  }
  rfps <- paste0(rfpDir, rfps)
  
  getCDS()
  getThreeUTRs()
  getLeaders()
  tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)
  
  cageFiveUTRs <- ORFik:::reassignTSSbyCage(fiveUTRs, standardCage, 1000, 1, cds)
  originalUorfsByTx <- getUnfilteredUORFs(cageFiveUTRs, assignRanges = F)
  gr <- unlist(originalUorfsByTx, use.names = F)
  grl <- groupGRangesBy(gr, gr$names)
  grl <- removeORFsWithinCDS(grl)
  tx <- ORFik:::extendLeaders(tx, extension = cageFiveUTRs)
 
  
  
  RFPShifted <- ORFik:::cageFromFile(rfps)
  RNA <- readGAlignments(rnas)
  
  dt <- ORFik:::allFeatures(grl = grl, RFP = RFPShifted, RNA = RNA,
                            fiveUTRs = fiveUTRs, cds = cds, tx = tx,
                            threeUTRs = threeUTRs, faFile = fa, riboStart = 26, riboStop = 34,
                            extension = 1000, orfFeatures = T,
                            cageFiveUTRs = cageFiveUTRs)
  if(i > 5)
    data.table::fwrite(dt, file = paste0("dt_HEK293", i-5, ".csv"))
  if(i < 6)
    data.table::fwrite(dt, file = paste0("dt_Brain", i, ".csv"))
}


stopCluster(cl)