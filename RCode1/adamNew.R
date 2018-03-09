rm(list=ls())
#.libPaths(c("/Home/ii/hakontj/R/x86_64-redhat-linux-gnu-library/3.3","/usr/lib64/R/library","/usr/share/R/library" ))
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./uorfomeGeneratorHelperFunctions.R")
library(doParallel)

# list of leaders
leadersList <- c("/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_03_64_cells.gff3"
                 ,"/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_04_512_cells.gff3"
                 ,"/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_07_sphere.gff3"
                 ,"/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_09_shield.gff3")
nLeadersList = length(leadersList)

# set cores to min of list and half max
maxCores = min(as.integer(detectCores()/2),as.integer(nLeadersList))  
cl <- makeCluster(maxCores)
registerDoParallel(cl)


#if libPaths are strange clusterCall(cl, function(x) .libPaths(x), .libPaths())

foreach(i=1:nLeadersList) %dopar% {
  source("./uorfomeGeneratorHelperFunctions.R")
  leadersList <- c("/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_03_64_cells.gff3"
                   ,"/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_04_512_cells.gff3"
                   ,"/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_07_sphere.gff3"
                   ,"/export/valenfs/projects/adam/leaders_Zv10/leaders_Jan_2018/GFF3s/zebrafish_09_shield.gff3")
  uorfFolder <- "/export/valenfs/projects/uORFome/adamUorfs/"
  getFasta("/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.dna.toplevel.fa")
  # Gtf <- makeTxDbFromGFF(file = leadersList[i], format = "gff3")
  # 
  # fiveUTRs <- GenomicFeatures::fiveUTRsByTranscript(Gtf, use.names = T)
  
  gff <- import.gff3(leadersList[i])
  leaders <-gff[gff$type == "five_prime_UTR"]
  leaders$score <- NULL
  leaders$phase <- NULL
  leaders$ID <- NULL
  leaders$type <- NULL
  names(leaders) <-  unlist(leaders$Parent)
  leaders$Parent <- NULL
  fiveUTRs <- groupGRangesBy(leaders)
  
  scanUORFs(fiveUTRs, outputName = leadersList[i], assignUorf = F, 
            outputFastaAndBed = F, filterORFs = F)
  
  print(paste("ok", i))
  
}

stopCluster(cl)
