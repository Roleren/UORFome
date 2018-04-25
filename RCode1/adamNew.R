rm(list=ls())


# this is adams first csv reader for orfik
otherSpeciesCageCSV = function(name){
  #pre loadings
  #name = "/export/valenfs/data/processed_data/CAGE/nepal_2013_zebrafish/called_peaks/leaders_zv10_zf_02_fertilized_egg.csv"
  cage = read.csv(name)
  cage = cage[cage$count_at_highest_peak != 0,]
  
  fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
  
  #make max peak object
  cageGR = GRanges(seqnames = gsub("chr","",as.character(cage$chr)), ranges = IRanges(start = cage$highest_peak,end = cage$highest_peak),strand = cage$dir, names = cage$X.gene_id)
  maxPeakPosition = as.data.table(findOverlaps(query = cageGR,subject = extendsTSSExons(fiveUTRs)))
  names(maxPeakPosition) = c("from","to")
  maxPeakPosition = maxPeakPosition[!duplicated(maxPeakPosition$to)]
  maxPeakPosition$strand = as.character(strand(cageGR[maxPeakPosition$from]))
  maxPeakPosition$start = start(cageGR[maxPeakPosition$from])
  maxPeakPosition$end = end(cageGR[maxPeakPosition$from])
  assign("maxPeakPosition",maxPeakPosition,envir = .GlobalEnv)
  
  #make new fiveUTRs
  fiveUTRs = makeGrlAndFilter(addNewTssOnLeaders(fiveUTRs), fiveUTRs)
  fiveUTRs = addFirstCdsOnLeaderEnds(fiveUTRs)
  
  unlistfgr = unlist(fiveUTRs)
  #On zebra fish, the genome and gtf had different seqname naming, so must match!
  seqnamesTransformed = as.character(seqnames(unlistfgr))
  seqnamesTransformed[nchar(seqnamesTransformed) > 4] = paste0("Un_",seqnamesTransformed[nchar(seqnamesTransformed) > 4])
  seqnamesTransformed = gsub("\\.","v",seqnamesTransformed)
  seqnamesTransformed = paste0("chr",seqnamesTransformed)
  unlistfgr = GRanges(seqnames = seqnamesTransformed, ranges = IRanges(start = start(unlistfgr),end = end(unlistfgr)),strand = strand(unlistfgr))
  names(unlistfgr) = names(unlist(fiveUTRs))
  fiveUTRs = groupGRangesBy(unlistfgr)
  
  return(fiveUTRs)
}

# this is adams second version, used for zebrafish
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
