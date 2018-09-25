# rm(list=ls())
#.libPaths(c("/Home/ii/hakontj/R/x86_64-redhat-linux-gnu-library/3.3","/usr/lib64/R/library","/usr/share/R/library" ))

setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path
source("./DataBaseSetup.R")
setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path


# set up multithreading options
pipelineCluster(4) #!! set number of cores, I use 75 usually on furu.

searchRegionList = list.files(regionUORFs)
nLeadersList = length(searchRegionList)
# clusterCall(cl, function(x) .libPaths(x), .libPaths())
out <- c()
out <- foreach(i=1:nLeadersList, .inorder = F, .verbose = T, .export = c("FaFile", "extractTranscriptSeqs", "stopSites", "groupGRangesBy")) %dopar% {
  usedCage = gsub(pattern = ".regionUORF.rdata",replacement = "",x = searchRegionList[i])
  saveName <- getUORFRDataName(usedCage)
  if(!file.exists(saveName)){
    load(p(regionUORFs,searchRegionList[i]))
    rangesOfuORFs <- getUnfilteredUORFs(uORFSeachRegion,assignRanges = F, isSorted = T,
                                    startCodons = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT")
    rangesOfuORFs <- filterORFs(rangesOfuORFs)
    save(rangesOfuORFs, file = saveName)
    return(i)
  }
  return(0)
}

stopCluster(cl)
# devtools::install_github("JokingHero/ORFik", ref = "orfFinderSpeedUp")
# i = 1800
# load(p(regionUORFs,searchRegionList[i]))
# 
# uORFSeachRegion <- uORFSeachRegion[1:10000]
# rangesOfuORFs <- getUnfilteredUORFs(uORFSeachRegion,assignRanges = F, isSorted = T,
#                                     startCodons = "ATG")
# # rangesOfuORFs <- getUnfilteredUORFs(uORFSeachRegion,assignRanges = F, isSorted = T,
# #                                     startCodons = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT")
# getSequencesFromFasta(uORFSeachRegion, isSorted = T)
# result <- ORFik:::orfs_as_List(fastaSeqs = as.character(seqs, use.names = FALSE), 
#                        startCodon = "ATG", stopCodon = stopDefinition(1), longestORF = F, 
#                        minimumLength = 2)
# 
# # rangesOfuORFs = ORFik:::mapToGRanges(uORFSeachRegion, result)
# grl <- uORFSeachRegion
# grl <- sortPerGroup(grl, ignore.strand = FALSE)
# ranges = IRanges(start = unlist(result$orf[1], use.names = FALSE), 
#                  end = unlist(result$orf[2], use.names = FALSE))
# genomicCoordinates <- ORFik:::pmapFromTranscriptF(ranges, grl, result$index)
# 
# ranks <- ORFik:::makeExonRanks(genomicCoordinates, byTranscript = TRUE)
# asGR <- unlistGrl(genomicCoordinates)
# mcols(x = asGR) <- DataFrame(row.names = names(asGR), names = paste0(names(asGR), "_", ranks))
# groupGRangesBy(asGR, asGR$names)
# 
# 
# names <- names(grl)
# names(grl) <- NULL
# genomicCoordinates <- pmapFromTranscripts(x = ranges,
#                                           transcripts = grl[result$index])