arcs = commandArgs(trailingOnly = T)
source("./uorfomeGeneratorHelperFunctions.R")


leaderFromBash = getRelativePathName(arcs[1])
print(arcs[1])

load(p(leadersFolder,leaderFromBash))
usedCage = gsub(pattern = ".leader.rdata",replacement = "",x = leaderFromBash)
cat("cage file now before is:\n",usedCage)
decideHowToGetUORFRanges(assignUorf = T,givenCage = usedCage)


