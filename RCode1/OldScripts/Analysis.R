setwd("/export/valenfs/projects/uORFome/RCode1/") #!! set this path as codeFolder
source("./pipelineSetup.R")
source(p(codeFolder,"/DataBaseSetup.R"))

grl <- getUorfsInDb()
# getTx(T)
getCageTx()
getCDS()


teNames <- readTable("cdsTETissueMean")
te <- readTable("cdsTETissueMean", with.IDs = F)
RNA <- readTable("RNAfpkm")
rnaFilter <- RNA[which(rowMeans(RNA[,-1]) > 1),]$txNames
load(paste0("forests/finalPrediction_filtered","all", ".rdata"))
uorfTxAll <- ORFik:::txNames(grl)[prediction$predict == 1]
uorfTx <- unique(uorfTxAll)

#1. TE variance of uORF existence

teUorf <-te[teNames$txNames %in% uorfTx & (teNames$txNames %in% rnaFilter),]
teNot <- te[!(teNames$txNames %in% uorfTx) & (teNames$txNames %in% rnaFilter),]

dim(teUorf)
dim(teNot)
summary(rowMeans(teUorf))
summary(rowMeans(teNot))

#2. How number of uORFs changes TE
uorfTxDT <- data.table(uorfTxAll)
numberOfUorfs <- uorfTxDT[, .N, by = uorfTxAll]
groupingsUORFs <- numberOfUorfs[,  .N, by = N]
groupingsUORFs <- groupingsUORFs[order(groupingsUORFs$N)]
storeHits <- c()
for(i in sort(unique(numberOfUorfs$N))) {
  if(i > 11) break
  print(i)
  hits <- unique(numberOfUorfs$uorfTxAll[numberOfUorfs$N == i])
  sums <- summary(rowMeans(te[teNames$txNames %in% hits & (teNames$txNames %in% rnaFilter),]))
  print(sums)
  storeHits <- c(storeHits, sums)
}
