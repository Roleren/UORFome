setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./DataBaseCreator.R")
setwd("/export/valenfs/projects/uORFome/dataBase/")

# load linking and ribo / rna

linking <- matchRNA_RFPInfo("linkRnaRfp")

RFP <- readTable("Ribofpkm")
RNA <- readTable("RNAfpkm")

# find number of linkings we have
nTE <- max(linking$matching)

# unfiltered without pseudoCounts
teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
  rows <- linking[linking$matching == i, c(Sample_Type, originalIndex)]
  if(length(rows) != 4) stop("something wrong te with nrows")
  type <- rows[1:2]
  indices <- as.integer(rows[3:4])
  
  if ((type[1] == "RNA") && (type[2] == "RPF")) {
    return(RFP[,(indices[2]+2), with = F] / RNA[,indices[1], with = F])
  } else if ((type[1] == "RPF") && (type[2] == "RNA")) {
    return(RFP[,(indices[1]+2), with = F] / RNA[,indices[2], with = F]) # + uorf id columns
  } else {
      stop("something wrong te with nrows")
  }
}


insertTable(teTable, "teUnfiltered")

# filtered with pseudoCounts

teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
  rows <- linking[linking$matching == i, c(Sample_Type, originalIndex)]
  if(length(rows) != 4) stop("something wrong te with nrows")
  type <- rows[1:2]
  indices <- as.integer(rows[3:4])
  
  if ((type[1] == "RNA") && (type[2] == "RPF")) {
    return( (RFP[,(indices[2]+2), with = F] + 1) / (RNA[,indices[1], with = F] + 1))
  } else if ((type[1] == "RPF") && (type[2] == "RNA")) {
    return((RFP[,(indices[1]+2), with = F] + 1)  / (RNA[,indices[2], with = F] + 1)) # + uorf id columns
  } else {
    stop("something wrong te with nrows")
  }
}

insertTable(teTable, "teFiltered")


