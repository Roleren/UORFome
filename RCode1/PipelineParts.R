
getLeadersFromCage <- function(nCageList){
  foreach(i=1:nCageList) %dopar% {
    
    source("./uorfomeGeneratorHelperFunctions.R")
    cageList = list.files(cageFolder)
    getLeaders(usingNewCage = T, cageName = p(cageFolder,cageList[i]) , assignLeader = F )
    
  }
}

getUorfsFromLeaders <- function(nLeadersList){
  foreach(i=1:nLeadersList) %dopar% {
    source("./uorfomeGeneratorHelperFunctions.R")
    leadersList = list.files(leadersFolder)
    
    usedCage = gsub(pattern = ".leader.rdata",replacement = "",x = leadersList[i]) 
    
    if(UorfRangesNotExists(assignUorf =  F,givenCage = usedCage)){
      load(p(leadersFolder,leadersList[i]))
      scanUORFs(fiveUTRs, outputName = usedCage, assignUorf = F)
      print("ok")
    }else{
      i
    }
  }
}

getAllFeaturesFromUorfs <- function(){
  rfpList <- grep(pattern = "merged",x = list.files(rfpFolder), value = T)
  nrfpList <- length(rfpList)
  # all rfp features
  foreach(i=1:nrfpList) %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./DataBaseCreator.R")
    setwd("/export/valenfs/projects/uORFome/dataBase/")
    rfpList <- grep(pattern = "merged",
                    x = list.files(rfpFolder), value = T)
    RFPPath <- p(rfpFolder, rfpList[i])
    grl <- getUorfsInDb()
    
    getAllFeatures(grl,RFPPath, i = i)
    print(i)
  }
  
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./DataBaseCreator.R")
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  
  uorfID <- getORFNamesDB(T, F, F)
  floss <- uorfID
  entropyRFP <- uorfID        
  disengagementScores <- uorfID
  RRS <- uorfID                
  RSS <- uorfID           
  fpkmRFP <- uorfID           
  ORFScores <- uorfID
  ioScore <- uorfID
  #rm(list = c("floss","entropyRFP","disengagementScores", "RRS", "RSS", "fpkmRFP", "ORFScores", "ioScore" ))
  
  # all dt features, split them, save them seperatly
  foreach(i=1:nrfpList) %do% {
    
    dtFolder <- "featureTablesTemp/"
    dtPath <- paste0(dtFolder, "dt_",i,".rdata")
    
    load(dtPath)
    if(ncol(dt) != 8) stop(paste("not correct ncol of", i))
    floss[, p("floss_",i)] <- dt$floss
    entropyRFP[, p("entropyRFP_",i)] <- dt$entropyRFP
    disengagementScores[, p("disengagementScores_",i)] <- dt$disengagementScores 
    RRS[, p("RRS_",i)] <- dt$RRS
    RSS[, p("RSS_",i)] <- dt$RSS
    fpkmRFP[, p("fpkmRFP_",i)] <- dt$fpkmRFP
    ORFScores[, p("ORFScores_",i)] <- dt$ORFScores
    ioScore[, p("ioScore_",i)] <- dt$ioScore
    print(i)
  }
  # insert all the ribo features tables
  insertTable(floss, "floss")
  insertTable(entropyRFP, "entropyRFP")
  insertTable(disengagementScores, "disengagementScores")
  insertTable(RRS, "RRS")
  insertTable(RSS, "RSS")
  insertTable(fpkmRFP, "Ribofpkm")
  insertTable(ORFScores, "ORFScores")
  insertTable(ioScore, "ioScore")
  
  # insert info
  getRiboInfoTable(rfpList = rfpList)
}

#' Make RNA-seq fpkm values for database
#' 
#' Since rna-seq fpkms are normalized by transcript, the size of this
#' tabke might be different than the ribo-seq. i.g. there can be several
#' uORFs per transcript.
getRNAFpkms <- function(){
  load("/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/unfilteredSpeciesGroup.rdata")
  rnaList <- SpeciesGroup[SpeciesGroup$Sample_Type == "RNA",]$RnaRfpFolders
  rnaList <- grep(pattern = ".bam",x = rnaList, value = T)
  nrnaList <- length(rnaList)
  rm(SpeciesGroup)
  # all rfp features
  
  rnaFPKMs <- foreach(i=1:nrnaList, .combine = 'cbind') %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./DataBaseCreator.R")
    setwd("/export/valenfs/projects/uORFome/dataBase/")
    load("/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/unfilteredSpeciesGroup.rdata")
    rnaList <- SpeciesGroup[SpeciesGroup$Sample_Type == "RNA",]$RnaRfpFolders
    rnaList <- grep(pattern = ".bam",x = rnaList, value = T)
    RNAPath <- rnaList[i]
    
    # get tx
    tx <- getTx()
    tx <- ORFik:::extendLeaders(tx)
    RNA <- readGAlignments(RNAPath)
    rnaFPKM <- ORFik:::fpkm(tx, reads = RNA)
  }
  tx <- getTx()
  rnaFPKMs <- as.data.table(rnaFPKMs)
  txNames <- names(tx)
  rnaFPKMs <- data.table(txNames, rnaFPKMs)
  
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./DataBaseCreator.R")
  insertTable(rnaFPKMs, "RNAfpkm")
  getRNASeqInfoTable(rnaList = rnaList)
}

getTeFeatures <- function(riboDbName = "Ribofpkm",
                          dbOutputNames = c("teUnfiltered", "teFiltered")){
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  if(length(dbOutputNames) != 2) stop("dbOutputNames must have 2 character elements")
  
  # load linking and ribo / rna
  
  linking <- matchRNA_RFPInfo("linkRnaRfp")
  
  RFP <- readTable(riboDbName)
  txNames <- RFP$txNames
  RNA <- readTable("RNAfpkm")
  RNA <- matchByTranscript(RNA, RFP)
  RFP <- removeIDColumns(RFP)
  RNA <- removeIDColumns(RNA)
  
  if(nrow(RFP) != nrow(RNA)) stop("riboseq and rnaseq tables have different # of rows")
  
  # find number of linkings we have
  nTE <- max(linking$matching)
  
  library(foreach)
  # unfiltered without pseudoCounts
  teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
    rows <- linking[linking$matching == i, c(Sample_Type, originalIndex)]
    if(length(rows) != 4) stop("something wrong te with nrows")
    type <- rows[1:2]
    indices <- as.integer(rows[3:4])
    
    if ((type[1] == "RNA") && (type[2] == "RPF")) {
      return(RFP[,(indices[2]), with = F] / RNA[,indices[1], with = F])
    } else if ((type[1] == "RPF") && (type[2] == "RNA")) {
      return(RFP[,(indices[1]), with = F] / RNA[,indices[2], with = F])
    } else {
      stop("something wrong te with nrows")
    }
  }
  
  teTable <- data.table(txNames, teTable)
  insertTable(teTable, dbOutputNames[1])
  
  # filtered with pseudoCounts
  
  teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
    rows <- linking[linking$matching == i, c(Sample_Type, originalIndex)]
    if(length(rows) != 4) stop("something wrong te with nrows")
    type <- rows[1:2]
    indices <- as.integer(rows[3:4])
    
    if ((type[1] == "RNA") && (type[2] == "RPF")) {
      return( (RFP[,(indices[2]), with = F] + 1) / (RNA[,indices[1], with = F] + 1))
    } else if ((type[1] == "RPF") && (type[2] == "RNA")) {
      return((RFP[,(indices[1]), with = F] + 1)  / (RNA[,indices[2], with = F] + 1))
    } else {
      stop("something wrong te with nrows")
    }
  }
  teTable <- data.table(txNames, teTable)
  insertTable(teTable, dbOutputNames[2])
}

getCDSFeatures <- function(){
  ## Riboseq fpkm
  if(tableNotExists("cdsRfpFPKMs")) {
    rfpList <- grep(pattern = "merged",x = list.files(rfpFolder), value = T)
    nrfpList <- length(rfpList)
    # all rfp features
    cdsRfpFPKMs <- foreach(i=1:nrfpList, .combine = 'cbind') %dopar% {
      setwd("/export/valenfs/projects/uORFome/RCode1/")
      source("./DataBaseCreator.R")
      setwd("/export/valenfs/projects/uORFome/dataBase/")
      rfpList <- grep(pattern = "merged",
                      x = list.files(rfpFolder), value = T)
      RFPPath <- p(rfpFolder, rfpList[i])
      RFP <- ORFik:::cageFromFile(RFPPath)
      getCDS()
      
      rfps <- fpkm(cds, RFP)
    }
   
    getCDS()
    txNames <- names(cds)
    cdsRfpFPKMs <- as.data.table(cdsRfpFPKMs)
    cdsRfpFPKMs <- data.table(txNames, cdsRfpFPKMs)
    
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./DataBaseCreator.R")
    insertTable(cdsRfpFPKMs, "cdsRfpFPKMs")
    
    riboAtlasFPKMTissue(riboDbName = "cdsRfpFPKMs",
                        dbOutputNames = c("cdsRiboByTissueTF", "cdsRiboByTissueMean"))
  }
  
  ## RNA fpkms already made, so go to TE:
  if(tableNotExists("cdsTeFiltered")) {
    getTeFeatures(riboDbName = "cdsRfpFPKMs",
                  dbOutputNames = c("cdsTeUnfiltered", "cdsTeFiltered"))
  }
  
  
}