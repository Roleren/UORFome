
getLeadersFromCage <- function(nCageList){
  foreach(i=1:nCageList, .inorder = F) %dopar% {
    
    source("./uorfomeGeneratorHelperFunctions.R")
    cageList = grep(pattern = ".bed", list.files(cageFolder), value = TRUE)
    getLeaders(cageName = p(cageFolder,cageList[i]),
               assignLeader = F , exportUorfRegions = T)
    print("ok")
  }
}

getUorfsFromLeaders <- function(nLeadersList){
  foreach(i=1:nLeadersList, .inorder = F, .verbose = T, .export = c("FaFile", "extractTranscriptSeqs", "stopSites", "groupGRangesBy")) %dopar% {
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
}

getIDsFromUorfs <- function(nuorfsList){
  foreach(i=1:nuorfsList, .inorder = F) %dopar% {
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./uorfomeGeneratorHelperFunctions.R")
    load(p(uorfFolder, list.files(uorfFolder)[i]))
    
    uorfID <- unique(ORFik:::orfID(rangesOfuORFs))
    saveName = paste0(resultsFolder,"/uorfIDs/",gsub("uorf.rdata","",list.files(uorfFolder)[i]),"uorfID.rdata")
    save(uorfID, file = saveName)
    print("ok")
  }
}

getAllFeaturesFromUorfs <- function(){
  if (tableNotExists("ioScore")) {
    nrfpList <- nrow(getRiboRNAInfoTable())
    # all rfp features
    foreach(i=1:nrfpList, .inorder = F) %dopar% {
      setwd("/export/valenfs/projects/uORFome/RCode1/")
      source("./DataBaseSetup.R")
      
      matching_rna_ribo <- getRiboRNAInfoTable()
      
      RFPPath <- p(rfpFolder, grep(pattern = matching_rna_ribo$ribo[i] ,x = list.files(rfpFolder), value = T))
      if(length(RFPPath) != 1) stop(paste("did not find unique RFP file for:", matching_rna_ribo$study[i],matching_rna_ribo$ribo[i]))
      
      getAllFeatures(grl = NULL, RFPPath, i = i)
      print(i)
    }
    
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./DataBaseSetup.R")
    #linkORFsToTxUnique
    uorfID <- getORFNamesDB(T, F, F, T)
    floss <- copy(uorfID)
    entropyRFP <- copy(uorfID)   
    disengagementScores <- copy(uorfID)
    RRS <- copy(uorfID)               
    RSS <- copy(uorfID)        
    fpkmRFP <- copy(uorfID)     
    ORFScores <- copy(uorfID)
    ioScore <- copy(uorfID)
    
    
    # all dt features, split them, save them seperatly
    foreach(i=1:nrfpList) %do% {
      
      dtFolder <- "featureTablesTemp/"
      dtPath <- paste0(dtFolder, "dt_",i,".rdata")
      
      load(dtPath)
      if(ncol(dt) != 8) stop(paste("not correct ncol of", i))
      floss[, p("floss_",i) := dt$floss]
      entropyRFP[, p("entropyRFP_",i) := dt$entropyRFP]
      disengagementScores[, p("disengagementScores_",i) := dt$disengagementScores ]
      RRS[, p("RRS_",i) := dt$RRS]
      RSS[, p("RSS_",i) := dt$RSS]
      fpkmRFP[, p("fpkmRFP_",i) := dt$fpkmRFP] 
      ORFScores[, p("ORFScores_",i) := dt$ORFScores]
      ioScore[, p("ioScore_",i) := dt$ioScore]
      print(i)
    }
    # insert all the ribo features tables
    insertTable(floss, "floss", rmOld = T)
    insertTable(entropyRFP, "entropyRFP", rmOld = T)
    insertTable(disengagementScores, "disengagementScores", rmOld = T)
    insertTable(RRS, "RRS", rmOld = T)
    insertTable(RSS, "RSS", rmOld = T)
    insertTable(fpkmRFP, "Ribofpkm", rmOld = T)
    insertTable(ORFScores, "ORFScores", rmOld = T)
    insertTable(ioScore, "ioScore", rmOld = T)
  } else {
    print("AllFeaturesFromUorfs exists in DB (ioScore), delete and run again if you want new")
  }
  
  # insert info
  #getRiboInfoTable(rfpList = rfpList)
  #rm(list = c("floss","entropyRFP","disengagementScores", "RRS",
  #"RSS", "fpkmRFP", "ORFScores", "ioScore" ))
}

#' Make RNA-seq fpkm values for database
#' 
#' Since rna-seq fpkms are normalized by transcript, the size of this
#' tabke might be different than the ribo-seq. i.g. there can be several
#' uORFs per transcript.
getRNAFpkms <- function(){
  if (!tableNotExists("RNAfpkm")) {
    print("RNAfpkm table exists, remove it if you want to rerun with new values")
    return(NULL)
  }
  matching_rna_ribo <- getRiboRNAInfoTable()
  nrnaList <- nrow(matching_rna_ribo)
  # all rfp features
  
  rnaFPKMs <- foreach(i=1:nrnaList, .combine = 'cbind', .inorder = T, .export = c("getCageTx"), .packages = c("GenomicFeatures", "ORFik")) %dopar% {
    load(file = "/export/valenfs/projects/uORFome/matching_rna_ribo.rdata")
    SRR <- matching_rna_ribo$rna[i]
    RNAPath <- grep(SRR,list.files(path = rnaFolder, all.files = TRUE, full.names = TRUE, recursive = TRUE), value=TRUE) 
    if (length(RNAPath) != 1) stop("did not find unique RNA path for SRR")
    # get tx
    getCageTx()
    RNA <- readGAlignments(RNAPath)
    rnaFPKM <- ORFik:::fpkm(tx, reads = RNA)
  }
  getCageTx()
  rnaFPKMs <- as.data.table(rnaFPKMs)
  txNames <- names(tx)
  rnaFPKMs <- data.table(txNames, rnaFPKMs)
  
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./DataBaseSetup.R")
  insertTable(rnaFPKMs, "RNAfpkm")
  return(NULL)
}

#' This function uses the fact that 1st col of ribo is connected to 1st col of RNA. 
getTeFeatures <- function(riboDbName = "Ribofpkm",
                          dbOutputNames = c("teUnfiltered", "teFiltered")){
  setwd("/export/valenfs/projects/uORFome/dataBase/")
  if(length(dbOutputNames) != 2) stop("dbOutputNames must have 2 character elements")
  
  # load linking and ribo / rna
  
  RFP <- readTable(riboDbName)
  txNames <- RFP$txNames
  RNA <- readTable("RNAfpkm")
  RNA <- matchByTranscript(RNA, RFP)
  RFP <- removeIDColumns(RFP)
  RNA <- removeIDColumns(RNA)
  
  if(nrow(RFP) != nrow(RNA)) stop("riboseq and rnaseq tables have different # of rows")
  if(ncol(RFP) != ncol(RNA)) stop("riboseq and rnaseq tables have different # of cols")
  # find number of linkings we have
  nTE <- ncol(RFP)
  
  library(foreach)
  # unfiltered without pseudoCounts
  teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
    print(i)
    return(RFP[,i, with = F] / RNA[,i, with = F])
  }
  
  teTable <- data.table(txNames, teTable)
  insertTable(teTable, dbOutputNames[1])
  
  # filtered with pseudoCounts
  
  teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
    print(i)
    return((RFP[,i, with = F] + 1) / (RNA[,i, with = F] + 1))
  }
  teTable <- data.table(txNames, teTable)
  insertTable(teTable, dbOutputNames[2])
}

getCDSFeatures <- function(){
  if(tableNotExists("cdsKozak")) {
    getCDS()
    getFasta()
    cdsKozak <- ORFik:::kozakSequenceScore(cds, fa)
    
    txNames <- names(cds)
    
    cdsKozak <- data.table(txNames, cdsKozak)
    
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./DataBaseSetup.R")
    insertTable(cdsKozak, "cdsKozak")
  }
  
  if(tableNotExists("cdsIos")) {
    
    rfpList <- grep(pattern = "merged",x = list.files(rfpFolder), value = T)
    nrfpList <- length(rfpList)
    getCDS()
    gr <- unlist(cds, use.names = T)
    cds <- groupGRangesBy(gr)
    getThreeUTRs()
    gr <- unlist(threeUTRs, use.names = T)
    threeUTRs <- groupGRangesBy(gr)
    rm(tx)
    tx <- getTx()
    tx <- tx[names(cds)]
    cageFiveUTRs <- leaderCage()
    tx[names(cageFiveUTRs)] <- ORFik:::extendLeaders(tx[names(cageFiveUTRs)], cageFiveUTRs)
    
    allRiboFeatures <- foreach(i=1:nrfpList, .combine = 'cbind', .export = c("cds", "threeUTRs", "tx")) %dopar% {
      setwd("/export/valenfs/projects/uORFome/RCode1/")
      source("./DataBaseSetup.R")
      
      rfpList <- grep(pattern = "merged",
                      x = list.files(rfpFolder), value = T)
      RFPPath <- p(rfpFolder, rfpList[i])
      RFP <- fread.bed(RFPPath)
      
      cdsDiseng <- disengagementScore(cds, RFP, tx)
      cdsIos <- insideOutsideORF(cds, RFP, tx, ds = cdsDiseng)
      cdsRRS <- ribosomeReleaseScore(cds, RFP, threeUTRs)
      cdsRSS <- ribosomeStallingScore(cds, RFP)
      cdsFloss <- floss(cds, RFP, cds)
      cdsEntropy <- entropy(cds, RFP)
      cdsORFScores <- orfScore(cds, RFP, T)$ORFScores
      cdsRfpFPKMs <- fpkm(cds, RFP)
      return(data.table(cdsDiseng,cdsIos, cdsRRS, cdsRSS, cdsORFScores, cdsFloss, cdsEntropy, cdsRfpFPKMs))
    }
    txNames <- names(cds)
    
    for(f in unique(colnames(allRiboFeatures))) {
      featu <-  allRiboFeatures[, which(colnames(allRiboFeatures) == f), with = F]
      colnames(featu) <- paste0(f,"_", 1:ncol(featu))
      featu <- data.table(txNames, featu)
      insertTable(featu, f)
    }
    
  }
  
  if(tableNotExists("cdsFractionLengths")) {
    getCDS()
    tx <- getTx()
    tx <- tx[names(cds)]
    cageFiveUTRs <- leaderCage()
    tx[names(cageFiveUTRs)] <- ORFik:::extendLeaders(tx[names(cageFiveUTRs)], cageFiveUTRs)
    cdsFrac <- fractionLength(cds, widthPerGroup(tx))
    
    txNames <- names(cds)
    
    cdsFrac <- data.table(txNames, cdsFrac)
    
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./DataBaseSetup.R")
    insertTable(cdsFrac, "cdsFractionLengths")
  }
  
  ## RNA fpkms already made, so go to TE:
  # if(tableNotExists("cdsTeFiltered")) {
  #   riboAtlasFPKMTissue(riboDbName = "cdsRfpFPKMs",
  #                       dbOutputNames = c("cdsRiboByTissueTF", "cdsRiboByTissueMean"))
  #   getTeFeatures(riboDbName = "cdsRfpFPKMs",
  #                 dbOutputNames = c("cdsTeUnfiltered", "cdsTeFiltered"))
  #   cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  #   cdsTETissue <- teAtlasTissueNew(cdsTEs, colExclusion = "result.")
  #   insertTable(cdsTETissue, "cdsTETissueMean")
  # }
}

#' downstream of three utrs is used as new 3'utre
getFeaturesThreeUTRs <- function(){
  if(tableNotExists("threeIos")) {
    
    getCDS()
    gr <- unlist(cds, use.names = T)
    cds <- groupGRangesBy(gr)
    getThreeUTRs()
    gr <- unlist(threeUTRs, use.names = T)
    threeUTRs <- groupGRangesBy(gr)
    threeUTRs <- threeUTRs[widthPerGroup(threeUTRs) > 5]
    threeWidth <- median(widthPerGroup(threeUTRs))
    fakeThree <- GRanges(seqnamesPerGroup(threeUTRs, F), IRanges(stopSites(threeUTRs, is.sorted = T) + 1, width = threeWidth), strand = strandPerGroup(threeUTRs, F))
    names(fakeThree) <- names(threeUTRs)
    fakeThree <- groupGRangesBy(fakeThree)
    rm(tx)
    tx <- getTx()
    tx <- tx[names(threeUTRs)]
    cageFiveUTRs <- leaderCage()
    matchName <- names(threeUTRs)[names(threeUTRs)%in% names(cageFiveUTRs)]
    tx[matchName] <- ORFik:::extendLeaders(tx[matchName], cageFiveUTRs[matchName])
    rfpList <- grep(pattern = "merged",x = list.files(rfpFolder), value = T)
    nrfpList <- length(rfpList)
    
    allRiboFeatures <- foreach(i=1:nrfpList, .combine = 'cbind', .export = c("cds", "threeUTRs", "tx", "fakeThree", "rfpList")) %dopar% {
      setwd("/export/valenfs/projects/uORFome/RCode1/")
      source("./DataBaseSetup.R")
      
      RFPPath <- p(rfpFolder, rfpList[i])
      RFP <- fread.bed(RFPPath)

      Diseng <- disengagementScore(threeUTRs, RFP, tx)
      Ios <- insideOutsideORF(threeUTRs, RFP, tx, ds = Diseng)
      RRS <- ribosomeReleaseScore(threeUTRs, RFP, fakeThree)
      RSS <- ribosomeStallingScore(threeUTRs, RFP)
      floss <- floss(threeUTRs, RFP, cds)
      Entropy <- entropy(threeUTRs, RFP)
      ORFScores <- orfScore(threeUTRs, RFP, T)$ORFScores
      RfpFPKMs <- fpkm(threeUTRs, RFP)
      data.table(Diseng, Ios, RRS, RSS, ORFScores, floss, Entropy, RfpFPKMs)
    }
    colnames(allRiboFeatures) <- paste0("three", colnames(allRiboFeatures))
    txNames <- names(threeUTRs)
    
    for(f in unique(colnames(allRiboFeatures))) {
      featu <-  allRiboFeatures[, which(colnames(allRiboFeatures) == f), with = F]
      colnames(featu) <- paste0(f,"_", 1:ncol(featu))
      featu <- data.table(txNames, featu)
      insertTable(featu, f)
    }
    
  }
  
  # if (tableNotExists("threeTeFiltered")) {
  #   getTeFeatures(riboDbName = "threeFpkm", dbOutputNames =  c("threeTeUnfiltered", "threeTeFiltered"))
  # }
  
  # sequence features
  if (tableNotExists("threeKozak")) {
    getThreeUTRs()
    gr <- unlistGrl(threeUTRs)
    threeUTRs <- groupGRangesBy(gr)
    threeUTRs <- threeUTRs[widthPerGroup(threeUTRs) > 5]
    getFasta()
    Kozak <- kozakSequenceScore(threeUTRs, fa)
    insertTable(Kozak, "threeKozak")
  }
  
}