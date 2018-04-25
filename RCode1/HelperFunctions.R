


#as number coverter function
an = function(fac){
  return(as.numeric(as.character(fac)))
}

#check that bamFile is paired or not, False = not paired
testBAM = function(name){
  testPairedEndBam(file = name,index = name)
}

#' remove non finite in data.table
#' 
#' Removes inf, NA, and NaN values
#' @param DT a data.table
#' @param replacement (0) what should the non finite values be changed to ?
removeNonFinite <- function(DT, replacement = 0){
  invisible(lapply(names(DT),function(.name)
    set(DT, which(!is.finite(DT[[.name]])), j = .name,value = replacement)))
  return(DT)
}


#Check if file format exists.
#if boolreturn == F, find the file name
#if boolreturn == T, return T.
#bamType, RNA or RFP/RPF
findFF = function(formatName, boolreturn = F, bamType = NULL){
  
  regEx = paste("\\.",formatName,"$",sep = "")
  ffList = list.files(pattern = regEx, ignore.case=TRUE)
  if(length(ffList) == 1){
    if(boolreturn)return(T)
    return(ffList)
    
  }else if((formatName == "bam") & (length(ffList) == 2) & (!is.null(bamType))){
    #find all sequences with bamtype in string  
    ffList1 = ffList[lapply(ffList,function(x) length(grep(bamType,x,value=F))) == 1]
    
    if(bamType == "RFP" & length(ffList1) == 0){
      ffList1 = ffList[lapply(ffList,function(x) length(grep("RPF",x,value=F))) == 1]
      if(length(ffList1) == 0){
        print("Error: rename the RFP.bam file to RFP.bam")
        stop()
      }
    }
    cat("using as ",bamType, "the file: ", unlist(ffList),"\n")
    if(length(ffList) > 0)
      return (ffList[1])
    else
      return(ffList1[1])
    
  } else{
    if(boolreturn)return(F)
    
    stop(cat("No ", formatName,"file in this folder or duplicates"))
  }
}

getRelativePathName = function(name){
  return (gsub(".*/", "", name))
}

rfe = function(name){#Only works on .bam right now!!!!
  return (gsub("*\\.bam", "", name))
}

loadBamFile <- function(fileLocation, typeOfBam, notPaired = T, destination = NULL){
  sortedBam <- paste0(bamFolder, getRelativePathName(fileLocation))
  if (!file.exists(p(sortedBam,".bai"))){
    if (is.null(destination)){
      sortBam(fileLocation,rfe(sortedBam)) #rfe - remove file extension
      indexBam(sortedBam)
    } else{
      sortBam(fileLocation,rfe(sortedBam)) #rfe - remove file extension
      indexBam(sortedBam)
    }
    cat("Created new ",typeOfBam,"-seq file, name:\n",sortedBam)
  }
  if (notPaired){
    if (testBAM(sortedBam)){ ##Check if this is realy necesary
      rseq <- readGAlignmentPairs(sortedBam)
    } else 
      rseq <- readGAlignments(sortedBam)
  } else{
    rseq <- readGAlignments(sortedBam)
  } 
  return(rseq)
}

#' fast import
#' 
#' versatile loader function, supports many formats
#' @param fileName the file path
fimport <- function(fileName){
  library(tools)
  if(!file.exists(fileName)) stop(paste(fileName, "does not name a file"))
  
  format <- file_ext(fileName)
  if(format == "") stop(paste("file", fileName,"have no format"))
  
  if (format == "csv") {
    return(fread(fileName))
  } else if (format == "tab") {
    return(fread(fileName, sep = "\t"))
  } else if (format == "bed") {
      return(ORFik:::cageFromFile(fileName))
  } else if (format == "bam") {
      return(readGAlignments(fileName))
  } else if (format == "rds") {
      return(readRDS(fileName))
  } else if (format == "rdata") {
    rdataname <- load(fileName)
    message(paste("file is loaded to global environment as:", rdataname))
    return(NULL)
  } else if (format == ".gtf") {
    return(makeTxDbFromGFF(fileName))
  } else {
    stop(paste("loading function for", fileName,
               "is not supported, do it manually."))
  }
}

#' Make directoy structure for orf finding
#' 
#' The main Path is ./.. relative to RCode1/ location
orfikDirs <- function(mainPath, makeDatabase = F){
  setwd(mainPath)
  print(paste("main path for project will be: ", mainPath))
  resultsLoc <- "test_results/"
  dir.create(resultsLoc)
  
  dir.create(p(resultsLoc,"New_Cage_Leaders"))
  dir.create(p(resultsLoc,"New_Cage_bedLeaders"))
  dir.create(p(resultsLoc,"rangesOfUORFs"))
  dir.create(p(resultsLoc,"bedUORFS"))
  dir.create(p(resultsLoc,"fasta"))
  
  if (makeDatabase) {
    dir.create(p(resultsLoc,"uorfIDs"))
    dir.create("dataBase")
  }
  
  print("directories created successfully")
}

#' set up cluster for pipeline
#' 
#' Cluster object saved as cl
#' @param maxCores what is max cores, if Null set to half of available cores
pipelineCluster <- function(maxCores = NULL){
  if(exists("cl") && any(class(cl) == "cluster")){
    message("cluster already exists, abort if not correct")
  } else {
    library(doParallel)
    if(is.null(maxCores)){
      maxCores = as.integer(detectCores()-(detectCores()/2)) # using half
    }
    cl <- makeCluster(maxCores)
    registerDoParallel(cl)
    assign("cl",cl,envir = .GlobalEnv)
  }
  
  message("running with number of threads: ", maxCores)
}


loadMatrix = function(mname,toGlobalEnv = T){
  return(fread(input = p(matrixFolder,mname),header = T))
}

#' save image of current session, makes name automaticly
saveRData <- function(){
  if(exists("generalName") == F){
    save.image(p(RdataFolder,"results.Rdata"))
  }else{
    save.image(paste0(RdataFolder,detailedFullName,".results.Rdata",sep=""))
  }
}

loadRData = function(rdataname,toGlobalEnv = T){
  if(toGlobalEnv)
    load(p(RdataFolder,rdataname),envir = .GlobalEnv)
  else
    load(p(RdataFolder,rdataname))
}

loadUorfID = function(rdataname){
  load(paste0(resultsFolder,"/uorfIDs/",rdataname), envir = .GlobalEnv)
}

#' From a gtf, make a new txdb with reassigned leaders from cage
#' 
#' 
makeTxDbFromGrangesLists <- function(oldGtfPath, fiveUTRs = NULL, cds = NULL, threeUTRs = NULL, tx = NULL){
  # import gtf as GRanges
  gr <- import(oldGtfPath)
  
  if (!is.null(fiveUTRs)) {
    # first find first exons of five UTRs, they have changed by cage
    firstExons <- ORFik:::firstExonPerGroup(fiveUTRs)
    posIndices <-  strandBool(firstExons)
    
    firstExons <- unlist(firstExons, use.names = FALSE)
    if (is.null(names(firstExons))) {
      firstExons <- unlist(firstExons, use.names = TRUE)
    }
    
    exonNames <- firstExons$exon_name
    
    posExonNames <- exonNames[posIndices]
    negExonNames <- exonNames[!posIndices]
    ########### For exons ##############
    # update leader starts for pos strand in gtf exons
    dt1 <- data.table(exon_names = gr2$exon_id, index = 1:length(gr))
    dt2 <- data.table(exon_names = posExonNames, indexPart = 1:length(posExonNames))
    merged <- merge(dt1, dt2, by = "exon_names", all.x = F, all.y = T)
    colnames(merged) <- c("exon_names", "index", "indexPart")
    # you will have duplicated of exon_names, but they all have same cage peak
    firstExonsPos <- firstExons[posIndices][merged$indexPart]
    
   
    start(ranges( gr[merged$index])) <-
      start(firstExonsPos)
    
    #update leaders start for neg strand in gtf exons
    dt2 <- data.table(exon_names = negExonNames, indexPart = 1:length(negExonNames))
    merged <- merge(dt1, dt2, by = "exon_names", all.x = F, all.y = T)
    colnames(merged) <- c("exon_names", "index", "indexPart")
    # you will have duplicated of exon_names, but they all have same cage peak
    firstExonsNeg <- firstExons[!posIndices][merged$indexPart]
    end(ranges( gr[merged$index])) <-
      end(firstExonsNeg)
    
    ########### For transcripts ##############
    # now update transcript starts to the biggest / lowest new start from cage
    txNames <- names(firstExons)
    
    posTxNames <- txNames[posIndices]
    negTxNames <- txNames[!posIndices]
    
    dt1 <- data.table(transcript_id = gr$transcript_id, type = gr$type, index = 1:length(gr))
    dt1 <- dt1[type == "transcript",]
    # for pos
    
    dt2 <- data.table(transcript_id = posTxNames, indexPart = 1:length(posTxNames))
    merged <- merge(dt1, dt2, by = "transcript_id", all.x = F, all.y = T)
    firstExonsPos <- firstExons[posIndices][merged$indexPart]
    start(ranges( gr[merged$index])) <-
      start(firstExonsPos)
    
    # for neg
    
    dt2 <- data.table(transcript_id = negTxNames, indexPart = 1:length(negTxNames))
    merged <- merge(dt1, dt2, by = "transcript_id", all.x = F, all.y = T)
    firstExonsNeg <- firstExons[!posIndices][merged$indexPart]
    end(ranges( gr[merged$index])) <-
      end(firstExonsNeg)
    
    ########### For genes ##############
    
    
    txdb <- makeTxDbFromGRanges(gr)
  }
}
