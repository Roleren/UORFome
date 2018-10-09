
#' Get the Genomic transcript format, currently using GRch38 data
getGTF = function(assignIt = T){
  if(exists("Gtf") == F){
    print("loading human GTF GRch38")
    library(AnnotationDbi)
    if(is.null(gtfdb)){
      if(file.exists(gtfName)){
        Gtf = makeTxDbFromGFF(gtfName)
      } else {
        stop("could not find gtf file, check path")
      }
    }
    Gtf = loadDb(gtfdb)
    if(assignIt)
      assign("Gtf",Gtf,envir = .GlobalEnv)
  } else if(!dbIsValid(Gtf$conn)) {
    Gtf = loadDb(gtfdb)
    assign("Gtf",Gtf,envir = .GlobalEnv)
  }
}

#' Get transcripts from gtf
getTx <- function(assignIt = F){
  if (exists("tx",mode = "S4") == F) {
    getGTF()
    tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)
    if (assignIt) {
      assign("tx",tx, envir = .GlobalEnv)
      return(tx)
    } else {
      return(tx)
    }
  } 
}

#Get the coding sequences from the gtf file
getCDS = function(assignIt = T){
  if (exists("cds",mode = "S4") == F) {
    getGTF()
    cds = cdsBy(Gtf,"tx", use.names = TRUE)
    if (assignIt) {
      assign("cds", cds, envir = .GlobalEnv)
      return(cds)
    } else {
      return(cds)
    }
  }
}

#Get the 3' sequences from the gtf file
getThreeUTRs = function(){
  if (!exists("threeUTRs")) {
    if (!exists(p(dataFolder, "/threeUTRs.rdata"))) {
      getGTF()
      threeUTRs = threeUTRsByTranscript(Gtf, use.names = TRUE)
      assign("threeUTRs", threeUTRs, envir = .GlobalEnv)
    } else {
      load(p(dataFolder, "/threeUTRs.rdata"), envir = .GlobalEnv)
    }
  }
}

#' Get the 5' leaders, either from gtf, cage data to reassign
#' the transcription start site(TSS), or load from existing data
#' either as .rdata or .bed (bed6)
getLeaders = function(cageName = NULL, assignLeader = T, exportUorfRegions = T){
  if(exists("fiveUTRs") == F){
    cat("creating leader from scratch\n")
    
    if (file.exists(p(dataFolder,"/leader.rdata"))) {
      load(p(dataFolder,"/leader.rdata"))
    } else {
      getGTF()
      cat("loading Leader from gtf\n")
      fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
      save(fiveUTRs, file = p(dataFolder,"/leader.rdata"))
    }
    
    if (!is.null(cageName)) {
      print("Using cage.. ")
      
      fiveUTRs = ORFik::reassignTSSbyCage(fiveUTRs,cageName)
      if (exportUorfRegions) {
        getCDS()
        uORFSeachRegion <- ORFik:::addCdsOnLeaderEnds(fiveUTRs, cds, onlyFirstExon = F)
        uORFSeachRegion <- sortPerGroup(uORFSeachRegion)
        print("exporting new uorf regions")
        exportNamerdata = paste0(regionUORFs,
                                 getRelativePathName(p(cageName,
                                                       ".regionUORF.rdata")))
        save(uORFSeachRegion, file = exportNamerdata)
      }
      
      print("exporting new leaders")
      exportNamerdata = paste0(leadersFolder,
                               getRelativePathName(p(cageName,
                                                     ".leader.rdata")))
      save(fiveUTRs,file = exportNamerdata)
      print("finished new 5' UTRs")
    }
  }
  else{
    print("fiveUTRs already exists! cancel if this is wrong!")
  }
  
  if(assignLeader)
    assign("fiveUTRs",fiveUTRs, envir = .GlobalEnv)
  
  print("finished loading leaders")
}

#' Get the leader that should span all uorfs
leaderAllSpanning <- function(){
  getGTF()
  getLeaders()
  getCDS()
  
  fiveUTRs <- ORFik:::addFirstCdsOnLeaderEnds(
    ORFik:::assignFirstExons(ORFik:::extendsTSSexons(fiveUTRs),
                             fiveUTRs), cds)
  fiveUTRs <- sortPerGroup(fiveUTRs)
  return(fiveUTRs)
}

leaderCage <- function(width.cds = TRUE){
  if(width.cds) {
    load(p(dataBaseFolder,"/CageFiveUTRsWithCDS.rdata"))
    return(CageFiveWithCDS)
  }
  load(p(dataBaseFolder,"/CageFiveUTRs.rdata"))
  
  return(CageFiveUTRs)
}


#' Get the fasta indexed file
#' 
#' if assignIt is TRUE, the object is not return to local scope
#' Only assigned to globalenvir
getFasta = function(filePath = NULL, assignIt = T){
  
  if(exists("fa") == F){ #index files
    if (is.null(filePath)){
      fa = FaFile(faiName)
    } else {
      fa = FaFile(filePath)
    }
    if (assignIt){
      assign("fa",fa,envir = .GlobalEnv)
    } else {
      return(fa)
    }
  }
}

#' get sequences from a GRangeslist
getSequencesFromFasta = function(grl, isSorted = F){
  getFasta() #get .fai
  if(!isSorted) grl <- ORFik:::sortPerGroup(grl)
  seqs = extractTranscriptSeqs(fa, transcripts = grl)
  assign("seqs",seqs,envir = .GlobalEnv)
}

getAll <- function(include.cage = T){
  getFasta()
  
  getCDS()
  getThreeUTRs()
  getLeaders()
  tx <- getTx()
  
  
  #or with extension
  if (include.cage) {
    cageFiveUTRs <- leaderCage()
    assign("cageFiveUTRs", cageFiveUTRs,  envir = .GlobalEnv)
    getCageTx()
  }
  return(NULL)
}

# delete all
da <- function(){
  if (exists("threeUTRs")) {
    rm(threeUTRs, envir = .GlobalEnv)
  }
  if (exists("cds", mode = "S4")){
    rm(cds, envir = .GlobalEnv)
  }
  if (exists("tx", mode = "S4")) {
    rm(tx, envir = .GlobalEnv)
  }
  if (exists("fiveUTRs", mode = "S4")) {
    rm(fiveUTRs, envir = .GlobalEnv)
  }
}

#' Convenience wrapper from findMapORFs
#' Get the upstream open reading frames from the 5' leader sequences, given as GRangesList 
getUnfilteredUORFs = function(uORFSeachRegion, assignRanges = T, isSorted = F,
                              startCodons = "ATG", groupByTx = F){
  
  getSequencesFromFasta(uORFSeachRegion, isSorted)
  
  
  rangesOfuORFs = ORFik::findMapORFs(grl = uORFSeachRegion,seqs = seqs, longestORF = F,
                                     minimumLength = 2, startCodon = startCodons,
                                     groupByTx = groupByTx)
  
  if(assignRanges)
    assign("rangesOfuORFs",rangesOfuORFs,envir = .GlobalEnv)
  print("finished unfiltered UORFs")
  return(rangesOfuORFs)
}

  
getAllTranscriptLengths = function(){
  load(p(dataFolder,"/transcriptLengths.rdata"))
  #transcriptLengths(Gtf,with.cds_len = T,with.utr5_len = T,with.utr3_len = T)
  assign("allLengths",allLengths,envir = .GlobalEnv)
}

#Get the number of overlaps between the upstream open reading frames and the coding sequences
getUOrfOverlaps = function(){
  overlap1 = findOverlaps(cds,rangesOfuORFs)
  overlapCount = countOverlaps(cds,rangesOfuORFs)
  numberOfOverlaps = sum(overlapCount >= 1)
  overlapHitsIndex = overlapCount[overlapCount == 1]
  return(numberOfOverlaps)
}

getCageTx <- function() {
  if (file.exists(p(dataBaseFolder, "/cageTx.rdata"))) {
    load(p(dataBaseFolder, "/cageTx.rdata"), envir = .GlobalEnv)
  } else {
    getTx()
    cageFiveUTRs <- leaderCage()
    tx <- ORFik:::extendLeaders(tx, cageFiveUTRs)
    assign("tx", tx,  envir = .GlobalEnv)
    save(tx, file = p(dataBaseFolder, "/cageTx.rdata"))
  }
  return(NULL)
}

#' get ribo seq and gtf parts
getRiboTest <- function(onlyGRL = FALSE, rfpIndex = 7) {
  if (!onlyGRL) {
    getAll()  
    assign(x = "RNA", NULL, envir = .GlobalEnv)
    #or with extension
    fiveUTRs <- leaderCage(F)
    tx <- extendLeaders(tx, fiveUTRs)
    assign("fiveUTRs", fiveUTRs, envir = .GlobalEnv)
    assign("tx", tx, envir = .GlobalEnv)
  }
  rfpList <- grep(pattern = "merged",x = list.files(rfpFolder), value = T)
  RFPPath <- p(rfpFolder, rfpList[rfpIndex])
  assign(x = "grl", getUorfsInDb(), envir = .GlobalEnv)
  assign(x = "RFP", fread.bed(RFPPath), envir = .GlobalEnv)
  
}