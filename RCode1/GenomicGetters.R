
#' Get riboseq file and read it
getRFP = function(rfpSeq){
  library(tools)
  if(!is.null(rfpSeq) && exists("RFP") == F){
    if (file_ext(rfpSeq) == "bam") {
      RFP <- loadBamFile(rfpSeq,"rfp")
      assign("RFP",RFP,envir = .GlobalEnv)
    }
    if (file_ext(rfpSeq) == "bed") {
      RFP <- ORFik:::cageFromFile(rfpSeq)
      assign("RFP",RFP,envir = .GlobalEnv)
    }
  }
}

#' Get rna seq file and read it
getRNAseq = function(rnaSeq){
  
  if(!is.null(rnaSeq) && exists("rna") == F){ #remember to fix this when bed files arrive!!!!
    rna = loadBamFile(rnaSeq,"rna")
    assign("rna",rna,envir = .GlobalEnv)
  }
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
  }
}

#' Get transcripts from gtf
getTx <- function(assignIt = F){
  if (exists("tx",mode = "S4") == F) {
    if (exists("Gtf") == F) {
      getGTF()
    }
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
    if (exists("Gtf") == F) {
      getGTF()
    }
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
    if (exists("Gtf") == F) {
      getGTF()
    }
    threeUTRs = threeUTRsByTranscript(Gtf, use.names = TRUE)
    assign("threeUTRs", threeUTRs, envir = .GlobalEnv)
  }
}

#' Get the 5' leaders, either from gtf, cage data to reassign
#' the transcription start site(TSS), or load from existing data
#' either as .rdata or .bed (bed6)
getLeaders = function(leaderBed = NULL, usingNewCage = F, cageName = NULL,
                      leader = NULL, assignLeader = T, exportAsBed = F){
  
  if(!is.null(cageName)){#check if leader is already made, either as rdata or bed
    if(file.exists(paste0(leadersFolder,getRelativePathName(cageName),".leader.rdata"))){
      leader = paste0(leadersFolder,getRelativePathName(cageName),".leader.rdata")
    }else if(file.exists(paste0(leadersbedFolder,getRelativePathName(cageName),".leader.bed"))){
      leaderBed = paste0(leadersbedFolder,getRelativePathName(cageName),".leader.bed")
    }
  }
  
  if(!is.null(leader)){ #load as rdata
    print("loading leader from pre-existing rdata")
    load(leader)
    if(!assignLeader)
      return(fiveUTRs)
  }
  else if(!is.null(leaderBed)){
    cat("retrieving 5utrs from bed file: ", leaderBed,"\n")
    fu = import.bed(leaderBed)
    vec = vector(mode="list",length = length(unique(fu$name)))
    names(vec) <- unique(fu$name)
    #combine by name, to make transcripts by exons
    for(i in unique(fu$name)){ 
      vec[[i]] <- fu[fu$name == i]
    }
    fiveUTRs = GRangesList(vec)
    
    exportNamerdata = paste0(leadersFolder,getRelativePathName(
      p(cageName,".leader.rdata")))
    save(fiveUTRs,file = exportNamerdata)
  }
  else{ #create from scratch
    
    if(exists("fiveUTRs") == F){
      cat("creating leader from scratch\n")
      getGTF()
      
      if (file.exists(p(dataFolder,"/leader.rdata"))) {
        load(p(dataFolder,"/leader.rdata"))
      } else {
        cat("loading Leader from gtf\n")
        fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
        save(fiveUTRs, p(dataFolder,"/leader.rdata"))
      }
      
      if (usingNewCage) {
        print("Using cage.. ")
        if (is.null(cageName)) {
          print("error no cageName, continueing without cage")
          assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
          return
        }
        getCDS()
        fiveUTRs = ORFik::reassignTSSbyCage(fiveUTRs,cageName, cds = cds)
        fiveUTRs <- ORFik:::sortPerGroup(fiveUTRs)
        print("exporting new leaders")
        exportNamerdata = paste0(leadersFolder,
                                 getRelativePathName(p(cageName,
                                                       ".leader.rdata")))
        save(fiveUTRs,file = exportNamerdata)
        if (exportAsBed) {
          #TODO add possibility to not save utrs, now it always saves
          exportNamebed = paste0(leadersbedFolder,
                                 getRelativePathName(p(cageName,
                                                       ".leader.bed")))
          export.bed(unlist(fiveUTRs),exportNamebed)
        }
        print("finished new 5' UTRs")
      }
    }
    else{
      print("fiveUTRs already exists! cancel if this is wrong!")
    }
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

getAll <- function(include.cage = T, extendTx = F){
  getFasta()
  getCDS()
  getThreeUTRs()
  getLeaders()
  tx <- getTx()
  
  
  #or with extension
  if (include.cage) {
    cageFiveUTRs <- leaderCage()
    assign("cageFiveUTRs", cageFiveUTRs,  envir = .GlobalEnv)
  }
  if (extendTx){
    tx <- ORFik:::extendLeaders(tx, cageFiveUTRs)
  }
  assign("tx", tx,  envir = .GlobalEnv)
  
  return(NULL)
}

# delete all
da <- function(){
  if (exists("threeUTRs")) {
    rm(threeUTRs)
  }
  if (exists("cds", mode = "S4")){
    rm(cds)
  }
  if (exists("tx", mode = "S4")) {
    rm(tx)
  }
  if (exists("fiveUTRs", mode = "S4")) {
    rm(fiveUTRs)
  }
}

#' Convenience wrapper from findMapORFs
#' Get the upstream open reading frames from the 5' leader sequences, given as GRangesList 
getUnfilteredUORFs = function(fiveUTRs, assignRanges = T, isSorted = F,
                              startCodons = "ATG"){
  
  getSequencesFromFasta(fiveUTRs, isSorted)
  
  
  rangesOfuORFs = ORFik::findMapORFs(grl = fiveUTRs,seqs = seqs,minimumLength = 2, startCodon = startCodons)
  
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