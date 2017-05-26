source("/export/valenfs/projects/uORFome/RCode1/HelperVariables.R")
##Get riboseq file and read it
getRFP = function(rfpSeq){
  if(!is.null(rfpSeq)){ #remember to fix this when bed files arrive!!!!
    if(testBAM(rfpSeq)){
      RFP = readGAlignmentPairs(rfpSeq)
    }else 
      RFP = readGAlignments(rfpSeq)
  }
  else if(findFF("bed",boolreturn = T)){
    RFP = import.bed(findFF("bed"))
  }else{ if(testBAM(findFF("bam",bamType = "RFP"))){
    RFP = readGAlignmentPairs(findFF("bam",bamType = "RFP"))
  }else
    RFP = readGAlignments(findFF("bam",bamType = "RFP"))
  }
  assign("RFP",RFP,envir = .GlobalEnv)
}
##Get rna seq file and read it
getRNAseq = function(rnaSeq){
  if(!is.null(rnaSeq)){ #remember to fix this when bed files arrive!!!!
    if(testBAM(rnaSeq)){
      rna = readGAlignmentPairs(rnaSeq)
    }else 
      rna = readGAlignments(rnaSeq)
  }
  else if(testBAM(findFF("bam",bamType = "RNA"))){
    rna = readGAlignmentPairs(findFF("bam",bamType = "RNA"))
  }else 
    rna = readGAlignments(findFF("bam",bamType = "RNA"))
  assign("rna",rna,envir = .GlobalEnv)
}

getFasta = function(){
  if(exists("fasta") == F){
    fasta =  readDNAStringSet(fastaName) ##Get fasta file
    assign("fasta",fasta,envir = .GlobalEnv)
  }
  #faiFileNotexists = !findFF("fai",T)## Get fasta indexed file
  #if(faiFileNotexists)indexFa(fastaName)
  #if(faiFileNotexists){indexFa(findFF("fa"))}
  if(exists("fa") == F){
    fa = FaFile(faiName)
    assign("fa",fa,envir = .GlobalEnv)
  }
}

getGTF = function(){
  if(exists("Gtf") == F){
    #load(file = "/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/Gtf.rdata")
    Gtf = makeTxDbFromGFF(gtfName)
    assign("Gtf",Gtf,envir = .GlobalEnv)
  }
}

getLeaders = function(leaderBed = NULL,usingNewCage = F, cageName = NULL){
  if(!is.null(leaderBed)){
    cat("retrieving 5utrs from bed file: ", leaderBed)
    fiveUTRs = import.bed(leaderBed)
    fiveUTRs = as.data.frame(fiveUTRs)
    fiveUTRstest <- split(fiveUTRs, fiveUTRs$name)
    fiveUTRstest = lapply(fiveUTRstest,function(x) as(x,"GRanges"))
    
    fiveUTRstest = GRangesList(fiveUTRstest)
    fiveUTRs = fiveUTRstest
  }
  else{
    getGTF()
    
    if(exists("fiveUTRs") == F){
      fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
      if(usingNewCage){
        cat("Using cage..")
        if(is.null(cageName)){
          cat("error no cageName, continueing without cage")
          assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
          return
        }
        fiveUTRs = getNewfivePrimeUTRs(fiveUTRs,cageName)
        cat("finished new 5' UTRs")
      }
    }
  }
  assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
}