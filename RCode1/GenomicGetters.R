
##Get riboseq file and read it
getRFP = function(rfpSeq){
  
  if(!is.null(rfpSeq) && exists("RFP") == F){ #remember to fix this when bed files arrive!!!!
    RFP = loadBamFile(rfpSeq,"rfp")
    assign("RFP",RFP,envir = .GlobalEnv)
  }
}
##Get rna seq file and read it
getRNAseq = function(rnaSeq){
  
  if(!is.null(rnaSeq) && exists("rna") == F){ #remember to fix this when bed files arrive!!!!
    rna = loadBamFile(rnaSeq,"rna")
    assign("rna",rna,envir = .GlobalEnv)
  }
}

#Get the fasta file, fasta indexed file and assign them
getFasta = function(){

  if(exists("fa") == F){ #index files
#       if(exists("fasta") == F){ #make possible to create fai from fa
#         print("loading fasta files")
#         fasta =  readDNAStringSet(fastaName) ##Get fasta file
#         assign("fasta",fasta,envir = .GlobalEnv)
#       }
    fa = FaFile(faiName)
    assign("fa",fa,envir = .GlobalEnv)
  }
}

#get sequences from a GRangeslist
getSequencesFromFasta = function(grl){
  getFasta() #get fasta and fai
  seqs = extractTranscriptSeqs(fa, transcripts = grl)
  assign("seqs",seqs,envir = .GlobalEnv)
}
#Get the Genomic transcript format, currently using GRch38 data
getGTF = function(){
  if(exists("Gtf") == F){
    print("loading human GTF GRch38")
    library(AnnotationDbi)
    Gtf = loadDb(gtfdb)
    #Gtf = makeTxDbFromGFF(gtfName) #Fix this!!!!!!!!!!!!
    assign("Gtf",Gtf,envir = .GlobalEnv)
  }
}
#Get the coding sequences from the gtf file
getCDS = function(assignIt = T){
  if(exists("cds",mode = "S4") == F){
    cds = cdsBy(Gtf,"tx",use.names = T)
    if(assignIt){
      assign("cds",cds,envir = .GlobalEnv)
      return(cds)
    }else{
      return(cds)
    }
  }
}
#Get the 3' sequences from the gtf file
getThreeUTRs = function(){
  if(!exists("threeUTRs")){
    threeUTRs = threeUTRsByTranscript(Gtf,use.names = T)
    assign("threeUTRs",threeUTRs,envir = .GlobalEnv)
  }
}

#Get the 5' leaders, either from gtf, cage data to reassign the transcription start site(TSS), or load from existing data either as .rdata or .bed (bed6)
getLeaders = function(leaderBed = NULL,usingNewCage = F, cageName = NULL,leader = NULL,assignLeader = T){
  
  if(!is.null(cageName)){#check if leader is already made, either as rdata or bed
    if(file.exists(paste0(leadersFolder,getRelativePathName(cageName),".leader.rdata"))){
      leader = paste0(leadersFolder,getRelativePathName(cageName),".leader.rdata")
    }else if(file.exists(paste0(leadersbedFolder,getRelativePathName(cageName),".leader.bed"))){
      leaderBed = paste0(leadersbedFolder,getRelativePathName(cageName),".leader.bed")
    }
  }
  
  if(!is.null(leader)){ #load as rdata
    print("loading leader from pre-existing rdata")
    if(assignLeader)
      load(leader)
  }
  else if(!is.null(leaderBed)){
    cat("retrieving 5utrs from bed file: ", leaderBed,"\n")
    fu = import.bed(leaderBed)
    vec = vector(mode="list",length = length(unique(fu$name)))
    names(vec) <- unique(fu$name)
    for(i in unique(fu$name)){ #combine by name, to make transcripts by exons
      vec[[i]] <- fu[fu$name == i]
    }
    fiveUTRs = GRangesList(vec)
    
    exportNamerdata = paste0(leadersFolder,getRelativePathName(p(cageName,".leader.rdata")))
    save(fiveUTRs,file = exportNamerdata)
  }
  else{ #create from scratch
    
    if(exists("fiveUTRs") == F){
      cat("creating leader from scratch\n")
      getGTF()
      cat("loading Leader from gtf\n")
      fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
      if(usingNewCage){
        print("Using cage.. ")
        if(is.null(cageName)){
          print("error no cageName, continueing without cage")
          assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
          return
        }
        getCDS()
        fiveUTRs = reassignTSSbyCage(fiveUTRs,cageName, cds = cds)
        if(1){ #TODO add possibility to not save utrs, now it always saves
          exportNamebed = paste0(leadersbedFolder,getRelativePathName(p(cageName,".leader.bed")))
          exportNamerdata = paste0(leadersFolder,getRelativePathName(p(cageName,".leader.rdata")))
          print("exporting new leaders")
          export.bed(unlist(fiveUTRs),exportNamebed)
          save(fiveUTRs,file = exportNamerdata)
        }
        print("finished new 5' UTRs")
      }
    }
    else{
      print("fiveUTRs already exists! cancel if this is wrong!")
    }
  }
  if(assignLeader)
    assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
  
  print("finished loading leaders")
}

#Get the upstream open reading frames from the 5' leader sequences, given as GRangesList
getUnfilteredUORFs = function(fiveUTRs, assignRanges = T){
  getSequencesFromFasta(fiveUTRs)
  rangesOfuORFs = find_in_frame_ORFs(grl = fiveUTRs,fastaSeqs = seqs,minimumLength = 2)
  
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