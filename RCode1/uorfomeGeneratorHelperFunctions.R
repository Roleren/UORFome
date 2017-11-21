source("./HelperVariables.R")
source("./NameCreator.R")
source("./GenomicGetters.R")
source("./HelperLibraries.R")
source("./CageDataIntegration.R")
source("./scanUORFs.R")
source("./Plotting&Statistics/PlotUORFome.R")
source("./HelperFunctions.R")
source("./uOrfFeatures.R")
source("./GRangesHelpers.R")


###Make the output matrix, containing normalizations, te's, lengths and names.
###Saves the matrix to inputfolder as matrix.csv
#Change ORFlengths to width ? 
makeMatrix = function(allLengths,teCDS,te5UTR,te3UTR,teUORF,transcriptNames){
  #Get the general values for the transcript, and te's for them.
  transcriptMatrix = as.data.table(cbind(allLengths,normCDSRNA,normCDSRFP,norm5UTRRNA,norm5UTRRFP,norm3UTRRNA,norm3UTRRFP,teCDS,te5UTR,te3UTR))
  #Get the uorf values
  uorfMatrix = getUOrfFeaturesMatrix(transcriptMatrix)
  
  #reduce transcriptMatrix to match uorfMatrix
  txNames1 = 1:length(transcriptMatrix$tx_name)
  names(txNames1) = transcriptMatrix$tx_name
  txNames1 = txNames1[uorfMatrix$transcriptNames]
  transcriptMatrix = transcriptMatrix[txNames1]
  
  matrix = as.data.table(cbind(transcriptMatrix,uorfMatrix))
  
  #set col classes
  class(matrix$width) = "integer"
  class(matrix$frame) = "integer"
  class(matrix$rank) = "integer"
  class(matrix$start) = "integer"
  class(matrix$end) = "integer"
  class(matrix$normUORFRNA) = "double"
  class(matrix$normUORFRFP) = "double"
  class(matrix$uorfName) = "character"
  class(matrix$overlap) = "logical"
  class(matrix$pass_filter) = "logical"
  class(matrix$UCdist) = "integer"
  class(matrix$teUORF) = "double"
  class(matrix$uorfID) = "character"
  class(matrix$tissue) = "character"
  class(matrix$seqnames) = "character"
  
  assign("matrix",matrix,envir = .GlobalEnv)
  if(exists("generalName") == F){
    write.csv(matrix, file = paste0(matrixFolder,"matrix.csv")) #filtered matrix
  }else{
    write.csv(matrix, file = paste0(matrixFolder,detailedFullName,".matrix.csv",sep="")) #filtered matrix
  }
}
#Decide name to use for rdata file and save it
saveRData = function(){
  if(exists("generalName") == F){
    save.image(p(RdataFolder,"results.Rdata"))
  }else{
    save.image(paste0(RdataFolder,detailedFullName,".results.Rdata",sep=""))
  }
}
##Get TE for leader, cds and 3' + all lengths of the different, then save them
getGeneralTEValues = function(usingNewCage,leaderBed){
  if(exists("te3UTR") == F){
    cat("finding all lengths\n")
    allLengths = getAllTranscriptLengths()
    cat("finding te's of RNA and UTRs\n")
    
    if(usingNewCage || !is.null(leaderBed)){
      print("redefining five-utr lengths for TE values, because of cage")
      new5Length = findCageUTRFivelen(fiveUTRs, allLengths$tx_name)
      allLengths$utr5_len = new5Length
      allLengths$tx_len = allLengths$utr5_len + allLengths$cds_len + allLengths$utr3_len
    }
    
    teCDS = getTE(cds,rna,RFP,allLengths$tx_len,allLengths$cds_len,allLengths,"teCDS")
    te5UTR = getTE(fiveUTRs,rna,RFP,allLengths$tx_len,allLengths$utr5_len,allLengths,"te5UTR")
    te3UTR = getTE(threeUTRs,rna,RFP,allLengths$tx_len,allLengths$utr3_len,allLengths,"te3UTR")
    
    
    assign("teCDS",teCDS,envir = .GlobalEnv)
    assign("te5UTR",te5UTR,envir = .GlobalEnv)
    assign("te3UTR",te3UTR,envir = .GlobalEnv)
    
  }
}
#check if uorf ranges already exists, if not make load them or make them from scratch
decideHowToGetUORFRanges = function(assignUorf = F,givenCage = NULL){ ###Watch out for assign uorf, might be buggy
  cat("started finding UORFS\n")
  if(UorfRangesNotExists(assignUorf,givenCage))
    rangesOfuORFs = scanUORFs(fiveUTRs,saveToFile = T, assignUorf = assignUorf)
}
#Check if uorfRanges exist already, or if must be created.
###########Should make this more failsafe!!!!!!!!!! add possibility to give ranges!!!!!!!!
UorfRangesNotExists = function(assignUorf = F, givenCage = NULL){
  if(exists("rangesOfuORFs") == F){
    if(!assignUorf){ #if not loading or assigning to global, we need cage name
      thisCage = givenCage
      assign("thisCage",thisCage,envir = .GlobalEnv)
    }
    if(file.exists(getUORFRDataName(givenCage))){#!!!Will not work for single run now!!!
      if(assignUorf){
        cat("loading rangesOfuorf from folder\n",uorfFolder,"\n")
        load(getUORFRDataName(givenCage),envir = .GlobalEnv)
      }
      return(F)
      #TODO add possibility to use bed files instead of just rdata files for uorfs
      #     }else if(file.exists(paste0(uorfBedFolder,generalName,".bed" ))){
      #       cat("loading rangesOfuorf from folder\n",uorfBedFolder)
      #       rangesOfuORFs1 = import.bed(paste0(uorfBedFolder,generalName,".bed" ))
      #       rangesOfuORFs = toGR(rangesOfuORFs)        
      #       assign("rangesOfuORFs",rangesOfuORFs,envir = .GlobalEnv)
      #       return(F)
    }else{return(T)}
  }
  return(F)
}


#####################################GETTESS##########################################
#seq: sequences of either: 5', 3', cds or uorfs
#dm1: detection method 1, RNA
#dm2: detection method 2, RFP
#seqLength: length of seq
#al: length of all sequences, might be NULL
#specificTE: what type you  want to get TE of, CDS, LEADER, 3', uorfs etc.
getTE = function(seq, dm1 = rna,dm2 = RFP,seqLength1,seqLength2,al = NULL,specificTE){
  
  libraryRna = length(dm1)
  libraryRPF = length(dm2)
  ###TE FOR UORF
  if(specificTE == "teUORF"){
    
    gr = unlist(seq, use.names = F)
    grl = GroupGRangesByOther(gr,gr$names)
    overlapRNA = countOverlaps(grl,dm1)
    overlapRFP = countOverlaps(grl,dm2)
    
    #Find lengths used for uorf fpkm values
    ORFLengths = widthPerGRangesGroup(grl,F)
    
    #find all tx-lengths used by orfs
    vector = allLengths$tx_len
    names(vector) = allLengths$tx_name
    index = match(gsub("_[0-9]*","", names(grl)),names(vector))
    vector = vector[index]
    
    norm1 = FPKMNorlization(overlapRNA,vector,libraryRna) #normalize by tx
    norm2 = FPKMNorlization(overlapRFP,ORFLengths,libraryRPF)# normalize by orf
    
    assign("normUORFRNA",norm1, envir = .GlobalEnv)
    assign("normUORFRFP",norm2, envir = .GlobalEnv)
    
    teUORF = as.matrix(norm2/norm1)
    assign("teUORF",teUORF, envir = .GlobalEnv)
    
    print("making matrix")
    makeMatrix(allLengths,teCDS,te5UTR,te3UTR,teUORF)
  }
  ###TE FOR CDS AND UTR
  else{
    overlap1 = countOverlaps(seq,dm1)
    overlap2 = countOverlaps(seq,dm2)
    
    
    norm1 = FPKMNorlization(cp(overlap1,al),seqLength1,libraryRna)
    norm2 = FPKMNorlization(cp(overlap2,al),unlist(seqLength2),libraryRPF)
    if(specificTE == "teCDS"){
      assign("normCDSRNA",norm1, envir = .GlobalEnv)
      assign("normCDSRFP",norm2, envir = .GlobalEnv)
    }else if(specificTE == "te5UTR"){
      assign("norm5UTRRNA",norm1, envir = .GlobalEnv)
      assign("norm5UTRRFP",norm2, envir = .GlobalEnv)
    }
    else{
      assign("norm3UTRRNA",norm1, envir = .GlobalEnv)
      assign("norm3UTRRFP",norm2, envir = .GlobalEnv)
    }
  }
  
  te = as.matrix(norm2/norm1)
  return(te)
}

###Script starting point
##Either run with reference to input folder, or use standard folder UORFome/test_data on kjempetuja
startUORFomeGenerator = function(arcs){
  
  if(lArcs == 4) #Run from kjempetuja@uib hakontj account or other people
    getMatrix(usingNewCage = as.logical(arcs[1]),cageName = arcs[2],rnaSeq = arcs[3],rfpSeq = arcs[4])
  if(lArcs == 5) #Run from kjempetuja@uib hakontj account or other people
    getMatrix(usingNewCage = as.logical(arcs[1]),cageName = arcs[2],rnaSeq = arcs[3],rfpSeq = arcs[4],tissueUsed = as.character(arcs[5]))
  if(lArcs == 1)
    getMatrix()
  print("script finished")
}
