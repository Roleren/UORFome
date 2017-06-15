source("./GenomicGetters.R")
source("./HelperLibraries.R")
source("./CageDataIntegration.R")
source("./scanUORFs.R")
source("./Plotting&Statistics/PlotUORFome.R")
source("./HelperFunctions.R")



### Print info about specific run
infoPrints = function(doubleBAM,usingNewCage,cageName,leaderBed,rnaSeq = NULL,rfpSeq = NULL){
  print("Estimated normal time is 5 hours\n")
  cat("folder used is:",normalizePath(dataFolder),"\n")
  if(doubleBAM == T)
    cat("using bam-file for RFP\n")
  if(usingNewCage)
    cat("using cage-file named: \n",cageName,"\n")
  
  makeGeneralName(cageName,leaderBed,rnaSeq,rfpSeq)
  cat("starting loading objects\n")
}

makeGeneralName = function(cageName = NULL,leaderBed = NULL,rnaSeq = NULL,rfpSeq = NULL){
  if(!is.null(cageName)){ #General name
    thisCage = getRelativePathName(cageName)
    generalName = strsplit(thisCage,".hg38*.")[[1]][1]
  }else if(!is.null(leaderBed)){
    generalName = getRelativePathName(leaderBed)
    generalName = strsplit(generalName,"*.Leader.bed")[[1]][1]
  }else{
    generalName = p("UORF run: ",Sys.time())
  }
  cat("generalName of run is: ",generalName,"\n")
  assign("generalName",generalName,envir = .GlobalEnv)
  
  
  if(!is.null(cageName) && !is.null(rnaSeq) && !is.null(rfpSeq)){
    detailedFullName = paste0(generalName,"_",getRelativePathName(rnaSeq),"_",getRelativePathName(rfpSeq))
    cat("detailed Full Name of run is: ",detailedFullName,"\n")
  }
}

###Make the output matrix, containing normalizations, te's, lengths and names.
###Saves the matrix to inputfolder as matrix.csv
makeMatrix = function(allLengths,teCDS,te5UTR,te3UTR,teUORF,transcriptNames){
  matrixA = cbind(allLengths,normCDSRNA,normCDSRFP,norm5UTRRNA,norm5UTRRFP,norm3UTRRNA,norm3UTRRFP,teCDS,te5UTR,te3UTR)
  
  matrixB = cbind(transcriptNames,agUORFRibo$Group.1, teUORF, normUORFRNA, normUORFRFP,ORFLengths)
  colnames(matrixB)[2] <- "namesOFuorfs"
  colnames(matrixB)[3] <- "teUORF"
  matrix = merge(matrixA, matrixB, by.x = "tx_name", by.y = "transcriptNames", all = T)
  write.csv(matrix, file = "matrixRaw.csv") #unfiltered matrix
  ###Index 1 remove all bad values from teCDS
  index = matrix[,"teCDS"] == 0 | is.na(matrix[,"teCDS"] ) | is.infinite(matrix[,"teCDS"] ) | is.nan(matrix[,"teCDS"] )
  matrix = matrix[!index,]
  
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
    save.image(paste0(RdataFolder,"results.Rdata"))
  }else{
    save.image(paste0(RdataFolder,detailedFullName,".results.Rdata",sep=""))
  }
}
##Get TE for leader, cds and 3' + all lengths of the different, then save them
getGeneralTEValues = function(usingNewCage,leaderBed){
  if(exists("te3UTR") == F){
    cat("finding all lengths\n")
    allLengths = transcriptLengths(Gtf,with.cds_len = T,with.utr5_len = T,with.utr3_len = T)
    cat("finding te's of RNA and UTRs\n")
    
    if(usingNewCage || !is.null(leaderBed)){
      new5Length = findCageUTRFivelen(fiveUTRs, allLengths$tx_name)
      allLengths$utr5_len = new5Length
      allLengths$tx_len = allLengths$utr5_len + allLengths$cds_len + allLengths$utr3_len
    }
    
    teCDS = getTE(cds,rna,RFP,allLengths$tx_len,allLengths$cds_len,allLengths,"teCDS")
    te5UTR = getTE(fiveUTRs,rna,RFP,allLengths$tx_len,allLengths$utr5_len,allLengths,"te5UTR")
    te3UTR = getTE(threeUTRs,rna,RFP,allLengths$tx_len,allLengths$utr3_len,allLengths,"te3UTR")
    
    assign("allLengths",allLengths,envir = .GlobalEnv)
    assign("teCDS",teCDS,envir = .GlobalEnv)
    assign("te5UTR",te5UTR,envir = .GlobalEnv)
    assign("te3UTR",te3UTR,envir = .GlobalEnv)
    
  }
}

decideHowToGetUORFRanges = function(){
  cat("started finding UORFS\n")
  if(UorfRangesNotExists())
    rangesOfuORFs = scanUORFs(fiveUTRs,saveToFile = T)
}

#####################################GETTESS##########################################
#seq: sequence, 
#dm1: detection method 1, RNA
#dm2: detection method 2, RFP
#seqLength: length of seq
#al: length of all sequences, might be NULL
#specificTE: what type you  want to get TE of, CDS, LEADER, 3' etc.
getTE = function(seq, dm1 = rna,dm2 = RFP,seqLength1,seqLength2,al = NULL,specificTE){
 
  libraryRna = length(dm1)
  libraryRPF = length(dm2)
  ###TE FOR UORF
  if(specificTE == "teUORF"){
    overlap1 = countOverlaps(unlist(seq),dm1)
    overlap2 = countOverlaps(unlist(seq),dm2)
    
    grl = as.data.frame(rangesOfuORFs)
    
    agUORFRNA = aggregate(overlap1,by = list(grl$names), sum)
    agUORFRibo = aggregate(overlap2,by = list(grl$names), sum)
    assign("agUORFRibo",agUORFRibo, envir = .GlobalEnv)
    
    ORFLengths = aggregate(grl[,"width"],by=list(grl[,"names"]),sum)
    assign("ORFLengths",ORFLengths, envir = .GlobalEnv)
    transcriptNames = gsub("_[0-9]*","", agUORFRNA$Group.1)
    assign("transcriptNames",transcriptNames, envir = .GlobalEnv)
    vector = allLengths$tx_len
    names(vector) = allLengths$tx_name
    index = match(transcriptNames,names(vector))
    vector = vector[index]
    
    norm1 = FPKMNorlization(agUORFRNA$x,vector,libraryRna)
    norm2 = FPKMNorlization(agUORFRibo$x,ORFLengths$x,libraryRPF)
    
    assign("normUORFRNA",norm1, envir = .GlobalEnv)
    assign("normUORFRFP",norm2, envir = .GlobalEnv)
    
    cat("making matrix\n")
    teUORF = norm2/norm1
    teUORF = as.matrix(teUORF)
    makeMatrix(allLengths,teCDS,te5UTR,te3UTR,teUORF,transcriptNames)
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
  te = norm2/norm1
  te = as.matrix(te)
  return(te)
}

###Script starting point
##Either run with reference to input folder, or use standard folder UORFome/test_data on kjempetuja
startUORFomeGenerator = function(arcs){
  
  if(lArcs == 4) #Run from kjempetuja@uib hakontj account or other people
    getMatrix(usingNewCage = as.logical(arcs[1]),cageName = arcs[2],rnaSeq = arcs[3],rfpSeq = arcs[4])
  if(lArcs == 1)
    getMatrix()
  print("script finished")
}
