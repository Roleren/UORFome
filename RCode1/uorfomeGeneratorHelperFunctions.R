source("/export/valenfs/projects/uORFome/RCode1/HelperLibraries.R")
source("/export/valenfs/projects/uORFome/RCode1/GenomicGetters.R")
source("/export/valenfs/projects/uORFome/RCode1/CageDataIntegration.R")
source("/export/valenfs/projects/uORFome/RCode1/scanUORFs.R")
source("/export/valenfs/projects/uORFome/RCode1/Plotting&Statistics/PlotUORFome.R")
source("/export/valenfs/projects/uORFome/RCode1/HelperFunctions.R")



### Print info about specific run
infoPrints = function(data,doubleBAM,usingNewCage,cageName,leaderBed,rnaSeq = NULL,rfpSeq = NULL){
  print("Estimated normal time is 5 hours\n")
  cat("folder used is:",data,"\n")
  if(doubleBAM == T)
    cat("using bam-file for RFP\n")
  if(usingNewCage)
    cat("using cage-file named: \n",cageName,"\n")
  setwd(data)
  makeGeneralName(cageName,leaderBed,rnaSeq,rfpSeq)
  cat("starting loading objects\n")
}

makeGeneralName = function(cageName = NULL,leaderBed = NULL,rnaSeq = NULL,rfpSeq = NULL){
  if(!is.null(cageName)){ #General name
    generalName = gsub(".*/", "", cageName)
    generalName = strsplit(generalName,".hg38*.")[[1]][1]
  }else if(!is.null(leaderBed)){
    generalName = gsub(".*/", "", leaderBed)
    generalName = strsplit(generalName,"*.Leader.bed")[[1]][1]
  }else{
    generalName = paste0("UORF run: ",Sys.time())
  }
  cat("Name of run is: ",generalName,"\n")
  assign("generalName",generalName,envir = .GlobalEnv)
  
  detailedFullName = ""
  if(!is.null(cageName) && !is.null(rnaSeq) && !is.null(rfpSeq)){
    detailedFullName = paste0(generalName,";",rnaSeq,";",rfpSeq)
  }
  assign("detailedFullName",detailedFullName,envir = .GlobalEnv)
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
    teCDS = getTE(cds,rna,RFP,allLengths$tx_len,allLengths$cds_len,allLengths,"teCDS")
    
    if(usingNewCage || !is.null(leaderBed)){
      testbest = findCageUTRFivelen(fiveUTRs, allLengths$tx_name)
      allLengths$utr5_len = testbest
      allLengths$tx_len = allLengths$utr5_len + allLengths$cds_len + allLengths$utr3_len
    }
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

startUORFomeGenerator = function(arcs){
  
    if(length(arcs) == 1)
      getMatrix(data = arcs[1])
    
    else if(length(arcs) == 2)
      getMatrix(arcs[1],  doubleBAM = arcs[2])
    
    else if(length(arcs) == 3)
      getMatrix(arcs[1],  doubleBAM = arcs[2], usingCage = arcs[3])
    
    else if(length(arcs) == 4)
      getMatrix(arcs[1],  doubleBAM = arcs[2], usingCage = arcs[3],cageName = arcs[4])
    
  
  print("script finished")
}
