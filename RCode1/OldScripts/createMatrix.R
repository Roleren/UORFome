arcs = commandArgs(trailingOnly = T)

####FOR FUNCTION, might need to remove nan in NORMRNA######!!

library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(Biostrings)


#3. creates normalzation of the count
FPKMNorlization = function(counts, lengthSize, librarySize){
  
  result = (as.numeric(counts)*(10^9)) / (as.numeric(lengthSize)*as.numeric(librarySize))
  return(result)
}

#Remove nan, 0 na and inf
removeNAN = function(vec,column = NULL, saveIndex = F){
  if(is.null(column)){
    index1 = vec == 0 | is.na(vec) | is.infinite(vec) | is.nan(vec)
    vec = vec[!index1]
  }
  else{
    index1 = vec[,column] == 0 | is.na(vec[,column] ) | is.infinite(vec[,column] ) | is.nan(vec[,column] )
    vec = vec[!index1,]
  }
  if(saveIndex){
    assign("index4", index1, envir = .GlobalEnv)
  }
  return(vec)
}

########################MAIN###################################
#produce matrix that contains Translation effiencies(TE's) of cds, utr and UORFs
getMatrix = function(data = "/export/valenfs/projects/uORFome/test_data", saveToGlobal = F){
  print("Estimated normal time is 15 minutes")
  setwd(data)
  
  if(exists("fiveUTRs") == F){
    cat("starting loading objects")
    rna = readGAlignmentPairs(findFF("bam"))
    Gtf = makeTxDbFromGFF(findFF("gtf"))
    RFP = import.bed(findFF("bed"))
    fasta =  readDNAStringSet(findFF("fa"))
    cds = cdsBy(Gtf,"tx",use.names = T)
    fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
    if(saveToGlobal){
      print("saving objects to global")
      assign("rna",rna,envir = .GlobalEnv)
      assign("Gtf",Gtf,envir = .GlobalEnv)
      assign("RFP",RFP,envir = .GlobalEnv)
      assign("fasta",fasta,envir = .GlobalEnv)
      assign("cds",cds,envir = .GlobalEnv)
      assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
    }
    cat("finished loading objects\n")
  }

  if(exists("teUTR") == F){
    cat("finding all lengths")
    allLengths = transcriptLengths(Gtf,with.cds_len = T,with.utr5_len = T,with.utr3_len = T)
    cat("finding te's of RNA and UTR\n")
    teRNA = getTE(cds,rna,RFP,allLengths$tx_len,allLengths$cds_len,allLengths)
    teUTR = getTE(fiveUTRs,rna,RFP,allLengths$tx_len,allLengths$utr5_len,allLengths)
    
    if(saveToGlobal){
      
      assign("allLengths",allLengths,envir = .GlobalEnv)
      assign("teRNA",teRNA,envir = .GlobalEnv)
      assign("teUTR",teUTR,envir = .GlobalEnv)
    }
  }
  cat("started finding UORFS\n")
  # teUORF = getTE(scanUORFs(fiveUTRs,fasta,saveToFile = T),rna,RFP,sapply(rangesOfuORFs, function(x) sum(width(x))))
  if(exists("rangesOfuORFs") == F){
    rangesOfuORFs = scanUORFs(fiveUTRs,fasta,saveToFile = T)
  }
  teUORF = getTE(rangesOfuORFs,rna,RFP,allLengths$tx_len,sapply(rangesOfuORFs, function(x) sum(width(x))))
  cat("making matrix\n")
  makeMatrix(allLengths,teRNA,teUTR,teUORF)
}
#Check file format exists
findFF = function(formatName, boolreturn = F){
  
  regEx = paste("\\.",formatName,"$",sep = "")
  ffList = list.files(pattern = regEx, ignore.case=TRUE)
  if(length(ffList) == 1){
    if(boolreturn)return(T)
    return(ffList)
    
  } else{
    if(boolreturn)return(F)
    cat("no ", formatName,"file in this folder")
    stop()
  }
}
#####################################GETTESS##########################################
#seq: sequence, 
#dm1: detection method 1,
#dm2: detection method 2
#seqLength: length of seq
#al: length of all sequences, might be NULL
getTE = function(seq, dm1,dm2,seqLength1,seqLength2,al = NULL){
  
  overlap1 = countOverlaps(unlist(seq),dm1)
  overlap2 = countOverlaps(unlist(seq),dm2)
  
  something1 = aggregate(overlap1,by = list(names(overlap1)), sum)
  something2 = aggregate(overlap2,by = list(names(overlap2)), sum)
  
  librarydm1 = length(dm1)
  librarydm2 = length(dm2)
  #norm1: rna normalization, norm2: rfp normalization
  if(is.null(al)){
    norm1 = FPKMNorlization(something1$x,seqLength1,librarydm1)
    norm2 = FPKMNorlization(something2$x,seqLength2,librarydm2)
    assign("normUORFRNA",norm1, envir = .GlobalEnv)
    assign("normUORFRFP",norm2, envir = .GlobalEnv)
  }else{
    norm1 = FPKMNorlization(cp(something1$x,al),unlist(seqLength1),librarydm1)
    norm2 = FPKMNorlization(cp(something2$x,al),unlist(seqLength2),librarydm2)
    if(length(seq)==length(cds)){
      assign("normCDSRNA",norm1, envir = .GlobalEnv)
      assign("normCDSRFP",norm2, envir = .GlobalEnv)
    }else{
      assign("normUTRRNA",norm1, envir = .GlobalEnv)
      assign("normUTRRFP",norm2, envir = .GlobalEnv)
    }
  }
  te = norm2/norm1
  names(te) = something1[,1]
  te = as.matrix(te)
  return(te)
}
#correct positions al: trainscript lengths
cp = function(overlap,al){
  return(unlist(overlap[match(as.character( al$tx_name),names(overlap))]))
}
makeMatrix = function(allLengths,teRNA,teUTR,teUORF){
  matrixA = cbind(allLengths,teRNA,teUTR)
  
  transcriptNames = gsub("_[0-9]*","", rownames(teUORF))
  matrixB = cbind(transcriptNames,teUORF)
  matrix = merge(matrixA, matrixB, by.x = "tx_name", by.y = "transcriptNames", all = T)
  
  index = matrix[,9] == 0 | is.na(matrix[,9] ) | is.infinite(matrix[,9] ) | is.nan(matrix[,9] )
  matrix = matrix[!index,]
  
  
  #matrix = removeNAN(matrix,"normRNA")
  #te remove
#   matrix = removeNAN(matrix,"teUTR")
#   matrix = removeNAN(matrix,"teRNA")
#   matrix = removeNAN(matrix,"teUTR",saveIndex = T)
  
  write.csv(matrix, file = "matrixtest.csv")
}


getMatrix(arcs[1],saveToGlobal = T)

