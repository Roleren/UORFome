#arcsCFB = commandArgs(trailingOnly = T)




#always define loadPath if rangesOfUorfs is not in global scope!
createFastaAndBedFile = function(loadPath =NULL,nameUsed = "rangesOfUorfs"){
  print("Starting making fasta and Bed files")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)
  
  getFasta() #get fasta and fai
  
  print("loading finished, now making fasta list")
  fastaList = lapply(X = 1:length(rangesOfuORFs), FUN = findFasta) #Get list of all fasta sequences of uorfs
  print("fasta list is finished")
  
  #make correct output name
  fName = getCleanName(loadPath,nameUsed)
  ffName = getFastaName(fName)
  bName =  getUORFBedName(fName) #for bed uorfs
  
  createUORFBedFile(bName) #create ranges bed file
  cat("new fasta name is: ", ffName)
  
  #write out fastalist
  
  stringSet = unlist(DNAStringSetList(fastaList))
  
  writeXStringSet(stringSet, ffName)
  
  return(stringSet)
}

findFasta = function(i){
  f = rangesOfuORFs[i]
  transcriptName = names(rangesOfuORFs[[i]])[1]
  dna <- as.character(unlist(getSeq(fa, unlist(f))))
  dna <- paste(dna,collapse = "")
  
  dna = DNAStringSet(dna)
  names(dna) = transcriptName
  dna
}

#bName is bed name
#chr, start,end,nameUORF,score = length? or ., strand
createUORFBedFile = function(bedName = NULL){
  cat("writing bed file with name: ", bedName)
  
  df.rangesOfuORFs = as.data.frame(rangesOfuORFs)
  
  bedColumns = df.rangesOfuORFs[,c("seqnames","start","end","names",score = "width","strand")]
  
  if(is.null(bedName))
    write.table(x = bedColumns,file = paste0(uorfBedFolder,"rangesOfUorfsLight.bed") ,sep = "\t",col.names = F,row.names = F, quote = F)
  else
    write.table(x = bedColumns,file = bedName ,sep = "\t",col.names = F,row.names = F, quote = F)
}

# if(length(arcsCFB) == 1)
#   createFastaAndBedFile(loadPath = normalizePath(arcsCFB[1]))
# if(length(arcsCFB) == 2)
#   createFastaAndBedFile(loadPath = normalizePath(arcsCFB[1]),nameUsed =  normalizePath(arcsCFB[2]))

