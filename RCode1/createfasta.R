#arcsCFB = commandArgs(trailingOnly = T)




#always define loadPath if rangesOfUorfs is not in global scope!
createFastaAndBedFile = function(rangesOfuORFs,loadPath =NULL,nameUsed = "rangesOfUorfs"){
  print("Starting making fasta and Bed files")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)

  getSequencesFromFasta(rangesOfuORFs)
  
  print("fasta list of uorfs is finished")
  
  #make correct output name
  fName = getCleanName(loadPath,nameUsed)
  ffName = getFastaName(fName)
  bName =  getUORFBedName(fName) #for bed uorfs
  
  createUORFBedFile(rangesOfuORFs,bName) #create ranges bed file
  cat("new fasta name is: ", ffName)
  
  #write out fastalist
  
  stringSet = unlist(DNAStringSetList(seqs))
  
  writeXStringSet(stringSet, ffName)
  
  return(stringSet)
}

#bName is bed name
#chr, start,end,nameUORF,score = length? or ., strand
createUORFBedFile = function(rangesOfuORFs,bedName = NULL){
  cat("writing bed file with name: ", bedName)
  
  df.rangesOfuORFs = as.data.frame(rangesOfuORFs)
  
  bedColumns = df.rangesOfuORFs[,c("seqnames","start","end","names",score = "width","strand")]
  
  if(is.null(bedName))
    write.table(x = bedColumns,file = paste0(uorfBedFolder,"rangesOfUorfsLight.bed") ,sep = "\t",col.names = F,row.names = F, quote = F)
  else
    write.table(x = bedColumns,file = bedName ,sep = "\t",col.names = F,row.names = F, quote = F)
}
