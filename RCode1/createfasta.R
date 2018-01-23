
#always define loadPath if rangesOfUorfs is not in global scope!
createFastaAndBedFile <- function(rangesOfuORFs,loadPath =NULL,nameUsed = "rangesOfUorfs"){
  print("Starting making fasta and Bed files")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)

  getSequencesFromFasta(rangesOfuORFs)
  
  print("fasta list of uorfs is finished")
  
  #make correct output name
  Name = getCleanName(loadPath,nameUsed)
  fName = getFastaName(Name)
  bName =  getUORFBedName(Name) #for bed uorfs
  
  bed6(rangesOfuORFs,bName) #create ranges bed file
  cat("new fasta name is: ", fName)
  
  #write out fastalist
  
  stringSet = unlist(DNAStringSetList(seqs))
  
  writeXStringSet(stringSet, fName)
  
  return(stringSet)
}


