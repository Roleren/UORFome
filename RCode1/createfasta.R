arcsCFB = commandArgs(trailingOnly = T)

source("/export/valenfs/projects/uORFome/RCode1/uorfomeGeneratorHelperFunctions.R")
source("/export/valenfs/projects/uORFome/RCode1/HelperVariables.R")

#always define loadPath if rangesOfUorfs is not in global scope!
createFastaAndBedFile = function(loadPath ="/export/valenfs/projects/uORFome/test_results/rangesOfUORFs/CD34%2b%20stem%20cells%20-%20adult%20bone%20marrow%20derived%2c%20donor1%2c%20tech_rep2.CNhs12553.12225-129F2.RangesUorf.rdata",nameUsed = "rangesOfUorfs"){
  print("Starting making fasta and Bed files")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)
  
  getFasta() #get fasta and fai
  
  print("loading finished, now making fasta list")
  fastaList = lapply(X = 1:length(rangesOfuORFs), FUN = findFasta) #Get list of all fasta sequences of uorfs
  print("fasta list is finished")
  #make correct output name
  if(!is.null(loadPath))
    fName = gsub(".*/", "", loadPath)
  else
    fName = gsub(".*/", "", nameUsed)
  fName = strsplit(fName,".rdata")[[1]][1]# used for uorf too
  ffName = paste0(fastaFolder,fName,".fasta",sep = "") #For fasta
  bName =  paste0(uorfBedFolder,fName,".bed",sep = "") #for bed uorfs
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
  
  bedColumns = df.rangesOfuORFs[,c("seqnames","start","end",score = "width","strand","names")]
  
  if(is.null(bedName))
    write.table(x = bedColumns,file = paste0(uorfBedFolder,"rangesOfUorfsLight.bed") ,sep = "\t",col.names = F,row.names = F, quote = F)
  else
    write.table(x = bedColumns,file = bedName ,sep = "\t",col.names = F,row.names = F, quote = F)
}

if(length(arcsCFB) == 1)
  createFastaAndBedFile(loadPath = normalizePath(arcsCFB[1]))
if(length(arcsCFB) == 2)
  createFastaAndBedFile(loadPath = normalizePath(arcsCFB[1]),nameUsed =  normalizePath(arcsCFB[2]))

