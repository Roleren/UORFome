###Only used in combination with uorfomeHelper



#'Make general name of experiment run
makeGeneralName = function(cageName = NULL,leaderBed = NULL,rnaSeq = NULL,rfpSeq = NULL){
  
  if(!is.null(cageName)){ #General name
    thisCage = getRelativePathName(cageName)
    assign("thisCage",thisCage,envir = .GlobalEnv) ###This is dangerous,something wrong here!
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
    assign("detailedFullName",detailedFullName,envir = .GlobalEnv)
    cat("detailed Full Name of run is: ",detailedFullName,"\n")
  }
}

### Print info about specific run
infoPrints = function(doubleBAM,usingNewCage,cageName,leaderBed,rnaSeq = NULL,rfpSeq = NULL){
  print("Estimated normal time is 5 hours\n")
  cat("data output folder used is:",normalizePath(dataFolder),"\n")
  if(doubleBAM == T)
    cat("using bam-file for RFP\n")
  if(usingNewCage)
    cat("using cage-file named: \n",cageName,"\n")
  
  makeGeneralName(cageName,leaderBed,rnaSeq,rfpSeq)
  cat("starting loading objects\n")
}

#'Get name without extension, update this!!
getCleanName = function(loadPath,nameUsed){
  if(!is.null(loadPath)){
    fName = gsub(".*/", "", loadPath)
  }else
    fName = gsub(".*/", "", nameUsed)
  return( gsub(pattern =".rdata",replacement = "", x = fName) )# get name without .rdata
}

getFastaName = function(cleanName){
  return( paste0(fastaFolder,cleanName,".fasta") ) #For fasta
}

getUORFBedName = function(cleanName){
  return( paste0(uorfBedFolder,cleanName,".bed") )# for uorfBed
}
getUORFRDataName = function(cleanName){
  return( paste0(uorfFolder,getRelativePathName(cleanName),".uorf.rdata") )# for uorfBed
}

chooseUORFName = function(loadPath = NULL, nameSave = NULL){
  if( !is.null(nameSave)){
    nameU = nameSave
  }else if( exists("thisCage") && !is.null(thisCage)){
    nameU = thisCage
  }else if(is.null(loadPath) && exists("generalName")){ #fix this!!!!
    nameU = generalName
  }else if(exists("cageName") &&is.null(cageName)){
    nameU = generalName
  }else if(!is.null(loadPath)){
    nameU = loadPath
  }else{
    stop("NO name for uorfs bed and fasta files!")
  }
  nameU
}

#remove extensions, _x and dots .
getTranscriptNames = function(txsnames){
  newNames =  gsub("_[0-9]*","", txsnames)
  return ( gsub("\\..*","", newNames) ) 
}
