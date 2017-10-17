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

chooseUORFName = function(loadPath){
  if( exists("thisCage") && !is.null(thisCage)){
    nameU = thisCage
  }else if(is.null(loadPath) && exists("generalName")){ #fix this!!!!
    nameU = generalName
  }else{
    nameU = loadPath
  }
  nameU
}

#remove extensions, _x and dots .
getTranscriptNames = function(txsnames){
  newNames =  gsub("_[0-9]*","", txsnames)
  return ( gsub("\\..*","", newNames) ) 
}
