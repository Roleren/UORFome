#Creates normalizations of the counts
FPKMNorlization = function(counts, lengthSize, librarySize){
  result = (as.numeric(counts)*(10^9)) / (as.numeric(lengthSize)*as.numeric(librarySize))
  return(result)
}

#as number coverter function
an = function(fac){
  return(as.numeric(as.character(fac)))
}

#check that bamFile is paired or not, False = not paired
testBAM = function(name){
  testPairedEndBam(name)
}

#correct positions al: trainscript lengths
cp = function(overlap,al){
  return(unlist(overlap[match(as.character( al$tx_name),names(overlap))]))
}

#Check if file format exists.
#if boolreturn == F, find the file name
#if boolreturn == T, return T.
#bamType, RNA or RFP/RPF
findFF = function(formatName, boolreturn = F, bamType = NULL){
  
  regEx = paste("\\.",formatName,"$",sep = "")
  ffList = list.files(pattern = regEx, ignore.case=TRUE)
  if(length(ffList) == 1){
    if(boolreturn)return(T)
    return(ffList)
    
  }else if((formatName == "bam") & (length(ffList) == 2) & (!is.null(bamType))){
    #find all sequences with bamtype in string  
    ffList1 = ffList[lapply(ffList,function(x) length(grep(bamType,x,value=F))) == 1]
    
    if(bamType == "RFP" & length(ffList1) == 0){
      ffList1 = ffList[lapply(ffList,function(x) length(grep("RPF",x,value=F))) == 1]
      if(length(ffList1) == 0){
        print("Error: rename the RFP.bam file to RFP.bam")
        stop()
      }
    }
    cat("using as ",bamType, "the file: ", unlist(ffList),"\n")
    if(length(ffList) > 0)
      return (ffList[1])
    else
      return(ffList1[1])
    
  } else{
    if(boolreturn)return(F)
    cat("No ", formatName,"file in this folder or duplicates")
    stop()
  }
}

getRelativePathName = function(name){
  return (gsub(".*/", "", name))
}

rfe = function(name){#Only works on .bam right now!!!!
  return (gsub("*\\.bam", "", name))
}

#Check if uorfRanges exist already, or if must be created.
###########Should make this more failsafe!!!!!!!!!!
UorfRangesNotExists = function(){
  if(exists("rangesOfuORFs") == F){
    
    if(findFF("rdata",boolreturn = T)){
      cat("loading rangesOfuorf from folder\n",getwd())
      load(findFF("rdata"),envir = .GlobalEnv)
      
      return(F)
    }else{return(T)}
  }
  return(F)
}



