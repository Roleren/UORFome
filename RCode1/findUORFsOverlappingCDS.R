#arcsRFU = commandArgs(trailingOnly = T)

source("./createfasta.R")
library(GenomicFeatures)
library(GenomicAlignments)
require(data.table)


###Find all uorf into the cds,filter bad ones

# If called without path, need rangesofUORFs in global scope!!

removeFalseUORFs = function(loadPath = NULL,saveToFile = F,outputFastaAndBed = F){
  
  ###########PRE LOADINGS#############
  print("starting to filter out bad ourfs...")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)
  if(exists("Gtf") == F){
    print("loading GTF")
    getGTF()
    cds = cdsBy(Gtf,"tx",use.names = T)
  }
 
  
  ################FILTERING############
  
  #get start cds, end uorf, find difference, check mod3 ?
  #check start upstream, ending downstream
  df.removedBadUORFs = removeUORFsThatAreAcctualyCDSs(rangesOfuORFs,cds)
  
  #Make object to GrangesList
  ####Check for frame shifted uorfs ???
  rangesOfuORFs = filterOutDuplicateUORFColumns(df.removedBadUORFs)
  print("finished filtering bad ourfs")
  
  ################SAVING###############
  ####Save overlaps of cds and uorfs for plotting later
  saveOverlaps(rangesOfuORFs,cds)
  
  
  
  if(saveToFile){
    cat("Saving uorfs rdata to: ",loadPath)
    save(rangesOfuORFs, file = loadPath)
  }
  if(outputFastaAndBed){
    nameU = chooseUORFName(loadPath)
    createFastaAndBedFile(loadPath = NULL,nameUsed = nameU)
  }
  return(rangesOfuORFs)
}
###TODO: Check if this method can be removed by countoverlaps with type = equal ? 
removeUORFsThatAreAcctualyCDSs = function(rangesOfuORFs,cds){
  #Split + / - strands
  df.rangesOfuORFs = as.data.frame(rangesOfuORFs)
  df.rangesOfuORFsOnlyFirstExon = setDT(df.rangesOfuORFs)[, if(.N>1) if(strand=="+") head(.SD, 1) else tail(.SD,1) else .SD , by = names]
  df.rangesOfuORFsOnlyFirstExon = as.data.frame(df.rangesOfuORFsOnlyFirstExon)
  df.pos = df.rangesOfuORFsOnlyFirstExon[df.rangesOfuORFsOnlyFirstExon$strand == "+",]
  df.neg = df.rangesOfuORFsOnlyFirstExon[df.rangesOfuORFsOnlyFirstExon$strand == "-",]
  df.pos.finished = df.pos[,cbind("names","start","strand","seqnames")]
  df.neg.finished = df.neg[,cbind("names","end","strand","seqnames")]
  
  #Set uorf names
  df.pos.finished$uorfnames = df.pos.finished$names
  df.neg.finished$uorfnames = df.neg.finished$names 
  
  #Set transcript names
  df.pos.finished$names = getTranscriptNames(df.pos.finished$names)
  df.neg.finished$names = getTranscriptNames(df.neg.finished$names)
  
  #same for cds
  df.cds = as.data.frame(cds)
  df.cds.onlyFirstExon = setDT(df.cds)[, if(.N>1) if(strand=="+") head(.SD, 1) else tail(.SD,1) else .SD , by = group_name]
  df.cds.onlyFirstExon = as.data.frame(df.cds.onlyFirstExon)
  
  df.cds.pos = df.cds.onlyFirstExon[df.cds.onlyFirstExon$strand == "+",]
  df.cds.neg = df.cds.onlyFirstExon[df.cds.onlyFirstExon$strand == "-",]
  
  ### Check for fake uorfs, that are acctualy just the first exon atg -> stop codon
  ### Only happens on 1 exon transcripts
  
  posMerge = merge(df.cds.pos,df.pos.finished[,c("names","start","strand","seqnames","uorfnames")],by.x="group_name",by.y="names")
  posMergesresult = posMerge[which(posMerge$start.x == posMerge$start.y),]
  
  negMerge = merge(df.cds.neg,df.neg.finished[,c("names","end","strand","seqnames","uorfnames")],by.x="group_name",by.y="names")
  negMergesresult = negMerge[which(negMerge$end.x == negMerge$end.y),]
  
  ### Remove fake uorfs
  posOrfsToBeRemoved = posMergesresult$uorfnames
  negOrfsToBeRemoved = negMergesresult$uorfnames
  
  df.removedBadUORFs  = df.rangesOfuORFs
  ###Remove positive ones
  for(i in posOrfsToBeRemoved){
    df.removedBadUORFs = df.removedBadUORFs[df.removedBadUORFs$names != i]
  }
  ###Remove negative ones
  for(i in negOrfsToBeRemoved){
    df.removedBadUORFs = df.removedBadUORFs[df.removedBadUORFs$names != i]
  }
  print("Removed uorfs that were acctualy cdss")
  return(df.removedBadUORFs)
}

filterOutDuplicateUORFColumns = function(df.removedBadUORFs){
  print("filtering out duplicate uorf columns")
  uorfs = df.removedBadUORFs[,-"group_name"]
  uorfs = uorfs[,-"group"]
  assign("uorfs",uorfs,envir = .GlobalEnv )
  
  uorfNames= gsub(".*\\.","", df.removedBadUORFs$names)
  assign("uorfNames",uorfNames,envir = .GlobalEnv )
  uniqueUORFs = unique(uorfNames)
  
  #Create uorf IDs
  uID = with(uorfs,createUorfIDs(seqnames,start,end,names,width,strand))
  uorfs$uID = uID
  
  uorfs = lapply(uniqueUORFs, function(x) getGRLbyName(x))
  uorfs = GRangesList(uorfs)
  print("finished filtering  duplicates")
  return( uorfs )
}

#this is not finished yet, will this be the final id ?
createUorfIDs = function(seqnames,start,end,names,width,strand){
  paste0(seqnames,";",start,";",end,";",names,";",width,";",strand)
}

#for each unique
#extract 
#return granges
getGRLbyName = function(x){
  tempEquals = uorfs[uorfNames == x]
  a = GRanges( tempEquals )
  names(a) = rep(x,nrow(tempEquals))
  return(a)
}

saveOverlaps = function(rangesOfuORFs,cds){
  numberOfOverlaps = getUOrfOverlaps()
  assign("numberOfOverlaps",numberOfOverlaps,envir = .GlobalEnv)
}




#if(length(arcsRFU) == 2){
#  removeFalseUORFs(loadPath = normalizePath(arcsRFU[1]),saveToFile = as.logical(arcsRFU[2]))
#}else if(length(arcsRFU) == 3){
#  removeFalseUORFs(loadPath = normalizePath(arcsRFU[1]),saveToFile = as.logical(arcsRFU[2]),outputFastaAndBed = as.logical(arcsRFU[3]))
#}

