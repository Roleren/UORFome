arcsRFU = commandArgs(trailingOnly = T)

source("./createfasta.R")
library(GenomicFeatures)
library(GenomicAlignments)
require(data.table)


###Find all uorf into the cds,filter bad ones
#To bed added!: check if reading frame changes. %3 = 0
# If called without path, need rangesofUORFs in global scope!!

removeFalseUORFs = function(loadPath = NULL,saveToFile = F,outputFastaAndBed = F){
  #rm(list = ls())
  print("starting to filter out bad ourfs...")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)
  if(exists("Gtf") == F){
    print("loading GTF")
    getGTF()
    cds = cdsBy(Gtf,"tx",use.names = T)
  }
  ####Save overlaps of cds and uorfs for plotting later
  overlap1 = findOverlaps(cds,rangesOfuORFs)
  overlapCount = countOverlaps(cds,rangesOfuORFs)
  numberOfOverlaps = sum(overlapCount >= 1)
  overlapHitsIndex = overlapCount[overlapCount == 1]
  assign("numberOfOverlaps",numberOfOverlaps,envir = .GlobalEnv)
  
  #get start cds, end uorf, find difference, check mod3
  #check start upstream, ending downstream
  df.rangesOfuORFs = as.data.frame(rangesOfuORFs)
  df.rangesOfuORFsOnlyFirstExon = setDT(df.rangesOfuORFs)[, if(.N>1) if(strand=="+") head(.SD, 1) else tail(.SD,1) else .SD , by = names]
  df.rangesOfuORFsOnlyFirstExon = as.data.frame(df.rangesOfuORFsOnlyFirstExon)
  df.pos = df.rangesOfuORFsOnlyFirstExon[df.rangesOfuORFsOnlyFirstExon$strand == "+",]
  df.neg = df.rangesOfuORFsOnlyFirstExon[df.rangesOfuORFsOnlyFirstExon$strand == "-",]
  df.pos.finished = df.pos[,cbind("names","start","strand","seqnames")]
  df.neg.finished = df.neg[,cbind("names","end","strand","seqnames")]
  
  df.pos.finished$uorfnames = df.pos.finished$names
  df.neg.finished$uorfnames = df.neg.finished$names 
  
  transcriptNames.pos = gsub("_[0-9]*","", df.pos.finished$names)
  transcriptNames.neg = gsub("_[0-9]*","", df.neg.finished$names)
  transcriptNames.pos = gsub("\\..*","", df.pos.finished$names)
  transcriptNames.neg = gsub("\\..*","", df.neg.finished$names)
  
  df.pos.finished$names = transcriptNames.pos
  df.neg.finished$names = transcriptNames.neg
  
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
  
  test = df.removedBadUORFs
  
  test = test[,-"group_name"]
  test = test[,-"group"]
  uorfName.bad = gsub(".*\\.","", df.removedBadUORFs$names)
 
  assign("uorfName.bad",uorfName.bad,envir = .GlobalEnv)
  assign("test",test,envir = .GlobalEnv)
  #Create uorf IDs
  uID = with(test,createUorfIDs(seqnames,start,end,names,width,strand))
  test$uID = uID
  
  #Make object to GrangesList
  uniqueUORFs = unique(uorfName.bad)
  
  test1 = lapply(uniqueUORFs, function(x) getGRLbyName(x))
  test2 = GRangesList(test1)
  rangesOfuORFs = test2
  
  
  if(saveToFile)
    save(rangesOfuORFs, file = loadPath)
  
  print("finished filtering bad ourfs")
  ####Check for frame shifted uorfs
  if(is.null(loadPath)){ #fix this!!!!
    nameU = generalName
  }else{
    nameU = loadPath
  }
  if(outputFastaAndBed)
    createFastaAndBedFile(loadPath = NULL,nameUsed = nameU)
  
  return(rangesOfuORFs)
  
}

#for each unique
#extract 
#return granges
getGRLbyName = function(x){
  tempEquals = test[uorfName.bad == x]
  a = GRanges(tempEquals)
  names(a) = rep(x,nrow(tempEquals))
  return(a)
}

createUorfIDs = function(seqnames,start,end,names,width,strand){
  paste0(seqnames,";",start,";",end,";",names,";",width,";",strand)
}


if(length(arcsRFU) == 2){
  removeFalseUORFs(loadPath = normalizePath(arcsRFU[1]),saveToFile = as.logical(arcsRFU[2]))
}else if(length(arcsRFU) == 3){
  removeFalseUORFs(loadPath = normalizePath(arcsRFU[1]),saveToFile = as.logical(arcsRFU[2]),outputFastaAndBed = as.logical(arcsRFU[3]))
}

