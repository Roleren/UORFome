#filter
#flank
#overlap
#for all that overlap, assign strongest peak as tss
#cageR package ?
#plot for 1 transcript, see peaks?

toGR=function(x,bed6=TRUE){
  require(GenomicRanges)
  if(!bed6){
    gr<- GRanges(x[,1],IRanges(x[,2]-1,x[,3]))
    return(gr)
  }
  starts<- ifelse(x[,6]=="+",x[,2]-1,x[,2])
  ends<- ifelse(x[,6]=="-",x[,3],x[,3]-1)
  gr<- GRanges(x[,1],IRanges(starts,ends))
  strand(gr)=x[,6]
  score(gr)=x[,5]
  if(ncol(x)>6)  mcols(gr)=x[,7:ncol(x)]
  return(gr)
}
#get CageData as granges
getFilteredCageData = function(dataName, filterValue = 1){
  rawCageData = toGR(as.data.frame(fread(paste("gunzip -c",dataName),sep = "\t")))
  getCDS()
  print("finished loading cage file")
  filteredrawCageData = rawCageData[rawCageData$score > filterValue,] #filter on score 1
  seqlevels(filteredrawCageData) = sub("chrY","Y",seqlevels(filteredrawCageData))
  seqlevels(filteredrawCageData) = sub("chrX","X",seqlevels(filteredrawCageData))
  seqlevels(filteredrawCageData) = sub("chrM","MT",seqlevels(filteredrawCageData))
  assign("filteredrawCageData",filteredrawCageData,envir = .GlobalEnv)
}

addFirstCdsOnLeaderEnds = function(fiveUTRs){
  getCDS()
  cdsForUTRs <- cds[names(fiveUTRs)] # get only the ones we need
  firstExons <- phead(cdsForUTRs, 1L) #select first in every, they must be sorted!
  gr = unlist(firstExons,use.names = F)
  gr$cds_id = NULL; gr$cds_name = NULL; gr$exon_rank = NULL 
  gr$exon_id = NA;  gr$exon_name = NA;  gr$exon_rank = NA
  grl = relist(gr,firstExons)
  fiveUTRsWithCdsExons <- pc(fiveUTRs, grl)
  reduce(fiveUTRsWithCdsExons)# ask gunnar why length != fiveUTRs
  
  return( reduce(fiveUTRsWithCdsExons) ) 
}

#Extend first exon of each transcript with 1000
# warning! does not check position < 1 and > seqlength
extendsTSSExons = function(fiveUTRs, extension = 1000){
  fiveAsgr = unlist(fiveUTRs)
  firstExons = fiveAsgr[fiveAsgr$exon_rank == 1]
  
  posIDs = firstExons[strand(firstExons) == "+"]
  minIDs = firstExons[strand(firstExons)  == "-"]
  
  firstExons[names(posIDs)] = resize(firstExons[names(posIDs)],  width = width(firstExons[names(posIDs)])+extension, fix = 'end')
  firstExons[names(minIDs)] = resize(firstExons[names(minIDs)],  width = width(firstExons[names(minIDs)])+extension, fix = 'start')
  return(firstExons)
}

#Find max peak for each transcript,
# returns as data.table, without names, but with index
findMaxPeaks = function(cageOverlaps){
  dt = as.data.table(filteredrawCageData)
  dt = dt[from(cageOverlaps)]
  dt$to = to(cageOverlaps)
  
  test = dt[,max(score), by = to]
  names(test) = c("to","score")
  res =  merge(  test, dt) 
  
  return(res[!duplicated(res$to)]) #check this line!!!
}
#fiveAsgr[fiveAsgr$exon_rank == 1] = firstExonsShifted
#return(relist(fiveAsgr,fiveUTRs))
#load leader, cds,
#extend leader 1000 upstream, and to end of 1st cds downstream
# look for cage peaks in leader upstream
findNewTSS = function(fiveUTRs,dataName){
  if(exists("maxPeakPosition") == F){
   
    getFilteredCageData(dataName) # get the cage data
    
    shiftedfiveUTRs = extendsTSSExons(fiveUTRs)
    assign("shiftedfiveUTRs",shiftedfiveUTRs,envir = .GlobalEnv)
    cageOverlaps = findOverlaps(query = filteredrawCageData,subject = shiftedfiveUTRs)
    
    maxPeakPosition = findMaxPeaks(cageOverlaps)
    assign("maxPeakPosition",maxPeakPosition,envir = .GlobalEnv)
    
    print("found new cage peaks")
  }
}

makeGrlAndFilter = function(firstExons, fiveUTRs){
  fiveAsgr = unlist(fiveUTRs)
  fiveAsgr[fiveAsgr$exon_rank == 1] = firstExons
  return(relist(fiveAsgr,fiveUTRs))
}

#add cage max peaks as new tss
addNewTssOnLeaders = function(fiveUTRs){
  fiveAsgr = unlist(fiveUTRs)
  firstExons = fiveAsgr[fiveAsgr$exon_rank == 1]
  
  maxPeakPosition$names = names(firstExons[maxPeakPosition$to])
  posIDs = maxPeakPosition$to[maxPeakPosition$strand == "+"]
  minIDs = maxPeakPosition$to[maxPeakPosition$strand == "-"]
  
  firstExons[posIDs] = resize(firstExons[posIDs], width = end(firstExons[posIDs]) - maxPeakPosition$start[maxPeakPosition$strand == "+"] + 1, fix = "end")
  firstExons[minIDs] = resize(firstExons[minIDs], width = end(firstExons[minIDs]) - maxPeakPosition$start[maxPeakPosition$strand == "-"] + 1, fix = "start")
  
  return( firstExons ) 
  
}

#find and add cage max peaks as new tss's
# if a max peak > filter_size is found, it is reassigned, else it is kept
###NB! Must have Gtf in global scope
getNewfivePrimeUTRs = function(fiveUTRs,dataName = standardCage){
  ###Read in cage files
  print(dataName)
  
  findNewTSS(fiveUTRs,dataName)
  fiveUTRs = makeGrlAndFilter(addNewTssOnLeaders(fiveUTRs), fiveUTRs)
  fiveUTRs = addFirstCdsOnLeaderEnds(fiveUTRs)
  return(fiveUTRs)
}
#Since cage redefines five utr lengths, the total lengths of tx must be updated
findCageUTRFivelen = function(fiveUTRs,oldTxNames){
  newfiveprimeLen = widthPerGRangesGroup(fiveUTRs)
  return( newfiveprimeLen[match(oldTxNames,names(newfiveprimeLen))])
}

otherSpeciesCageCSV = function(name){
  #pre loadings
  #name = "/export/valenfs/data/processed_data/CAGE/nepal_2013_zebrafish/called_peaks/leaders_zv10_zf_02_fertilized_egg.csv"
  cage = read.csv(name)
  cage = cage[cage$count_at_highest_peak != 0,]
  
  fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)

  #make max peak object
  cageGR = GRanges(seqnames = gsub("chr","",as.character(cage$chr)), ranges = IRanges(start = cage$highest_peak,end = cage$highest_peak),strand = cage$dir, names = cage$X.gene_id)
  maxPeakPosition = as.data.table(findOverlaps(query = cageGR,subject = extendsTSSExons(fiveUTRs)))
  names(maxPeakPosition) = c("from","to")
  maxPeakPosition = maxPeakPosition[!duplicated(maxPeakPosition$to)]
  maxPeakPosition$strand = as.character(strand(cageGR[maxPeakPosition$from]))
  maxPeakPosition$start = start(cageGR[maxPeakPosition$from])
  maxPeakPosition$end = end(cageGR[maxPeakPosition$from])
  assign("maxPeakPosition",maxPeakPosition,envir = .GlobalEnv)
  
  #make new fiveUTRs
  fiveUTRs = makeGrlAndFilter(addNewTssOnLeaders(fiveUTRs), fiveUTRs)
  fiveUTRs = addFirstCdsOnLeaderEnds(fiveUTRs)
  
  unlistfgr = unlist(fiveUTRs)
  #On zebra fish, the genome and gtf had different seqname naming, so must match!
  seqnamesTransformed = as.character(seqnames(unlistfgr))
  seqnamesTransformed[nchar(seqnamesTransformed) > 4] = paste0("Un_",seqnamesTransformed[nchar(seqnamesTransformed) > 4])
  seqnamesTransformed = gsub("\\.","v",seqnamesTransformed)
  seqnamesTransformed = paste0("chr",seqnamesTransformed)
  unlistfgr = GRanges(seqnames = seqnamesTransformed, ranges = IRanges(start = start(unlistfgr),end = end(unlistfgr)),strand = strand(unlistfgr))
  names(unlistfgr) = names(unlist(fiveUTRs))
  fiveUTRs = GroupGRangesByNames(unlistfgr)
  
  return(fiveUTRs)
}

