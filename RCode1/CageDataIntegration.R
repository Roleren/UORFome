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

#Find max peak for each transcript
findMaxPeaks = function(x){
  vec = (cageOverlaps@to == x)
  if(length(vec) == 0)
    return(NA)
  vec2 = (cageOverlaps@from[vec])
  df = filteredrawCageData[vec2]
  vec3 = which.max(score(df))
  return(start(df)[vec3])
}

 
optimizefiveUTRs = function(x){
  if(length(maxPeakPosition[[x]])  != 0){
    lapply(1:length(optimizedfiveUTRs[[x]]), function(x) modifyUTR)
    return(resize(optimizedfiveUTRs[x],  width = maxPeakPosition[[x]], fix = 'end')) 
  }
  else{
    return(NA)
  }
}

insertFirstCDS = function(newUTR,x){
  
  firstExon = unlist(cds[names(cds) == names(shiftedfiveUTRs[x])])[1]
  fiveTemp = newUTR
  #transcriptName = names(newUTR)
  mcols(fiveTemp) = NULL
  mcols(firstExon) = NULL
 
  combination = c(fiveTemp,firstExon)
  combination = sort(combination)
  
  newUTR = combination
  return(newUTR)
}

#set new tss to max peak
modifyUTR = function(x){
  if(length(maxPeakPosition[[x]])  != 0){
    indexOfOverlaptemp = findOverlaps(shiftedfiveUTRs[[x]],GRanges(seqnames = seqnames(shiftedfiveUTRs[[x]][1]),IRanges(start =maxPeakPosition[[x]], end = maxPeakPosition[[x]]+1)))
    indexOfOverlap = indexOfOverlaptemp@from
    newUTR = shiftedfiveUTRs[[x]]
    
    if(as.character(newUTR@strand)[1] == "+"){
      start(newUTR[indexOfOverlap]) = maxPeakPosition[[x]]
      if(indexOfOverlap > 1)
        newUTR = newUTR[indexOfOverlap:length(newUTR)]
      }
    else{
      end(newUTR[indexOfOverlap]) = maxPeakPosition[[x]]
      if(indexOfOverlap > 1)
        newUTR = newUTR[indexOfOverlap:length(newUTR)]
    }
    newUTR = insertFirstCDS(newUTR,x)
    return(newUTR)
  }
  return(NA)
}
#load leader, cds,
#extend leader 1000
# look for cage peaks in leader upstream
findNewTSS = function(fiveUTRs,dataName){
  if(exists("optimizedfiveUTRs") == F){
    #rawCageData = read.table(dataName,sep = "\t")
    rawCageData = as.data.frame(fread(paste("gunzip -c",dataName),sep = "\t"))
    rawCageData = toGR(rawCageData)
    getCDS()
    print("finished loading cage file")
    filteredrawCageData = rawCageData[rawCageData$score > 1,] #filter on score 1
    assign("filteredrawCageData",filteredrawCageData,envir = .GlobalEnv)
    
    namesUTRs = names(fiveUTRs)
    shiftedfiveUTRs <- resize(fiveUTRs,  width = width(fiveUTRs)+1000, fix = 'end')
    #shiftedfiveUTRs <- resize(fiveUTRs,  width = width(fiveUTRs)+200, fix = 'start')
    assign("shiftedfiveUTRs",shiftedfiveUTRs,envir = .GlobalEnv)
    cageOverlaps = findOverlaps(query = filteredrawCageData,subject = shiftedfiveUTRs)
    assign("cageOverlaps",cageOverlaps,envir = .GlobalEnv)
    
    maxPeakPosition = lapply(1:length(shiftedfiveUTRs),function(x) findMaxPeaks(x))
    
    assign("maxPeakPosition",maxPeakPosition,envir = .GlobalEnv)
    optimizedfiveUTRs = shiftedfiveUTRs
    assign("optimizedfiveUTRs",optimizedfiveUTRs,envir = .GlobalEnv)
    print("found new cage peaks")
  }
}
#add cage max peaks as new tss
addNewTssOnLeaders = function(fiveUTRs){
  print("now assigning new tss's")
  goal = lapply(1:length(shiftedfiveUTRs), function(x) modifyUTR(x))
  print("finished assigning new tss's")
  assign("goal",goal,envir = .GlobalEnv)
  names(goal) = names(fiveUTRs)
  goal2 = !is.na(goal)
  
  newUTRs = GRangesList(goal[goal2])
  # Filter out UTRs with length 0
  newUTRs = newUTRs[width(newUTRs) > 0]
  return(newUTRs)
}
#find and add cage max peaks as new tss's

###NB! Must have Gtf in global scope
getNewfivePrimeUTRs = function(fiveUTRs,dataName = standardCage){
  ###Read in cage files
  print(dataName)
  
  findNewTSS(fiveUTRs,dataName)
  
  newUTRs = addNewTssOnLeaders(fiveUTRs)
  return(newUTRs)
}
#needed to update lengths of utrs
findCageUTRFivelen = function(fiveUTRs,oldTxNames){
  newfiveprimeLen = sapply(fiveUTRs, function(x) sum(width(x)))
  newfiveprimeLen1 = newfiveprimeLen[match(oldTxNames,names(newfiveprimeLen))]
}

