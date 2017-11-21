getUOrfFeaturesMatrix = function(transcriptMatrix){
  #stop("not working now!!!")
  
  grORFs = unlist(rangesOfuORFs, use.names = F)
  grlByORF = GroupGRangesByOther(grORFs, grORFs$names)
  txCoords = getUORFTranscriptCoordinates()
  
  #starts, remake these to faster versions
  # starts = aggregate(grl$start,by = list(grl$names),min)
  starts = GroupGRangesFirstExon(txCoords)
  
  ends = GroupGRangesLastExon(txCoords)
  
  #distance between uorf and cds
  distUC = distancebetweenUORFandCds(ends)
  
  inFrame = inFrameWithCDS(distUC)
  #overlaping cds
  overlapUOrfCds = getOverlappingCds(distUC)
  
  #filtered name
  uorfName = unique(getUORFnames(grORFs$names))
  
  #Get rank order of uorf
  rank = getUOrfRankOrder(uorfName)
  
  pass_filter = getPassFilter(normUORFRNA,normUORFRFP)
  
  tx.widths = widthPerGRangesGroup(txCoords)
  
  transcriptNames =gsub("_[0-9]*","", names(grlByORF)) 
  gen.starts = GroupGRangesFirstStart(grlByORF, F)
  gen.ends = GroupGRangesLastEnd(grlByORF,F)
    
  seqnames = as.character(seqnames(phead(grlByORF,1L)))
  strands = as.character(strand(phead(grlByORF,1L)))
  #TODO: ORF score, too slow!
  # orfScores = ORFScores(unlist(rangesOfuORFs[1:50]))
  
  uorfID = paste(transcriptNames,seqnames, gen.starts,gen.ends,strands)
  tissue = rep(tissueUsed,length(gen.starts)) #this is scary I think!
  matrixB = cbind(transcriptNames,uorfName,seqnames,start = gen.starts,end = gen.ends,
                  width = tx.widths,teUORF,
                  normUORFRNA, normUORFRFP,UCdists = distUC,frame = inFrame,
                  overlapUOrfCds,rank,pass_filter,uorfID,tissue)
  
  colnames(matrixB)[7] <- "teUORF"
  colnames(matrixB)[10] <- "UCdist"
  
  matrixB = as.data.table(as.matrix(matrixB))
  class(matrixB$transcriptNames) = "character"
  return(matrixB)
}

getUORFTranscriptCoordinates = function(){
  n = unlist(rangesOfuORFs, use.names = F)
  #names of fives go to seqnames, so need to be removed
  txCoordUORF = mapToTranscripts(n,fiveUTRs)
  txCoordUORF = txCoordUORF[names(n[txCoordUORF$xHits]) == seqnames(txCoordUORF)]
  txCoordUORF$xHits = NULL;txCoordUORF$transcriptsHits = NULL
  #sort them by uorfs
  txCoordUORF$names = n$names
  return(GroupGRangesByOther(txCoordUORF, txCoordUORF$names))
}

#get distances between uorf end and start of that transcripts cds
distancebetweenUORFandCds = function(ends){
  cdsFirstExons = phead(cds,1L)
  namesToUse = as.character( unlist(unique( seqnames(ends))))
  cdsToUse = cdsFirstExons[namesToUse]
  c = unlist(cdsToUse, use.names = F)
  cdsToTranscript = mapToTranscripts(c,fiveUTRs)
  c = unlist(cdsToUse)
  cdsToUse = cdsToTranscript[names(c[cdsToTranscript$xHits]) == seqnames(cdsToTranscript)]
  
  endsPos = ends[as.character(strand(ends)) == "+"]
  endsMin = ends[as.character(strand(ends)) == "-"]
  
  distPos = start(cdsToUse[as.character(strand(cdsToUse)) == "+"]) - as.integer(end(endsPos))
  distMin = as.integer(start(endsMin)) - end(cdsToUse[as.character(strand(cdsToUse)) == "-"])
  dists = rep(NA,length(names(unlist(ends, use.names = F))))
  dists[as.character(strand(unlist(ends))) == "+"] = distPos  
  dists[as.character(strand(unlist(ends))) == "-"] = distMin
  return(dists)
}

#check if reading frame changes. %3 = 0
inFrameWithCDS = function(distUC){
  #for each distance found beween uorfEnd and cds start, do %3
  frame = cbind(unlist(distUC) %% 3)
  colnames(frame) = "frame"
  as.integer(frame)
}
#check if uorf is overlapping cds
getOverlappingCds = function(distUC){
  overlap = cbind(distUC < 0)
  colnames(overlap) = "overlap"
  overlap
}
#Get the uorf number in transcript from left, ig. second uorf -> 2
getUOrfRankOrder = function(uorfName){
  n = cbind(as.integer(gsub(".*_","", uorfName)))
  colnames(n) = "rank"
  n
}

getTissue = function(){
  name = NA
  if(exists("tissue")){
    name = tissue
  }else if(exists("cageName")){
    
  }
  cbind(rep(1:length(transcriptName),name))
}

getPassFilter = function(normUORFRNA,normUORFRFP){
  pass_filter = normUORFRNA > 0.1 & normUORFRFP > 0.1
  pass_filter[is.na(pass_filter)] = F
  pass_filter
}

getUORFnames = function(unfilteredNames){
  gsub(".*\\.","", unfilteredNames)
}

getKozacSequenceScore = function(grl, fastaSeq){
  #reassign start of + strands, and restrict end
  #reassign end of - strands, and restrict start
  
  #get all sequences
  
  #score
  
}

ORFScores = function(ORFs = NULL){
  
  
  #ORFscore
  tilex <- tile(ORFs, width=1)
  
  tilex1 <- lapply(tilex, function(x){
    if(as.vector(strand(x) == "+")[1]){
      x[seq(1, length(x), 3)]
    }else{
      x[seq(length(x), 1, -3)]
    }})
  tilex2 <- lapply(tilex, function(x){
    if(as.vector(strand(x) == "+")[1]){
      x[seq(2, length(x), 3)]
    }else{
      x[seq(length(x)-1, 1, -3)]
    }})
  tilex3 <- lapply(tilex, function(x){
    if(as.vector(strand(x) == "+")[1]){
      x[seq(3, length(x), 3)]
    }else{
      x[seq(length(x)-2, 1, -3)]
    }})
  
  countsTile1 <- countOverlaps(GRangesList(tilex1), RFP)#my edits!!
  countsTile2 <- countOverlaps(GRangesList(tilex2), RFP)#my edits!!
  countsTile3 <- countOverlaps(GRangesList(tilex3), RFP)#my edits!!
  
  RP = countsTile1 + countsTile2 + countsTile3
  
  Ftotal <- RP/3
  
  tile1 <- (countsTile1 - Ftotal)^2 / Ftotal
  tile2 <- (countsTile2 - Ftotal)^2 / Ftotal
  tile3 <- (countsTile3 - Ftotal)^2 / Ftotal
  
  dfORFs <- NULL
  dfORFs$frame_zero_RP <- countsTile1
  dfORFs$frame_one_RP <- countsTile2
  dfORFs$frame_two_RP <- countsTile3
  
  ORFscore <- log2(tile1 + tile2 + tile3 + 1)
  revORFscore <-  which(tile1 < tile2 | tile1 < tile3)
  ORFscore[revORFscore] <- -1 * ORFscore[revORFscore]
  ORFscore[is.na(ORFscore)] <- 0
  dfORFs$ORFscore <- ORFscore
}