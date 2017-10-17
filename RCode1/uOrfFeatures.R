getUOrfFeaturesMatrix = function(transcriptMatrix){
  transcriptNames = gsub("\\..*","", transcriptNames)
  grl = as.data.frame(rangesOfuORFs)
  
  #starts
  starts = aggregate(grl$start,by = list(grl$names),min)
  names(starts) = c("names","start")
  
  #ends
  ends = aggregate(grl$end,by = list(grl$names),min)
  names(ends) = c("names","end")
  
  #filtered name
  uorfName = getUORFnames(agUORFRibo$Group.1)
  
  #distance between uorf and cds
  distUC = distancebetweenUORFandCds(ends)
  
  inFrame = inFrameWithCDS(distUC)
  #overlaping cds
  overlapUOrfCds = getOverlappingCds(distUC)
  
  #Get rank order of uorf
  rank = getUOrfRankOrder(uorfName)
  
  pass_filter = getPassFilter(normUORFRNA,normUORFRFP)
  #TODO: ORF score, too slow!
  # orfScores = ORFScores(unlist(rangesOfuORFs[1:50]))
  
  
  matrixB = cbind(transcriptNames,uorfName,start = starts$start,end = ends$end,
                  width = ORFLengths$x,teUORF,
                  normUORFRNA, normUORFRFP,UCdists = distUC,frame = inFrame,
                  overlapUOrfCds,rank,pass_filter)
  
  colnames(matrixB)[6] <- "teUORF"
  colnames(matrixB)[9] <- "UCdist"
  
  matrixB = as.data.frame(as.matrix(matrixB))
  class(matrixB$transcriptNames) = "character"
  matrixB
}

#get distances between uorf end and start of that transcripts cds
distancebetweenUORFandCds = function(ends){
  ends$txNames = transcriptNames
  
  df.cds = as.data.frame(cds)
  #only first start exon
  df.cdsa = merge(aggregate( exon_rank ~ group_name, df.cds,min ),df.cds)
  
  #Merge to get correct order
  df.merge = merge(x = ends,y = df.cdsa, by.x = "txNames",by.y = "group_name",all.x = T)
  
  cdsStarts = df.merge[,c("txNames","names","end.x","start","strand")]

  #match distance by strand, different for positive and negative strand
  UCdist = lapply(1:nrow(cdsStarts),
                 function(i) ifelse(test = (cdsStarts[i,5] == "+"), yes = cdsStarts[i,4] - cdsStarts[i,3],no = cdsStarts[i,3] - cdsStarts[i,4]))
  cbind(UCdist)
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