
source("./features.R")


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
  
  transcriptNames =OrfToTranscriptNames(grlByORF)
  gen.starts = GroupGRangesFirstStart(grlByORF, F)
  gen.ends = GroupGRangesLastEnd(grlByORF,F)
    
  seqnames = GroupGRangesSeqnames(grlByORF, F)
  strands = GroupGRangesStrand(grlByORF, F)
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

