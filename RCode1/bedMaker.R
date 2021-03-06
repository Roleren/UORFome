library(data.table)

#' bedName is bed name
#' chr, start,end,nameUORF,score = length? or ., strand
#' uscs bed 6 format
bed6 <- function(grl,bedName = NULL){
  if(class(grl) != "GRangesList") stop("grl, must be of class GRangesList")
  cat("writing bed file with name: ", bedName)
  
  dt.grl <- as.data.table(grl)
  
  bedColumns <- dt.grl[,c("seqnames","start","end","names",score = "width","strand")]
  
  if (is.null(bedName))
    write.table(x = bedColumns,file = paste0(uorfBedFolder,"rangesOfUorfsLight.bed") ,sep = "\t",col.names = F,row.names = F, quote = F)
  else
    write.table(x = bedColumns,file = bedName ,sep = "\t",col.names = F,row.names = F, quote = F)
}
#' grl A GRangesList
#' bedName is bed name
#' uscs bed 12 format
bed12 <- function(grl, bedName, fixChromoNaming = F){
  if(!ORFik:::is.grl(class(grl))) stop("grl, must be of class GRangesList")
  if (fixChromoNaming) print(seqlevels(grl))
  grl <- sortPerGroup(grl,ignore.strand = T) # <- sort bed way!
  
  dt.grl <- data.table(seqnames = ORFik:::seqnamesPerGroup(grl, F))
  dt.grl$start <- as.integer(ORFik:::firstStartPerGroup(grl,keep.names = F) -1)
  dt.grl$end <- ORFik:::lastExonEndPerGroup(grl,keep.names = F) #non inclusive end
  dt.grl$name <- names(grl)
  dt.grl$score <- ORFik:::widthPerGroup(grl, keep.names = F)
  dt.grl$strand <- ORFik:::strandPerGroup(grl, F)
  dt.grl$thickStart <- dt.grl$start
  dt.grl$thickEnd <- dt.grl$end
  dt.grl$rgb <- rep(0,length(grl))
  dt.grl$blockCount <- ORFik:::numExonsPerGroup(grl)
  blockSizes <- paste(width(grl), collapse = ",")
  names(blockSizes) <- NULL
  dt.grl$blockSizes <- blockSizes
  relativeStarts <- (start(grl) -1) - dt.grl$start
  blockStarts <- paste(relativeStarts, collapse = ",")
  names(blockStarts) <- NULL
  dt.grl$blockStarts <- blockStarts
  
  #chromStart + chromStarts[last] + blockSizes[last])
  #must equal chromEnd. 
  data.table::fwrite(x = dt.grl, file = bedName,
                     sep = "\t", col.names = F, row.names = F, quote = F)
  return(NULL)
}