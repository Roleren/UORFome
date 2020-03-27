
toUniqueIDFromGR <- function(grl, with.tx = FALSE){
  seqnames = as.character(seqnames(phead(grl,1L)))
  strands = ORFik:::strandPerGroup(grl,F)
  
  exonInfo <- paste(start(grl),width(grl))
  exonInfo = paste(exonInfo, sep = '', collapse = ';')
  names(exonInfo) <- NULL
  
  uorfID <- paste(seqnames, strands, exonInfo, sep = ",")
  if (with.tx) {
    uorfID <- paste(uorfID, OrfToTxNames(grl))
  }
  return(uorfID)
}


toGRFromUniqueID <- function(uniqueIDs){
 
  if (class(uniqueIDs)[1] != "data.table"){
    if (class(uniqueIDs)[1] == "character"){
      uniqueIDs <- data.table(uorfID = uniqueIDs)
      
    } else stop("uniqueIDs must either be data.table or character!")
  }else if (ncol(uniqueIDs) != 1) stop("must only be 1 column")
  colnames(uniqueIDs) = "uorfID"
  
  library(splitstackshape)
  library(stringr)
  # not working anymore, update to cSplit
  # splitList <- splitstackshape::cSplit_f(indt = uniqueIDs,
  #                         stripWhite = F, sep = ",", splitCols = "uorfID")
  splitList <- splitstackshape::cSplit(indt = uniqueIDs,
                                         stripWhite = F, sep = ",", splitCols = "uorfID")
  
  if (ncol(splitList) == 3){
    a <- splitstackshape::cSplit(splitList,splitCols = "uorfID_3", 
                                 sep = ";",drop = T)
    rm(splitList)
    seqnamesUsed <- as.character(unlist(a[,1]))
    strands <- as.character(unlist(a[,2]))
    a <- a[,-c(1,2)]
    
    # Split on exons
    exonsList <- as.data.table(matrix(ncol = ncol(a)*2, nrow = nrow(a)))
    j = 1
    for(i in 1:ncol(a)){
      
      exons <- splitstackshape::cSplit(splitCols = 1,indt = a[,i,with=F],  sep = " ", drop = T)
      exonsList[,j] <- exons[,1]
      exonsList[,ncol(a)+j] <- exons[,2]
      j <- j + 1
    }
    
    starts <- exonsList[, 1:ncol(a)]
    widths <- exonsList[, (ncol(a)+1):(ncol(a)*2)]
    rm(exonsList)
    rm(exons)
    rm(a)
    counts <- rowSums(!is.na(starts))
    #counts <-  unlist(lapply(1:nrow(starts), function(x) sum(!is.na(starts[x]))))
    t <- unlist(lapply(seq.int(length(counts)), function(x) {
      rep.int(x, counts[x])
    }))
    
    #IntegerLists
    starts <- c(t(starts))
    starts <- starts[!is.na(starts)]
    widths <- c(t(widths))
    widths <- widths[!is.na(widths)]
    
    gr = GRanges(seqnames = seqnamesUsed[t],
                 ranges = IRanges(start = starts, width = widths),
                 strand = strands[t])
    
    grl = groupGRangesBy(gr,t)
    grl <- sortPerGroup(grl)
  } else stop("not correct ncols in uniqueIDs")
  #test the orfs
  widths <- ORFik:::widthPerGroup(grl,F)
  if(sum(widths %% 3) != 0) stop("widths of uorfs are not %3 = 0!")
  if (length(grl) != nrow(uniqueIDs)) stop("not all ranges was reconstructed properly, check data!")
  return(grl)
}

getNegExonGrl <- function(grl) {
  grl <- grl[!strandBool(grl)]
  grl <- grl[ORFik:::numExonsPerGroup(grl) > 1]
  return(grl)
}
