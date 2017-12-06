

GroupGRangesByNames <- function(gr){
  l <- Rle(names(gr))
  t <- unlist(lapply(1:length(l@lengths),function(x){ rep(x,l@lengths[x])}))
  grl <- split(gr,t)
  names(grl) <- unique(names(gr))
  return(grl)
}

GroupGRangesByOther <- function(gr, other){
  if (length(gr) != length(other)){stop(" in GroupGRangesByOther: lengths of gr and other does not match")}
  l <- Rle(other)
  t <- unlist(lapply(1:length(l@lengths),function(x){ rep(x,l@lengths[x])}))
  grl <- split(gr,t)
  names(grl) <- unique(other)
  return(grl)
}

widthPerGRangesGroup = function(grl, keep.names = T){
  if(keep.names){
    return(sum(width(grl)))
  }else{
    return(as.integer(sum(width(grl))))
  }
}
GroupGRangesSeqnames <- function(grl, keep.names = T){
  if(keep.names){
    return(seqnames(phead(grl,1L)))
  }else{
    return(as.character(seqnames(phead(grl,1L))))
  }
}

GroupGRangesStrand <- function(grl, keep.names = T){
  if(keep.names){
    return(strand(phead(grl,1L)))
  }else{
    return(as.character(strand(phead(grl,1L))))
  }
}

GroupGRangesFirstExon = function(grl){
  return(phead(grl,1L))
}
GroupGRangesLastExon = function(grl){
  return(ptail(grl,1L))
}

GroupGRangesFirstStart = function(grl, keep.names = T){
  if(keep.names){
    return( start(GroupGRangesFirstExon(grl)))
  }else{
    return( as.integer(start(GroupGRangesFirstExon(grl))))
  }
}

GroupGRangesExonsPerGroup <- function(grl, keep.names = T){
  if (!is.null(names(unlist(grl, use.names = F)))){
    exonsPerGroup <- Rle(names(unlist(grl, use.names = F)))
    return(runLength(exonsPerGroup))
  } else if (!is.null(names(unlist(grl)))){
    exonsPerGroup <- Rle(names(unlist(grl)))
    return(runLength(exonsPerGroup))
  } else stop("no names to group exons")
}

GroupGRangesLastEnd = function(grl,keep.names = T){
  
  if(keep.names){
    return(end(ptail(grl,1L)))
  }else{
    return( as.integer(end(ptail(grl,1L))))
  }
}

GroupGRangesFixSeqnames <- function(grl){
  temp <- unlist(grl)
  seqnamesTransformed <- as.character(seqnames(temp))
  indexes <- which(nchar(seqnamesTransformed) < 6)
  temp <- temp[indexes]
  seqlevels(temp) <- sub(replacement = "chrY",pattern = "Y",seqlevels(temp))
  seqlevels(temp) <- sub(replacement = "chrX",pattern = "X",seqlevels(temp))
  seqlevels(temp) <- as.character(unique(seqnames(unlist(temp))))
  return(GroupGRangesByNames(temp))
}

GrangesSplitByExonSkeleton = function(grl, skeleton){
  
  unlNEW <- unlist(grl, use.names = F)
  if (!is.null(names(grl))){
    unlSkel <- unlist(skeleton[names(grl)], use.names = F)
  } else{
    unlSkel <- unlist(skeleton, use.names = F)
  }
  ol <- findOverlaps(query = unlNEW, subject = unlSkel, type = "any")
  
  if (is.null(names(unlSkel))){
    if (!is.null(names(grl))){
      unlSkel <- unlist(skeleton[names(grl)], use.names = T)
    } else{
      unlSkel <- unlist(skeleton, use.names = T)
    }
  } 
  if(is.null(names(unlNEW)))
    unlNEW <- unlist(grl, use.names = T)
  
  c <- ol[names(unlNEW[from(ol)]) == names(unlSkel[to(ol)])]
  N <- unlNEW[from(c)]
  ff <- unlSkel[to(c)]
  dups <- duplicated(from(c))
  start(N[dups]) <- start(ff[dups])
  dupsR <- duplicated(rev(from(c)))
  rN <- rev(N)
  end(rN[dupsR]) <- end(rev(ff)[dupsR])
  N <- rev(rN)
  
  l <- Rle(names(N))
  t <- unlist(lapply(1:nrun(l), function(x) {
    rep(x, runLength(l)[x])
  }))
  froms <- from(c)
  Inds <- rep(1, length(N))
  for (x in 2:length(N)) {
    if (t[x] != t[x - 1]) {
      Inds[x] <- 1
    }
    else {
      if (froms[x] != froms[x - 1]) {
        Inds[x] <- Inds[x - 1] + 1
      }
      else {
        Inds[x] <- Inds[x - 1]
      }
    }
  }
  N$names <- paste0(names(N), "_", Inds)
  newGRL <- split(N, t)
  names(newGRL) <- unique(names(N))
  return(newGRL)
}

toUniqueIDFromGR <- function(grl){
  seqnames = as.character(seqnames(phead(grl,1L)))
  strands = GroupGRangesStrand(grl,F)
  
  exonInfo <- paste(start(grl),width(grl))
  exonInfo = paste(exonInfo, sep = '', collapse = ';')
  names(exonInfo) <- NULL
  
  uorfID <- paste(seqnames, strands, exonInfo, sep = ",")
  return(uorfID)
}


toGRFromUniqueID <- function(uniqueIDs){
  if (ncol(uniqueIDs) != 1) stop("must only be 1 column")
  library(splitstackshape)
  library(stringr)
  splitList <- splitstackshape::cSplit_f(indt = uniqueIDs,
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
    t <- unlist(lapply(1:length(counts), function(x) {
      rep(x, counts[x])
    }))
    
    #IntegerLists
    starts <- c(t(starts))
    starts <- starts[!is.na(starts)]
    widths <- c(t(widths))
    widths <- widths[!is.na(widths)]
    
    gr = GRanges(seqnames = seqnamesUsed[t],
                 ranges = IRanges(start = starts, width = widths),
                 strand = strands[t])
    
    grl = GroupGRangesByOther(gr,t)
  } else {stop("not correct ncols")}
  
  return(grl)
}
