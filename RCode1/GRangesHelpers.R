

GroupGRangesByNames = function(gr){
  l = Rle(names(gr))
  t = unlist(lapply(1:length(l@lengths),function(x){ rep(x,l@lengths[x])}))
  grl = split(gr,t)
  names(grl) = unique(names(gr))
  return(grl)
}
GroupGRangesByOther = function(gr, other){
  if(length(gr) != length(other)){stop(" in GroupGRangesByOther: lengths of gr and other does not match")}
  l = Rle(other)
  t = unlist(lapply(1:length(l@lengths),function(x){ rep(x,l@lengths[x])}))
  grl = split(gr,t)
  names(grl) = unique(other)
  return(grl)
}

widthPerGRangesGroup = function(grl, keep.names = T){
  if(keep.names){
    return(sum(width(grl)))
  }else{
    return(as.integer(sum(width(grl))))
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

GroupGRangesLastEnd = function(grl,keep.names = T){
  
  if(keep.names){
    return(end(ptail(grlByORF,1L)))
  }else{
    return( as.integer(end(ptail(grlByORF,1L))))
  }
}
