
#Since cage redefines five utr lengths, the total lengths of tx must be updated
findCageUTRFivelen = function(fiveUTRs, oldTxNames){
  newfiveprimeLen <- widthPerGRangesGroup(fiveUTRs)
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
  fiveUTRs = groupGRangesBy(unlistfgr)
  
  return(fiveUTRs)
}

