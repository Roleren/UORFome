

#####################################GETTESS##########################################
#seq: sequences of either: 5', 3', cds or uorfs
#dm1: detection method 1, RNA
#dm2: detection method 2, RFP
#seqLength: length of seq
#al: length of all sequences, might be NULL
#specificTE: what type you  want to get TE of, CDS, LEADER, 3', uorfs etc.
RiboFPKM = function(seq, riboFilePath){
  if (class(riboFilePath) != "character" 
      && class(riboFilePath) != "GAlignments"
      && class(riboFilePath) != "GAlignmentPairs" ){
    stop("riboFilePath must be either character, Galignments or GalignmentPairs")
  }
  if (class(riboFilePath) != "GAlignments"){
    riboGAllignment <- readGAlignments(riboFilePath)
  }
  
  libraryRPF <- length(riboGAllignment)
  
  # should check here that grl grouping is correct
  #gr <- unlist(seq, use.names = F)
  #grl <- GroupGRangesByOther(gr, gr$names)
  
  overlapRFP <- countOverlaps(seq, riboGAllignment)
  
  # Find lengths used for uorf fpkm values
  ORFLengths <- widthPerGRangesGroup(seq, F)
  
  fpkm <- fpkm(overlapRFP, ORFLengths, libraryRPF)# normalize by orf
  return(fpkm)
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

getORFnames = function(unfilteredNames){
  gsub(".*\\.","", unfilteredNames)
}

getAllFeatures <- function(){
  getCDS()
  getThreeUTRs()
  RFP <- readGAlignments("/export/valenfs/data/processed_data/Ribo-seq/gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/Gonzalez_C_2014.Human.brain.RPF.GRCh38.SRR1562540.bam")
  RNA <- readGAlignments("/export/valenfs/data/processed_data/RNA-seq/gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/Gonzalez_C_2014.Human.brain.RNA.GRCh38.SRR1562546.bam")
  
  dt <- ORFik:::allFeatures(grl = grl,orfFeatures = T, RFP = RFP, RNA = RNA,
                             Gtf = Gtf, fiveUTRs = fiveUTRs, cds = cds,
                            threeUTRs = threeUTRs, faFile = fa, riboStart = 26, riboStop = 34,
                            extension = 1000)
}




