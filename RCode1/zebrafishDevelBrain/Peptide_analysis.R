#' Get peptides from singalp or tmhmm(short form)
#' @param outputSingalp
#' @param inputFa
#' @return character, a list of gene names
getPeptides <- function(outputSingalp = "/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/signal_peptide_prediction/danio_rerio_peptides__summary.signalp5", 
                        inputFa = "/export/valenfs/projects/Hakon/ZF_RNA-seq_neuropep/signal_peptide_prediction/Danio_rerio.GRCz10.pep.all.fa") {
  fileExt <- tools::file_ext(inputFa)
  if(fileExt %in% c("fa", "fasta")) { # Reference
    faa = FaFile(inputFa)
    peps <- getSeq(faa)
    splits <- tstrsplit(names(peps), " ", fixed=TRUE, fill="<NA>")
    split <- splits[[1]]
    refgenes <- splits[[4]]
  } else stop("Bad input fasta!")
  
  # proteins 
  fileExt <- tools::file_ext(outputSingalp)
  if (fileExt == "signalp5") {
    tab <- read.table(outputSingalp, sep = "\t")
    tab <- tab[tab$V2 != "OTHER",]$V1 
  } else {
    hmm <- read.table(outputSingalp)
    colnames(hmm) <- c("pname", "len", "ExpAA", "first60", "PredHel", "Topology")
    hmm$PredHel <-  as.integer(sub(pattern = p("PredHel="), "", hmm$PredHel))
    tab <- hmm[hmm$PredHel == 0,]$pname # 0 transmemebrane helices
  }
  valid <- split %in% tab
  # genes
  genes <- gsub(pattern = "gene:", x = refgenes, replacement = "")
  genes <- genes[valid]
  return(genes)
}