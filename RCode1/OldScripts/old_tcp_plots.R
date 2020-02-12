#' heatmap
tcpHeatMap <- function(txdb, df = getTCPdf(), outdir = p(mainFolder, "/tcp_plots/mir430/normal_"), 
                       shifting = NULL, scores = c("sum", "zscore"), zeroPosition = 100) {
  outdir <- paste0(outdir, shifting, "_")
  
  names2 <- filterTranscripts(txdb, 100, 50, 0)
  cds <- cdsBy(txdb, use.names = TRUE)[names2]
  tx <- exonsBy(txdb, use.names = TRUE)[names2]
  
  tcpHeatMap_int(cds, tx, df = df, outdir = outdir, 
                 shifting = shifting, scores = scores, upstream = 100, downstream = 49, 
                 zeroPosition = zeroPosition)
}
#TODO FIX!
#' heatmap
tcpHeatMap_int <- function(region, tx, df = getTCPdf(), outdir = p(mainFolder, "/tcp_plots/mir430/normal_"), 
                           shifting = NULL, scores = c("sum", "zscore"), upstream, downstream, 
                           zeroPosition = upstream, returnCoverage = FALSE) {
  stop("Not fixed!")
  window <- ORFik:::startRegion(region, tx, TRUE, upstream, downstream)
  
  outdir <- paste0(outdir, shifting, "_")
  libTypes <- libraryTypes(df)
  
  for (i in 1:nrow(df)) { # For each stage
    print(i)
    for (lib in libTypes) { # For each library of that stage (SSU, LSU, RNA-seq, RIBO-seq)
      reads <- readGAlignments(df[i, lib]) ;seqlevelsStyle(reads) <- seqlevelsStyle(window)[1]
      if (!is.null(shifting)) {
        reads <- ORFik:::convertToOneBasedRanges(reads, method = shifting, addSizeColumn = TRUE)
      } 
      
      
      all_lengths <- sort(unique(ORFik:::readWidths(reads)))
      dt <- data.table()
      for(l in all_lengths) {
        dt <- data.table::rbindlist(list(dt, ORFik:::metaWindow(x = reads[ORFik:::readWidths(reads) == l], windows = window,
                                                                zeroPosition = zeroPosition,
                                                                forceUniqueEven = T, fraction = l)))
      }
      if (!is.null(dt$count)) {
        dt$score <- dt$count
        dt$count <- NULL
      }
      dt$score <- log2(dt$score)
      
      cov <- plotHelper(dt, df[i,], paste0(outdir, "_", lib), scores, plotFunction = "coverageHeatMap", 
                        returnCoverage)
    }
  }
  return(cov)
}

# Fastp, peptide analysis
if (0) {
  getFasta("/export/valenfs/data/references/Zv10_zebrafish/tRNAscan-SE-2.0/GRcZ10_tRNA_scan_output.fa")
  trna <- getSeq(fa)
  # Test
  bamFiles <- allBamFilesInFolder("/export/valenfs/data/processed_data/TCP-seq/valen_2018_zebrafish_6_trim3_15nt_tRNAscan/tRNAscan_match_GRCz10")
  df6 <- tRNA_helper(trna, bamFiles)
  # Result, equal to Adams
  # New ones
  
  bamFiles <- allBamFilesInFolder("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_15_trim3_15nt_tRNAscan/tRNA_match_GRCz10")
  df15 <- tRNA_helper(trna, bamFiles)
  bamFiles <- allBamFilesInFolder("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_14_trim3_15nt_tRNAscan/tRNA_match_GRCz10")
  df14 <- tRNA_helper(trna, bamFiles)
  bamFiles <- allBamFilesInFolder("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_13_trim3_15nt_tRNAscan/tRNA_match_GRCz10")
  df13 <- tRNA_helper(trna, bamFiles)
  bamFiles <- allBamFilesInFolder("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_12_trim3_15nt_tRNAscan/tRNA_match_GRCz10")
  df12 <- tRNA_helper(trna, bamFiles)
  bamFiles <- allBamFilesInFolder("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_11_trim3_15nt_tRNAscan/tRNA_match_GRCz10")
  bamFiles <- bamFiles[-(2:4)]
  df11 <- tRNA_helper(trna, bamFiles)
  bamFiles <- allBamFilesInFolder("/export/valenfs/data/processed_data/TCP-seq/valen_2018_zebrafish_10_trim3_15nt_tRNAscan/tRNA_match_GRCz10")
  df10 <- tRNA_helper(trna, bamFiles)
}

#' Get peptides from singalp or tmhmm(short form)
#' @param outputSingalp
#' @param inputFa
#' @return character, a list of gene names
getPeptides <- function(outputSingalp = "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/signal_peptide_prediction/danio_rerio_peptides__summary.signalp5", 
                        inputFa = "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/signal_peptide_prediction/Danio_rerio.GRCz10.pep.all.fa") {
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