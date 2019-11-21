library(ORFik)
library(ggplot2)


# Input of this are reads (bam files) after alignin process

# STEPS:
#1: Find most highly expressed isoform of each gene (based on RNA-seq coverage)  
#2: Use cage to update leaders
#3: Extend improperly trimmed regions (while trimmed reads exactly match transcripts)
#4: Remove high coverage peaks in transcripts. Reads with same start + stop coord >= 200X mean transcript coverage.
#5: Remove ambiguous reads from TCP-seq 40S complexes. Ambiguous reads are those with the same length as translating 80S complexes (25-35nt in length). 
#6: Do what you need to do

#' After aligning with STAR, run this
tcpPipeLine <- function() {
  # 1:
  tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
  names <- filterTranscripts(txdb)
  for(stage in stages) {
    
    tx[order(countOverlaps(tx, RNA), decreasing = TRUE)]
  }
  
  
  # 2:
  
  for(stage in stages) {
    
    # reassignTSSbyCage() #  or txdb ?
  }
  
  
  # 3: Trimming upwards, skip this now
  
  # 4: Remove all with extreme peaks (200x)
  # get coverage
  # find mean
  # find all with >= 200x, that has same start stop
  # remove those
  for(bam in files) {
    reads <- readGAlignments(bam);seqlevelsStyle(reads) <- seqlevelsStyle(tx);
    
    cov <- ORFik:::coveragePerTiling(tx, reads, TRUE, TRUE)
    means <- mean(cov)
    hits <- max(cov) > (200*means) & means > 1 
    # here are hits
    valid <- cov[hits]
    
    
    outdir <- p(mainFolder, "/")
    outNameSSU <- paste0(outdir, "SSU_", getRelativePathName(bam))
    outNameLSU <- paste0(outdir, "LSU_", getRelativePathName(bam))
    
    export(LSU, con = outNameLSU)
  }
  
  # coverage = something
  
  # 5: Split in SSU/LSU: (Ambigious read lengths in SSU)
  # the real split is done in the fractions, etc. 12, 13, 14, LSU: 18, 19 etc. 
  for(bam in files) {
    reads <- readGAlignments(bam);seqlevelsStyle(reads) <- "NCBI";
    
    hitsLSU <- ORFik:::readWidths(reads) >= 25 & ORFik:::readWidths(reads) <= 35 
    SSU <- reads[!hitsLSU] 
    LSU <- reads[hitsLSU]
    
    outdir <- p(mainFolder, "/")
    outNameSSU <- paste0(outdir, "SSU_", getRelativePathName(bam))
    outNameLSU <- paste0(outdir, "LSU_", getRelativePathName(bam))
    
    export(SSU, con = outNameSSU)
    export(LSU, con = outNameLSU)
  }
}

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

#' heatmap
heatMapL <- function(region, tx, df, outdir = p(mainFolder, "/tcp_plots/mir430/hm_"), 
                        scores = "sum", upstream, downstream,  zeroPosition = upstream,
                        logIt = FALSE, acLen = NULL, legendPos = "right", colors = NULL,
                        addFracPlot = TRUE, location = "TIS", shifting = NULL, skip.last = FALSE) {
  up <- upstream; down <- downstream
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, region)
    
    for (i in varNames) { # For each stage
      for (score in scores) {
        for (s in 1:length(shifting)) {
          if (s == 0) next
          shift <- shifting[s]
          if (length(upstream) > 1) {
            up <- upstream[s]
          }
          if (length(downstream) > 1) {
            down <- downstream[s]
          }
          format <- ".png"
          out <- ifelse(!is.null(shifting), paste0(outdir, df@experiment,"_hm_", location, "_",i , 
                                                   "_", shift, "_", score, format), 
                        paste0(outdir, df@experiment,"_hm_", location, "_", i,"_", score, format))
          
          tcpHeatMap_single(region, tx, reads = get(i), outdir = out, 
                            shifting = shift, scores = score, upstream = up, downstream = down, 
                            zeroPosition = zeroPosition, logIt = logIt, acLen = acLen, 
                            legendPos = legendPos, colors = colors, addFracPlot = addFracPlot, 
                            location = location, skip.last = skip.last)
        }
      }
    }
  }
}

tcpHeatMap_single <- function(region, tx, reads, outdir = p(mainFolder, "/tcp_plots/mir430/hm_"), 
                              scores = "sum", upstream, downstream,  zeroPosition = upstream,
                              returnCoverage = FALSE, logIt = FALSE, acLen = NULL, legendPos = "right", 
                              colors = NULL, addFracPlot = TRUE, location = "start site", shifting = NULL, 
                              skip.last = FALSE) {
  if (length(scores) != 1) stop("scores must exactly only 1 score type")
  if (!is.null(shifting)) {
    reads <- convertToOneBasedRanges(reads, method = shifting, addSizeColumn = TRUE)
  } 
  if (skip.last) {
    all_lengths <- sort(unique(readWidths(reads)))
    if (!is.null(acLen)) 
      all_lengths <- all_lengths[all_lengths %in% acLen]
    acLen <- all_lengths[-c((length(all_lengths)-0):length(all_lengths))]
  }
  dt <- windowPerReadLength(region, tx, reads, upstream = upstream, downstream = downstream, 
                            zeroPosition = zeroPosition, scoring = scores[1], acceptedLengths = acLen)
  
  if (logIt) dt$score <- log2(dt$score)
  if (scores[1] == "log2sum") scores <- "log2(sum)"
  
  dt$fraction <- factor(dt$fraction, levels = unique(dt$fraction), 
                              labels = unique(dt$fraction))
  if (is.null(colors))
    colors <- c("white", "yellow2", "yellow3", "lightblue", "blue", "navy")
  
  if (all(colors == "high")) 
    colors <- c("white", "yellow2", "yellow3", "lightblue", "blue", "blue", "blue", "navy", "black")
  
  
  plot <- ggplot(dt, aes(x = position, y = fraction, 
                               fill = score)) + geom_tile() + 
    scale_fill_gradientn(colours = colors, name = ORFik:::prettyScoring(scores[1])) + 
    xlab(paste("Position relative to", location)) + 
    ylab("Protected fragment length") + scale_x_continuous(breaks = c(seq.int(-upstream, downstream, by = 10), downstream)) + 
    theme_bw(base_size = 15) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_y_discrete(breaks = yAxisScaler(levels(dt$fraction))) + theme(legend.position = legendPos)
  
  if (addFracPlot) {
    plot2 <- pSitePlot(dt, forHeatmap = TRUE)
    plot <- gridExtra::grid.arrange(plot2, plot + theme(legend.position = "bottom", legend.key.width = unit(2, "lines")),
                         heights = c(1, 4))
  }
  
  ggsave(filename = outdir, plot = plot, width = 350, height = 180, units = "mm",
         dpi = 300, limitsize = FALSE)  
  if (returnCoverage) {
    return(dt)
  } else return(plot)
}

#' Make 100 bases size meta window for all libraries in input data.frame
transcriptWindow <- function(leaders, cds, trailers, df = getTCPdf(), 
                             outdir,
                             scores = c("sum", "zscore"), allTogether = FALSE, 
                             colors = rep("skyblue4", nrow(df)), 
                             windowSize = min(100, 
                                              min(widthPerGroup(leaders, F)),
                                              min(widthPerGroup(cds, F)),
                                              min(widthPerGroup(trailers, F)))
                             , returnPlot = FALSE) {
  if (windowSize != 100) message(paste0("NOTE: windowSize is not 100!
                                        It is ", windowSize))
  
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, leaders)
    coverage <- data.table()
    if (!allTogether) {
      stop("fix!")
      libTypes <- libraryTypes(df)
      j <- 0
      for (i in 1:nrow(df)) { # For each stage
        j <- j + 1
        print(i)
        readsList <- list()
        for (lib in libTypes) {
          # For each library of that stage (SSU, LSU, RNA-seq, RIBO-seq)
          readsList <- list(readsList, get(varNames[j]))
        }
        readsList <- readsList[-1]
        transcriptWindowPer(leaders, cds, trailers, df[i,], outdir, scores,
                            fractions, readsList)
      }
    } else { # all combined
        coverage <- data.table()
        for (i in varNames) { # For each stage
          print(i)
          coverage <- rbindlist(list(coverage,
                                     splitIn3Tx(leaders, cds, trailers, 
                                                get(i), fraction = i,
                                                windowSize = windowSize)))
        }
        for(s in scores) {
          a <- windowCoveragePlot(coverage, scoring = s, colors = colors)
          ggsave(pasteDir(outdir, paste0(df@experiment,"_cp_all_", s, ".png"))
                 , a, height = 10)
        }
    }
  }
  if (returnPlot) return(a)
}

transcriptWindowPer <- function(leaders, cds, trailers, df, 
                                outdir, scores = c("sum", "zscore"),
                                reads, returnCoverage = FALSE,
                                windowSize = 100) {
  
  libTypes <- libraryTypes(df)
  if (is(reads, "list") | is(reads, "GAlignmentsList") |
      is(reads, "GRangesList")) {
    if (length(libTypes) != length(reads)) 
      stop("not matching length of reads and lib types in df!")
  } else if(!(is(reads, "GRanges") | is(reads, "GAlignments"))) {
    stop("reads must be GRanges or GAlignments")
  }
  coverage <- data.table()
  
  for(i in 1:length(reads)) {
    coverage <- rbindlist(list(coverage, 
                               splitIn3Tx(leaders, cds, trailers, 
                                          unlist(reads[[i]]), 
                                          fraction = libTypes[i],
                                          windowSize = windowSize)))
  }
  
  return(plotHelper(coverage, df, outdir, scores, returnCoverage))
}

#' TIS window 150 bases
regionWindow <- function(region, tx, df = getTCPdf(), 
                         outdir = p(mainFolder, "/tcp_plots/TIS_region/target_"),
                         upstream = 75, downstream = 74, scores = c("sum", "zscore"),
                         title = "cds metacoverage") {
  # find most translated genes
  windowsOne <- startRegion(region, tx, TRUE, upstream = upstream, downstream = downstream)
  
  libTypes <- libraryTypes(df)
  coverage <- data.table()
  for (i in 1:nrow(df)) { # For each stage
    print(i)
    stage <- df$stage[i]
    type <- df$type[i] # Add type to fraction
    
    for (lib in libTypes) { # For each library of that stage (SSU, LSU, RNA-seq, RIBO-seq)
      reads <- readGAlignments(df[i, lib]) ;seqlevelsStyle(reads) <- seqlevelsStyle(windowsOne)[1]
      hitMap <- metaWindow(reads, windowsOne, scoring = "meanPos", withFrames = FALSE, returnAs = "data.table", 
                              fraction = lib, feature = stage)
      
      coverage <- rbindlist(list(coverage, hitMap))
    }
  }
  
  return(plotHelper(coverage, df, outdir, scores, returnCoverage = FALSE, title = title))
}

regionWindowAll <- function(one, two, tx, df = getTCPdf(), 
                             outdir = p(mainFolder, "/tcp_plots/TIS_region/target_"),
                             upstream = 75, downstream = 74, scores = c("sum", "zscore"),
                             seperations = c("lowly", "highly"), title = "cds metacoverage") {
  # find most translated genes
  windowsOne <- startRegion(one, tx, TRUE, upstream = upstream, downstream = downstream)
  windowTwo <- startRegion(two, tx, TRUE, upstream = upstream, downstream = downstream)
  
  libTypes <- libraryTypes(df)
  coverage <- data.table()
  for (i in 1:nrow(df)) { # For each stage
    print(i)
    stage <- df$stage[i]
    type <- df$type[i] # Add type to fraction
    
    for (lib in libTypes) { # For each library of that stage (SSU, LSU, RNA-seq, RIBO-seq)
      reads <- readGAlignments(df[i, lib]) ;seqlevelsStyle(reads) <- seqlevelsStyle(windowsOne)[1]
      
      hitMapOne <- metaWindow(reads, windowsOne, scoring = "meanPos", withFrames = FALSE,
                              fraction = paste0(seperations[1], "_", lib), feature = stage)
      hitMapTwo <- metaWindow(reads, windowTwo, scoring = "meanPos", withFrames = FALSE, 
                              fraction = paste0(seperations[2], "_", lib), feature = stage)
      
      coverage <- rbindlist(list(coverage, hitMapOne, hitMapTwo))
    }
  }
  return(plotHelper(coverage, df, outdir, scores, returnCoverage = FALSE, title = title, 
                    colors =  c("skyblue4", "orange")[c(1,1,2,2)]))
}

regionBarPlot <- function(region, tx, reads, 
                          outdir = p(mainFolder, "/tcp_plots/TIS_region/target_"),
                          upstream = 75, downstream = 74, scores = "sum", proportion = TRUE) {
  windows <- startRegion(region, tx, TRUE, upstream = upstream, downstream = downstream)
  hitMap <- metaWindow(reads, windows, zeroPosition = upstream, feature = "1", fraction = "123", 
                       scoring = scores, forceUniqueEven = F)
  if (proportion) {
    hitMap[, score := score / sum(score)]
    scores <- "Proportion of reads"
  }
  plot <- ggplot(hitMap, aes(x = position, y = score)) + 
    geom_bar(stat = "identity", color = 'white', width=1, position = position_dodge(width=0.099))  + 
    xlab("Position relative to start of transcript") + ylab(scores) +
    scale_x_continuous(breaks = seq.int(-upstream, downstream, by = 10)) +  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    
  if (!is.null(outdir)) {
    ggsave(filename = outdir, plot = plot, width = 350, height = 180, units = "mm",
           dpi = 300, limitsize = FALSE) 
  }
  
  return(plot)
}

#' Find ratio of iMet to met amino acids
#' Also get basic tRNA statistics
#' @param fa location of fasta file with tRNA aligned to, or a biostring object of preloaded
#' @param bamFiles paths to bamFiles matched to fa. A character vector of path to bam files aligned to fa.
tRNA_stats <- function(fa, bamFiles) {
  if (is(fa, "FaFile")) {
    getFasta("/export/valenfs/data/references/Zv10_zebrafish/tRNAscan-SE-2.0/GRcZ10_tRNA_scan_output.fa")
    trna <- getSeq(fa)
  } else trna <- fa
  
  
  met <- trna[grep(" Met ", names(trna), value = FALSE, ignore.case = TRUE)]
  imet <- trna[grep(" iMet ", names(trna), value = FALSE, ignore.case = TRUE)]
  
  seqsMet <- gsub(pattern = " .*", replacement = "", x = names(met))
  seqsiMet <- gsub(pattern = " .*", replacement = "", x = names(imet))
  df <- data.frame(file = bamFiles, total_tRNA_reads = NA, met = NA, imet = NA, imet_met_ratio = NA,imet_total_percent = NA)
  
  for(i in seq(length(bamFiles))) {
    ga <- ORFik:::readBam(bamFiles[i])
    seqsGA <- as.character(seqnames(ga))
    
    df[i,]$total_tRNA_reads <- length(ga)
    df[i,]$met <- sum(!is.na(chmatch(seqsGA, seqsMet)))
    df[i,]$imet <- sum(!is.na(chmatch(seqsGA, seqsiMet)))
    df[i,]$imet_met_ratio <-  df[i,]$imet / df[i,]$met
    df[i,]$imet_total_percent <-  (df[i,]$imet / df[i,]$total_tRNA_reads)*100
  }
  
  return(df)
}

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

removeBadTxByRegion <- function(tx, reads, upstream, downstream, value = 5, extension = 100) {
  region <- ORFik:::startRegion(tx, extendLeaders(tx, extension = extension), upstream = upstream, downstream = downstream)
  coverage <- coveragePerTiling(region, reads, T, as.data.table = T)
  #counts <- countOverlaps(region, reads);print(summary(counts)); print(sum(counts > value)); print(which(counts > value))
  coverage[, median_per_gene := median(count), by = genes]
  toFilter <- coverage[count > median_per_gene*value + 25]
  print(summary(coverage$count))
  print(toFilter)
  print(names(tx[unique(toFilter$genes)]))
  return(tx[-unique(toFilter$genes)])
}

removeRepeatRegionTx <- function(tx, reads) {
  if (is.character(reads) | is.factor(reads)) 
    reads <- ORFik:::readBam(reads, tx)
  
  cov <- ORFik:::coveragePerTiling(tx, reads, TRUE, TRUE)
  means <- mean(cov)
  hits <- max(cov) > (200*means) & means > 1 
  return(tx[!hits])
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

plotHelper <- function(coverage, df, outdir, scores, returnCoverage = FALSE,
                       title = "coverage metaplot", colors = c("skyblue4", "orange"),
                       plotFunction = "windowCoveragePlot") {
  if (!is.null(outdir)) {
    for(s in scores) {
      stage <- df$stage[1]
      type <- df$type[1]
      sample_name <- paste0(stage, "_", type, "_", s) # What happens on cds ?
      outName <- paste0(outdir, sample_name, ".png")
      if (plotFunction == "windowCoveragePlot") {
        windowCoveragePlot(coverage, output = outName, scoring = s, title = title, colors = colors)
      } else if (plotFunction == "coverageHeatMap") {
        plot <- ORFik:::coverageHeatMap(coverage = coverage, scoring = s)
        ggsave(outName, plot = plot, width = 350, height = 180, units = "mm",
               dpi = 300, limitsize = FALSE)
      } else stop("invalid plot name")
        
    }
  }
  if (returnCoverage) return(coverage)
  return(NULL)
}

yAxisScaler <- function(covPos) {
  covPos <- as.integer(covPos)
  pos <- length(covPos)
  min <- min(covPos)
  max <- max(covPos)
  
  by <- ifelse(pos > 40, ifelse(pos > 70, ifelse(pos > 120, ifelse(pos > 300, 100, 50), 20), 10), 1)
  if (max > 100) { 
    return(as.character(c(50, 100, max)))
  } else {
    if (by < 10) return(as.character(seq.int(10, max, by)))
    return(as.character(seq.int(10, max, 10)))
  }
}
