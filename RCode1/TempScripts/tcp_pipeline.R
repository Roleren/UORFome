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

#' heatmap
tcpHeatMap_int <- function(region, tx, df = getTCPdf(), outdir = p(mainFolder, "/tcp_plots/mir430/normal_"), 
                       shifting = NULL, scores = c("sum", "zscore"), upstream, downstream, 
                       zeroPosition = upstream, returnCoverage = FALSE) {
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

# TEMP VERSION
# tcpHeatMap_single <- function(region, tx, reads, outdir = p(mainFolder, "/tcp_plots/mir430/normal_"), 
#                               scores = "sum", upstream, downstream,  zeroPosition = upstream,
#                               returnCoverage = FALSE, logIt = FALSE, acLen = NULL, legendPos = "right") {
#   if (length(scores) != 1) stop("scores must exactly only 1 score type")
#   dt <- windowPerReadLength(region, tx, reads, upstream = upstream, downstream = downstream, 
#                             zeroPosition = zeroPosition, scoring = scores[1], acceptedLengths = acLen)
#   
#   if (logIt) dt$score <- log2(dt$score)
#   if (scores[1] == "log2sum") scores <- "log2(sum)"
#   if (scores[1] == "sum") scores <- "Count\nSum"
#   
#   #plot <- ORFik:::coverageHeatMap(coverage = dt, scoring = scores[1])
#   
#   dt$fraction <- factor(dt$fraction, levels = unique(dt$fraction), 
#                         labels = unique(dt$fraction))
#   
#   plot <- ggplot(dt, aes(x = position, y = fraction, 
#                          fill = score)) + geom_tile() + scale_fill_gradientn(colours = c("white", 
#                                                                                          "yellow2","yellow3", "lightblue", "blue", "navy", "navy", "black"), 
#                                                                              name = scores[1]) + xlab("Position relative to start of transcript") + 
#     ylab("Protected fragment length") + scale_x_continuous(breaks = seq.int(-upstream, downstream, by = 10)) + 
#     theme_bw(base_size=15) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_y_discrete(breaks = yAxisScaler(levels(dt$fraction))) + theme(legend.position = legendPos)
#   
#   ggsave(filename = outdir, plot = plot, width = 200, height = 100, units = "mm",
#          dpi = 300, limitsize = FALSE)  
#   if (returnCoverage) {
#     return(dt)
#   } else return(plot)
# }

tcpHeatMap_single <- function(region, tx, reads, outdir = p(mainFolder, "/tcp_plots/mir430/normal_"), 
                              scores = "sum", upstream, downstream,  zeroPosition = upstream,
                              returnCoverage = FALSE, logIt = FALSE, acLen = NULL, legendPos = "right", 
                              colors = NULL, addFracPlot = TRUE) {
  if (length(scores) != 1) stop("scores must exactly only 1 score type")
  dt <- windowPerReadLength(region, tx, reads, upstream = upstream, downstream = downstream, 
                            zeroPosition = zeroPosition, scoring = scores[1], acceptedLengths = acLen)
 
  if (logIt) dt$score <- log2(dt$score)
  if (scores[1] == "log2sum") scores <- "log2(sum)"
  
  dt$fraction <- factor(dt$fraction, levels = unique(dt$fraction), 
                              labels = unique(dt$fraction))
  if (is.null(colors)) colors <- c("white", "yellow2", "yellow3", "lightblue", "blue", "navy")
  if (colors == "high") colors <- c("white", "yellow2", "yellow3", "lightblue", "blue", "blue", "blue", "navy", "black")
  
  plot <- ggplot(dt, aes(x = position, y = fraction, 
                               fill = score)) + geom_tile() + 
    scale_fill_gradientn(colours = colors, name = scores[1]) + 
    xlab("Position relative to start of transcript") + 
    ylab("Protected fragment length") + scale_x_continuous(breaks = seq.int(-upstream, downstream, by = 10)) + 
    theme_bw(base_size = 15) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_y_discrete(breaks = yAxisScaler(levels(dt$fraction))) + theme(legend.position = legendPos)
  
  if (addFracPlot) {
    plot2 <- pSitePlot(dt, forHeatmap = TRUE)
    plot <- gridExtra::grid.arrange(plot2, plot + theme(legend.position = "bottom"),
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
                             outdir = p(mainFolder, "/tcp_plots/mir430/normal_"),
                             scores = c("sum", "zscore"), allTogether = FALSE) {
  
  varNames <- bamVarName(df)
  outputBams(df, leaders)
  coverage <- data.table()
  if (!allTogether) {
    libTypes <- libraryTypes(df)
    j <- 0
    for (i in 1:nrow(df)) { # For each stage
      j <- j + 1
      print(i)
      readsList <- list()
      for (lib in libTypes) { # For each library of that stage (SSU, LSU, RNA-seq, RIBO-seq)
        readsList <- list(readsList, get(varNames[j]))
      }
      readsList <- readsList[-1]
      transcriptWindowPer(leaders, cds, trailers, df[i,], outdir, scores, fractions,
                          readsList)
    }
  } else { # all combined
      coverage <- data.table()
      for (i in varNames) { # For each stage
        print(i)
        coverage <- rbindlist(list(coverage, ORFik:::splitIn3Tx(leaders, cds, trailers, 
                                                                get(i), fraction = i)))
      }
      for(s in scores) {
        a <- windowCoveragePlot(coverage, scoring = s)
        ggsave(paste0(outdir, "all_", s, ".png"), a, height = 10)
      }
  }
}

transcriptWindowPer <- function(leaders, cds, trailers, df = getTCPdf()[1,], 
                                outdir = p(mainFolder, "/tcp_plots/mir430/normal_"),
                                scores = c("sum", "zscore"),
                                reads, returnCoverage = FALSE) {
  
  libTypes <- libraryTypes(df)
  if (is(reads, "list") | is(reads, "GAlignmentsList") | is(reads, "GRangesList")) {
    if (length(libTypes) != length(reads)) stop("not matching length of reads and lib types in df!")
  } else if(!(is(reads, "GRanges") | is(reads, "GAlignments"))) {
    stop("reads must be GRanges or GAlignments")
  }
  coverage <- data.table()
  
  for(i in 1:length(reads)) {
    coverage <- rbindlist(list(coverage, ORFik:::splitIn3Tx(leaders, cds, trailers, 
                                                            unlist(reads[[i]]), fraction = libTypes[i])))
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
                          upstream = 75, downstream = 74, scores = "transcriptNormalized") {
  windows <- startRegion(region, tx, TRUE, upstream = upstream, downstream = downstream)
  hitMap <- metaWindow(reads, windows, zeroPosition = upstream, feature = "1", fraction = "123", 
                       scoring = scores, forceUniqueEven = F)
  hitMap[, score := score / sum(score)]
  if (scores == "transcriptNormalized") scores <- "Proportion of reads"
  plot <- ggplot(hitMap, aes(x = position, y = score)) + 
    geom_bar(stat = "identity", color = 'white', width=1, position = position_dodge(width=0.099))  + 
    xlab("Position relative to start of transcript") + ylab(scores) +
    scale_x_continuous(breaks = seq.int(-upstream, downstream, by = 10)) +  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    
  
  ggsave(filename = outdir, plot = plot, width = 350, height = 180, units = "mm",
         dpi = 300, limitsize = FALSE) 
  return(plot)
}

#' Experiment info table
getTCPMZ <- function(stage = c("64", "64", "shield", "shield"), 
                     type = c("MZ", "WT", "MZ", "WT"),
                     SSU = c("/export/valenfs/projects/adam/TCP_seq/valen_11/64cell_mzdicer/processed_29_11_18/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_SSU_reps_1_2_peaks_removed_translating_filter.bam",
                              "/export/valenfs/projects/adam/TCP_seq/valen_11/Shield_mzdicer/processed_29_11_18/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_SSU_reps_1_2_3_peaks_removed_translating_filter.bam"),
                     LSU = c("/export/valenfs/projects/adam/TCP_seq/valen_11/64cell_mzdicer/processed_29_11_18/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_LSU_reps_1_2_peaks_removed_translating_filter.bam",
                              "/export/valenfs/projects/adam/TCP_seq/valen_11/Shield_mzdicer/processed_29_11_18/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_LSU_reps_1_2_3_peaks_removed_translating_filter.bam")) {
  
  return(data.frame(SSU, LSU, stage, type, stringsAsFactors = FALSE))
}

getTCPdf <- function(stage = c("64", "sphere", "shield"), 
                     type = c("WT", "WT", "WT"),
                     SSU = c("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_SSU_reps_1_2_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/sphere_SSU_reps_1_2_3_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_SSU_reps_1_2_3_peaks_removed_translating_filter.bam"),
                     RFP = c("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/256Cell_trimmed.bam",
                             "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Dome_trimmed.bam",
                             "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Shield_trimmed.bam")) {
  
  return(data.frame(SSU, RFP, stage, type, stringsAsFactors = FALSE))
}

getTCPdfAll <- function(stage = c("64", "sphere", "shield"), 
                     type = c("WT", "WT", "WT"),
                     SSU = c("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_SSU_reps_1_2_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/sphere_SSU_reps_1_2_3_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_SSU_reps_1_2_3_peaks_removed_translating_filter.bam"),
                     RFP = c("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/256Cell_trimmed.bam",
                             "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Dome_trimmed.bam",
                             "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Shield_trimmed.bam"),
                     RNA = c("/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_2hpf_Tota_mRNA.bam",
                             "/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_4hpf_Total_mRNA.bam",
                             "/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_6hpf_Total_mRNA_merged.bam"),
                     LSU = c("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_LSU_reps_1_2_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/sphere_LSU_reps_1_2_3_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_LSU_reps_1_2_3_peaks_removed_translating_filter.bam")) {
  
  return(data.frame(SSU, RFP, RNA, LSU, stage, type, stringsAsFactors = FALSE))
}

getTCPNew <- function(){
  mergedF <- "/export/valenfs/projects/uORFome/withrRNA/aligned/"
  p <- paste0
  df2 <- data.frame(LSU = c(p(mergedF, "64_cell_LSU_V7.bam"), p(mergedF, "64_cell_LSU_V8.bam"),
                            p(mergedF, "shield_V5_merged_LSU.bam"), p(mergedF, "shield_V6_merged_LSU.bam"),
                            p(mergedF, "shield_V15_merged_LSU.bam"), p(mergedF, "shield_all_merged_LSU.bam"),
                            p(mergedF, "sphere1_V7_merged_LSU.bam"), p(mergedF, "sphere2_V7_merged_LSU.bam"),
                            p(mergedF, "sphere3_V7_merged_LSU.bam"), p(mergedF, "64_LSU_V12_4Ei.bam")),
                    stage = c(rep("64", 2), rep("shield", 4), rep("sphere", 3), "64"),
                    type = c(rep("WT", 2), rep("WT", 3), "WT_ALL", rep("WT", 3), "4Ei"), stringsAsFactors = FALSE)
  return(df2)
}



allBamFilesInFolder <- function(dir) {
  files <- grep(pattern = ".bam", x = list.files(dir, full.names = T), value = T)
  bai <- -grep(pattern = ".bai", x = files)
  if (length(bai)) {
    return(files[bai])
  } 
  return(files)
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

removeBadTxByRegion <- function(tx, reads, upstream, downstream, value = 200, extension = 100) {
  region <- ORFik:::startRegion(tx, extendLeaders(tx, extension = extension), upstream = upstream, downstream = downstream)
  counts <- countOverlaps(region, reads);summary(counts); sum(counts > value)
  return(leadersShield[!(counts > value)])
}

removeRepeatRegionTx <- function(tx, reads) {
  reads <- readGAlignments(bam);seqlevelsStyle(reads) <- seqlevelsStyle(tx);
  
  cov <- ORFik:::coveragePerTiling(tx, reads, TRUE, TRUE)
  means <- mean(cov)
  hits <- max(cov) > (200*means) & means > 1 
  return(tx[!hits])
}

#' Get peptides from singalp
#' @param outputSingalp
#' @param inputFa
#' @return character, a list of gene names
getPeptides <- function(outputSingalp = "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/signal_peptide_prediction/danio_rerio_peptides__summary.signalp5", 
                        inputFa = "/export/valenfs/projects/Håkon/ZF_RNA-seq_neuropep/signal_peptide_prediction/Danio_rerio.GRCz10.pep.all.fa") {
  faa = FaFile(inputFa)
  peps <- getSeq(faa)
  splits <- tstrsplit(names(peps), " ", fixed=TRUE, fill="<NA>")
  # proteins
  split <- splits[[1]]
  tab <- read.table(outputSingalp, sep = "\t")
  tab <- tab[tab$V2 != "OTHER",]$V1 
  valid <- split %in% tab
  # genes
  split <- splits[[4]]
  genes <- gsub(pattern = "gene:", x = split, replacement = "")
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
  
  by <- ifelse(pos > 70, ifelse(pos > 120, ifelse(pos > 300, 100, 50), 20), 10)
  if (max > 100) { 
    return(as.character(c(50, 100, max)))
  } else {
    return(as.character(seq.int(10, max, 10)))
  }
}
