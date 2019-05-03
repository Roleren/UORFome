# Input of this are reads (bam files) after alignin process

# STEPS:
#1: Find most highly expressed isoform of each gene (based on RNA-seq coverage)  
#2: Use cage to update leaders
#3: Extend improperly trimmed regions (while trimmed reads exactly match transcripts)
#4: Remove high coverage peaks in transcripts. Reads with same start + stop coord >= 200X mean transcript coverage.
#5: Remove ambiguous reads from TCP-seq 40S complexes. Ambiguous reads are those with the same length as translating 80S complexes (25-35nt in length). 
#6: 

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
    reads <- readGAlignments(bam);seqlevelsStyle(reads) <- "NCBI";
    
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
  
  # 5: Split in SSU/LSU:
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



# Standard
# SSU = "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/combined_time_points/all_three_SSU_peaks_removed_200_removed_translating_lengths_25_35.bam"
# LSU = "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/combined_time_points/all_three_LSU_peaks_removed_200_selected_translating_lengths_25_35.bam"
# stage <- "64"
# type <- "WT"


# score <- c("sum", "zscore")
# 
# gtfPath <- "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81.gtf"
# txdb <- makeTxDbFromGFF(file = gtfPath)
# seqlevelsStyle(txdb) <- "NCBI"
# 
# leaders = fiveUTRsByTranscript(txdb,use.names = T)
# cds <- cdsBy(txdb,"tx", use.names = TRUE)
# trailers = threeUTRsByTranscript(txdb, use.names = TRUE)
# tx <- exonsBy(txdb, by = "tx", use.names = TRUE)

# names <- names(leaders)[(names(leaders) %in% names(cds)) & (names(leaders) %in% names(trailers))]
# leaders <- leaders[names]
# cds <- cds[names]
# trailers <- trailers[names]




#' heatmap
tcpHeatMap <- function(txdb, df = getTCPdf()) {
  names2 <- filterTranscripts(txdb, 100, 50, 0)
  cds <- cds[names2]
  window <- ORFik:::startRegion(cds,tx, TRUE, 100, 49)
  
  for(i in 1:nrow(df)) {
    print(i)
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    #readsSSU <- GRanges(readGAlignments(SSU)) ;seqlevelsStyle(readsSSU) <- "NCBI";
    readsLSUbam <- readGAlignments(LSU)
    readsLSU <- GRanges(readGAlignments(LSU)) ;seqlevelsStyle(readsLSU) <- "NCBI";
    
    reads <- ORFik:::convertToOneBasedRanges(readsLSU)
    reads$score <- ORFik:::readWidths(readsLSUbam)
    
    all_lengths <- sort(unique(reads$score))
    dt <- data.table()
    for(l in all_lengths){
      dt <- data.table::rbindlist(list(dt, ORFik:::metaWindow(x = reads[reads$score == l], windows = window, withFrames = T, zeroPosition = 100,
                                                              returnAs = "data.table", forceUniqueEven = T, fraction = l)))
    }
    
    
    for(s in score) {
      outdir <- p(mainFolder, "/tcp_plots/")
      sample_name <- paste0("heatmap_coverage_", stage, "_",type, "_", s) # What happens on cds ?
      outName <- paste0(outdir, sample_name, ".pdf")
      coverageHeatMap(coverage = dt, output = outName, scoring = s)
    }
  }
}


getTCPdf <- function() {
  stage <- c("64", "64", "shield", "shield")
  type <- c("MZ", "WT", "MZ", "WT")
  SSU <- c("/export/valenfs/projects/adam/TCP_seq/valen_11/64cell_mzdicer/processed_29_11_18/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
           "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/64cell_SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
           "/export/valenfs/projects/adam/TCP_seq/valen_11/Shield_mzdicer/processed_29_11_18/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
           "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/shield_only_SSU_peaks_removed_200_removed_translating_lengths_25_35.bam")
  
  LSU <- c("/export/valenfs/projects/adam/TCP_seq/valen_11/64cell_mzdicer/processed_29_11_18/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
           "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/64cell_LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
           "/export/valenfs/projects/adam/TCP_seq/valen_11/Shield_mzdicer/processed_29_11_18/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
           "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/shield_only_LSU_peaks_removed_200_selected_translating_lengths_25_35.bam")
  
  df <- data.frame(SSU, LSU, stage, type, stringsAsFactors = FALSE)
  return(df)
}

#' Make 100 bases size meta window
transcriptWindow <- function(leaders, cds, trailers, df = getTCPdf(), 
                             outdir = p(mainFolder, "/tcp_plots/mir430/normal_"),
                             scores = c("sum", "zscore"), fractions = c("SSU", "LSU")) {
  
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    
    
    readsSSU <- GRanges(readGAlignments(SSU)) ;seqlevelsStyle(readsSSU) <- "NCBI";
    readsLSU <- GRanges(readGAlignments(LSU)) ;seqlevelsStyle(readsLSU) <- "NCBI";
    
    transcriptWindowPer(leaders, cds, trailers, df[1,], outdir, scores, fractions,
                        readsSSU, readsLSU)
  }
}

transcriptWindowPer <- function(leaders, cds, trailers, df = getTCPdf()[1,], 
                                outdir = p(mainFolder, "/tcp_plots/mir430/normal_"),
                                scores = c("sum", "zscore"), fractions = c("SSU", "LSU"), 
                                read1, read2 = NULL, plot = TRUE, returnCoverage = FALSE) {
  stage <- df$stage[1]
  type <- df$type[1]
  
  seqlevelsStyle(read1) <- seqlevelsStyle(leaders);
  
  leaderSSU <- scaledWindowPositions(leaders, read1)
  leaderSSU[, `:=` (fraction = fractions[1], feature = "leaders")]
  cdsSSU <- scaledWindowPositions(cds, read1)
  cdsSSU[, `:=` (fraction = fractions[1], feature = "cds")]
  trailerSSU <- scaledWindowPositions(trailers, read1)
  trailerSSU[, `:=` (fraction = fractions[1], feature = "trailers")]
  
  if (!is.null(read2)) {
    leaderLSU <- scaledWindowPositions(leaders, read2)
    leaderLSU[, `:=` (fraction = fractions[2], feature = "leaders")]
    cdsLSU <- scaledWindowPositions(cds, read2)
    cdsLSU[, `:=` (fraction = fractions[2], feature = "cds")]
    trailerLSU <- scaledWindowPositions(trailers, read2)
    trailerLSU[, `:=` (fraction = fractions[2], feature = "trailers")]
    coverage <- rbind(leaderSSU, cdsSSU, trailerSSU, leaderLSU, cdsLSU, trailerLSU)
  } else {
    coverage <- rbind(leaderSSU, cdsSSU, trailerSSU)
  }
  
  if (plot) {
    for(s in scores) {
      sample_name <- paste0(stage, "_", type, "_", s) # What happens on cds ?
      outName <- paste(outdir, sample_name, ".pdf", sep="")
      windowCoveragePlot(coverage, output = outName, scoring = s)
    }
  }
  if (returnCoverage) return(coverage)
  return(NULL)
}

#' TIS window 150 bases
regionWindow <- function(cds, tx, df = getTCPdf(), 
                         outdir = p(mainFolder, "/tcp_plots/TIS_region/target_"),
                         upstream = 75, downstream = 74) {
  # find most translated genes
  windowsStart <- startRegion(cds, tx, TRUE, upstream = upstream, downstream = downstream)
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    # readsSSU <- readGAlignments(SSU) ;seqlevelsStyle(readsSSU) <- "NCBI";
    readsLSU <- readGAlignments(LSU) ;seqlevelsStyle(readsLSU) <- seqlevelsStyle(windowsStart)[1];
    
    
    hitMapStart <- metaWindow(readsLSU, windowsStart, withFrames = TRUE)
    
    sample_name <- paste0(stage, "_", type)
    output <- paste0(outdir, sample_name, ".pdf")
    ORFik:::pSitePlot(hitMapStart, length = paste0("TIS region ", sample_name), output = output)
  }
  return(NULL)
}

#' TIS window 150 bases
regionWindowBoth <- function(cds, tx, df = getTCPdf(), 
                         outdir = p(mainFolder, "/tcp_plots/TIS_region/target_"),
                         upstream = 75, downstream = 74, scores = c("sum", "zscore")) {
  # find most translated genes
  windowsStart <- startRegion(cds, tx, TRUE, upstream = upstream, downstream = downstream)
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    readsSSU <- readGAlignments(SSU) ;seqlevelsStyle(readsSSU) <- seqlevelsStyle(windowsStart)[1]
    readsLSU <- readGAlignments(LSU) ;seqlevelsStyle(readsLSU) <- seqlevelsStyle(windowsStart)[1]
    
    hitMapStartSSU <- metaWindow(readsSSU, windowsStart, withFrames = FALSE, returnAs = "data.table")
    hitMapStartSSU[, `:=` (fraction = "SSU", feature = "max structure")]
    hitMapStartLSU <- metaWindow(readsLSU, windowsStart, withFrames = FALSE, returnAs = "data.table")
    hitMapStartLSU[, `:=` (fraction = "LSU", feature = "max structure")]
    
    coverage <- rbind(hitMapStartSSU, hitMapStartLSU)
    
    
    for(s in scores) {
      sample_name <- paste0(stage, "_", type, "_", s) # What happens on cds ?
      outName <- paste(outdir, sample_name, ".pdf", sep="")
      windowCoveragePlot(coverage, output = outName, scoring = s)
    }
  }
  return(NULL)
}


countsPerLibraryOverTranscript <- function(tx, df = getTCPdf(), 
                                           output = p(mainFolder, "/tcp_plots/countsPerTranscript_Targets.csv")) {
  # find most translated genes
  overlaps <- data.table(txNames = names(tx))
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    readsSSU <- readGAlignments(SSU) ;seqlevelsStyle(readsSSU) <- seqlevelsStyle(tx)
    readsLSU <- readGAlignments(LSU) ;seqlevelsStyle(readsLSU) <- seqlevelsStyle(tx)
    
    overlaps <- cbind(overlaps, data.table(counts = countOverlaps(tx, readsSSU)))
    colnames(overlaps)[i + 1] <- paste0(type,"_",stage)
  }
  
  write.csv(overlaps, file = output)
  return(NULL)
}

# Input of this are reads (bam files) after alignin process

# STEPS:
#1: Find most highly expressed isoform of each gene (based on RNA-seq coverage)  
#2: Use cage to update leaders
#3: Extend improperly trimmed regions (while trimmed reads exactly match transcripts)
#4: Remove high coverage peaks in transcripts. Reads with same start + stop coord >= 200X mean transcript coverage.
#5: Remove ambiguous reads from TCP-seq 40S complexes. Ambiguous reads are those with the same length as translating 80S complexes (25-35nt in length). 
#6: 

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
    reads <- readGAlignments(bam);seqlevelsStyle(reads) <- "NCBI";
    
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
  
  # 5: Split in SSU/LSU:
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



# Standard
# SSU = "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/combined_time_points/all_three_SSU_peaks_removed_200_removed_translating_lengths_25_35.bam"
# LSU = "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/combined_time_points/all_three_LSU_peaks_removed_200_selected_translating_lengths_25_35.bam"
# stage <- "64"
# type <- "WT"


# score <- c("sum", "zscore")
# 
# gtfPath <- "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81.gtf"
# txdb <- makeTxDbFromGFF(file = gtfPath)
# seqlevelsStyle(txdb) <- "NCBI"
# 
# leaders = fiveUTRsByTranscript(txdb,use.names = T)
# cds <- cdsBy(txdb,"tx", use.names = TRUE)
# trailers = threeUTRsByTranscript(txdb, use.names = TRUE)
# tx <- exonsBy(txdb, by = "tx", use.names = TRUE)

# names <- names(leaders)[(names(leaders) %in% names(cds)) & (names(leaders) %in% names(trailers))]
# leaders <- leaders[names]
# cds <- cds[names]
# trailers <- trailers[names]




#' heatmap
tcpHeatMap <- function(txdb, df = getTCPdf()) {
  names2 <- filterTranscripts(txdb, 100, 50, 0)
  cds <- cds[names2]
  window <- ORFik:::startRegion(cds,tx, TRUE, 100, 49)
  
  for(i in 1:nrow(df)) {
    print(i)
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    readsLSUbam <- readGAlignments(LSU)
    readsLSU <- GRanges(readGAlignments(LSU)) ;seqlevelsStyle(readsLSU) <- seqlevelsStyle(window)
    
    reads <- ORFik:::convertToOneBasedRanges(readsLSU)
    reads$score <- ORFik:::readWidths(readsLSUbam)
    
    all_lengths <- sort(unique(reads$score))
    dt <- data.table()
    for(l in all_lengths){
      dt <- data.table::rbindlist(list(dt, ORFik:::metaWindow(x = reads[reads$score == l], windows = window, withFrames = T, zeroPosition = 100,
                                                              returnAs = "data.table", forceUniqueEven = T, fraction = l)))
    }
    
    
    for(s in score) {
      outdir <- p(mainFolder, "/tcp_plots/")
      sample_name <- paste0("heatmap_coverage_", stage, "_",type, "_", s) # What happens on cds ?
      outName <- paste0(outdir, sample_name, ".pdf")
      coverageHeatMap(coverage = dt, output = outName, scoring = s)
    }
  }
}


getTCPdf <- function(stage = c("64", "64", "shield", "shield"), 
                     type = c("MZ", "WT", "MZ", "WT"),
                     SSU = c("/export/valenfs/projects/adam/TCP_seq/valen_11/64cell_mzdicer/processed_29_11_18/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/64cell_SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/valen_11/Shield_mzdicer/processed_29_11_18/tidy_bams/SSU_peaks_removed_200_removed_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/shield_only_SSU_peaks_removed_200_removed_translating_lengths_25_35.bam"),
                     LSU = c("/export/valenfs/projects/adam/TCP_seq/valen_11/64cell_mzdicer/processed_29_11_18/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/64cell_LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/valen_11/Shield_mzdicer/processed_29_11_18/tidy_bams/LSU_peaks_removed_200_selected_translating_lengths_25_35.bam",
                              "/export/valenfs/projects/adam/TCP_seq/combined_bams_most_highly_expressed/shield_only_LSU_peaks_removed_200_selected_translating_lengths_25_35.bam")) {
  
  return(data.frame(SSU, LSU, stage, type, stringsAsFactors = FALSE))
}

#' Make 100 bases size meta window
transcriptWindow <- function(leaders, cds, trailers, df = getTCPdf(), 
                             outdir = p(mainFolder, "/tcp_plots/mir430/normal_"),
                             scores = c("sum", "zscore")) {
  
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    readsSSU <- GRanges(readGAlignments(SSU)) ;seqlevelsStyle(readsSSU) <- "NCBI";
    readsLSU <- GRanges(readGAlignments(LSU)) ;seqlevelsStyle(readsLSU) <- "NCBI";
    
    leaderSSU <- scaledWindowPositions(leaders, readsSSU)
    leaderSSU[, `:=` (fraction = "SSU", feature = "leaders")]
    cdsSSU <- scaledWindowPositions(cds, readsSSU)
    cdsSSU[, `:=` (fraction = "SSU", feature = "cds")]
    trailerSSU <- scaledWindowPositions(trailers, readsSSU)
    trailerSSU[, `:=` (fraction = "SSU", feature = "trailers")]
    
    leaderLSU <- scaledWindowPositions(leaders, readsLSU)
    leaderLSU[, `:=` (fraction = "LSU", feature = "leaders")]
    cdsLSU <- scaledWindowPositions(cds, readsLSU)
    cdsLSU[, `:=` (fraction = "LSU", feature = "cds")]
    trailerLSU <- scaledWindowPositions(trailers, readsLSU)
    trailerLSU[, `:=` (fraction = "LSU", feature = "trailers")]
    
    
    coverage <- rbind(leaderSSU, cdsSSU, trailerSSU, leaderLSU, cdsLSU, trailerLSU)
    
    for(s in scores) {
      sample_name <- paste0(stage, "_", type, "_", s) # What happens on cds ?
      outName <- paste(outdir, sample_name, ".pdf", sep="")
      windowCoveragePlot(coverage, output = outName, scoring = s)
    }
  }
}

#' TIS window 150 bases
regionWindow <- function(cds, tx, df = getTCPdf(), 
                         outdir = p(mainFolder, "/tcp_plots/TIS_region/target_")) {
  # find most translated genes
  windowsStart <- startRegion(cds, tx, TRUE, upstream = 75, 74)
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    # readsSSU <- readGAlignments(SSU) ;seqlevelsStyle(readsSSU) <- "NCBI";
    readsLSU <- readGAlignments(LSU) ;seqlevelsStyle(readsLSU) <- "NCBI";
    
    
    hitMapStart <- metaWindow(readsLSU, windowsStart, withFrames = TRUE)
    
    sample_name <- paste0(stage, "_", type)
    output <- paste0(outdir, sample_name, ".pdf")
    ORFik:::pSitePlot(hitMapStart, length = paste0("TIS region ", sample_name), output = output)
  }
  return(NULL)
}

countsPerLibraryOverTranscript <- function(tx, df = getTCPdf(), 
                                           output = p(mainFolder, "/tcp_plots/countsPerTranscript_Targets.csv")) {
  # find most translated genes
  overlaps <- data.table(txNames = names(tx))
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    readsSSU <- readGAlignments(SSU) ;seqlevelsStyle(readsSSU) <- seqlevelsStyle(tx);
    # readsLSU <- readGAlignments(LSU) ;seqlevelsStyle(readsLSU) <- "NCBI";
    
    overlaps <- cbind(overlaps, data.table(counts = countOverlaps(tx, readsSSU)))
    colnames(overlaps)[i + 1] <- paste0(type,"_",stage)
  }
  
  write.csv(overlaps, file = output)
  return(NULL)
}

countsPerLibraryOverTranscriptPerSubunit <- function(tx, df = getTCPdf(), 
                                           output = p(mainFolder, "/Repeats/countsPerRepeat.csv")) {
  # find most translated genes
  overlaps <- data.table(txNames = names(tx))
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    readsSSU <- readGAlignments(SSU) ;seqlevelsStyle(readsSSU) <- seqlevelsStyle(tx);
    readsLSU <- readGAlignments(LSU) ;seqlevelsStyle(readsLSU) <- seqlevelsStyle(tx);
    
    current <- data.table(counts1 = countOverlaps(tx, readsSSU), counts2 = countOverlaps(tx, readsLSU))
    colnames(current)[1] <- paste0(type,"_",stage, "_", "SSU")
    colnames(current)[2] <- paste0(type,"_",stage, "_", "LSU")
    
    overlaps <- cbind(overlaps, current)
  }
  
  write.csv(overlaps, file = output)
  return(NULL)
}

readLengthsPerLibrary <- function(df = getTCPdf()) {
  # find most translated genes
  overlaps <- data.table(id = 1)
  for (i in 1:nrow(df)) {
    print(i)
    SSU <- df$SSU[i]
    LSU <- df$LSU[i]
    stage <- df$stage[i]
    type <- df$type[i]
    
    readsSSU <- readGAlignments(SSU)
    readsLSU <- readGAlignments(LSU)
    
    current <- data.table(counts1 = length(readsSSU), counts2 = length(readsLSU))
    colnames(current)[1] <- paste0(type,"_",stage, "_", "SSU")
    colnames(current)[2] <- paste0(type,"_",stage, "_", "LSU")
    
    overlaps <- cbind(overlaps, current)
  }
  overlaps[, id := NULL]
  return(overlaps)
}
