#' Create heatmap of specified region
#'
#' Simplified input space for easier abstraction
#' @param region a region as character (default: TIS), or TSS
heatMapRegion <- function(df, region = "TIS",
                          outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/"),
                          scores = c("transcriptNormalized", "sum"), cage = NULL, format = ".png") {
  if (!(any(region %in% c("TIS", "TSS", "TTS")))) stop("region must be either TSS, TIS or TTS")
  dir.create(outdir, showWarnings = FALSE,
             recursive = TRUE)
  txdb <- loadTxdb(df)
  acLen <- 21:75
  if ("TIS" %in% region) {
    message("TIS")
    txNames <- filterTranscripts(txdb, 51, 70, 0)
    center <- loadRegion(txdb, "cds")[txNames]
    mrna <- loadRegion(txdb, "mrna")[txNames]
    upstream <- c(50, 30)
    downstream <- c(29, 69)
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TIS", shifting = c("5prime", "3prime"),
             skip.last = FALSE, acLen = acLen)
  }
  if ("TSS" %in% region) {
    message("TSS")
    txNames <- filterTranscripts(txdb, 70, 0, 0)
    center <- loadRegion(txdb, "leaders")[txNames]
    mrna <- loadRegion(txdb, "mrna")
    len <- startSites(center, keep.names = TRUE, is.sorted = TRUE) > 51
    center <- center[len]
    mrna <- mrna[names(center)]
    if (!is.null(cage)) {
      message("Using cage file to update TSS")
      center <- reassignTSSbyCage(center, cage)
      mrna <- ORFik:::downstreamFromPerGroup(mrna, startSites(center, is.sorted = TRUE))
    }
    mrna <- extendLeaders(mrna, 51)
    upstream <- c(50, 30)
    downstream <- c(29, 69)
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TSS", shifting = c("5prime", "3prime"),
             skip.last = FALSE, acLen = acLen)
  }
  if ("TTS" %in% region) {
    message("TTS")
    txNames <- filterTranscripts(txdb, 0, 51, 70)
    center <- loadRegion(txdb, "trailers")[txNames]
    mrna <- loadRegion(txdb, "mrna")[txNames]
    upstream <- c(50, 30)
    downstream <- c(29, 69)
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TTS", shifting = c("5prime", "3prime"),
             skip.last = FALSE, acLen = acLen)
  }
  return(NULL)
}

#' heatmap
#' Old default path:  = p(mainFolder, "/tcp_plots/mir430/hm_")
heatMapL <- function(region, tx, df, outdir, scores = "sum", upstream, downstream,
                     zeroPosition = upstream, acLen = NULL,
                     legendPos = "right", colors = NULL, addFracPlot = TRUE,
                     location = "TIS", shifting = NULL, skip.last = FALSE, format = ".png") {
  up <- upstream; down <- downstream
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, region, type = "bedo")

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
          if (length(zeroPosition) > 1) {
            zero <- zeroPosition[s]
          }
          print(paste(i, shift, score))
          out <- paste0(outdir, df@experiment,"_hm_", location, "_",i , "_")
          out <- ifelse(!is.null(shifting),
                        paste0(out, shift, "_", score, format),
                        paste0(out, score, format))

          tcpHeatMap_single(region, tx, reads = get(i), outdir = out,
                            shifting = shift, scores = score, upstream = up, downstream = down,
                            zeroPosition = zero, acLen = acLen,
                            legendPos = legendPos, colors = colors, addFracPlot = addFracPlot,
                            location = location, skip.last = skip.last)
        }
      }
    }
  }
}

tcpHeatMap_single <- function(region, tx, reads, outdir,
                              scores = "sum", upstream, downstream,  zeroPosition = upstream,
                              returnCoverage = FALSE, acLen = NULL, legendPos = "right",
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
                            zeroPosition = zeroPosition, scoring = scores, acceptedLengths = acLen)

  if (scores == "log2sum") scores <- "log2(sum)"
  dt$fraction <- factor(dt$fraction, levels = unique(dt$fraction),
                              labels = unique(dt$fraction))
  if (is.null(colors))
    colors <- c("white", "yellow2", "yellow3", "lightblue", "blue", "navy")

  if (all(colors == "high"))
    colors <- c("white", "yellow2", "yellow3", "lightblue", "blue", "blue", "blue", "navy", "black")

  internal_xscaler <- function(upstream, downstream, by = 10) {
    if ((upstream %% by) < by / 2) {
      return(c(seq.int(-upstream, downstream, by = by)))
    }
    return(c(seq.int(-upstream, downstream, by = by), downstream))
  }
  plot <- ggplot(dt, aes(x = position, y = fraction,
                               fill = score)) + geom_tile() +
    scale_fill_gradientn(colours = colors, name = ORFik:::prettyScoring(scores), na.value = "white") +
    xlab(paste("Position relative to", location)) +
    ylab("Protected fragment length") +
    scale_x_continuous(breaks = internal_xscaler(upstream, downstream)) +
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

#' Meta coverage of region
#' @param dfr RNA seq ORFik experiment to normalize against
regionWindow <- function(region, tx, df = getTCPdf(),
                         outdir = p(mainFolder, "/tcp_plots/TIS_region/target_"),
                         upstream = 75, downstream = 74, scores = c("sum", "zscore"),
                         colors = c("skyblue4", "orange"), allTogether = TRUE,
                         title = "cds metacoverage", dfr = NULL) {
  # find most translated genes
  windowsOne <- startRegion(region, tx, TRUE, upstream = upstream, downstream = downstream)
  size <- min(widthPerGroup(windowsOne))
  message(paste("min window size", size))
  message(paste("all window size", unique(sum(width(windowsOne)))))
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, windowsOne)
    coverage <- data.table()
    if (!allTogether) {
      stop("fix!")
    } else { # all combined
      for (i in varNames) { # For each stage
        print(i)
        hitMapOne <- metaWindow(get(i), windowsOne, scoring = "meanPos", withFrames = FALSE,
                                fraction = i, feature = "")
        coverage <- rbindlist(list(coverage, hitMapOne))
      }

      if (!is.null(dfr)) {
        coverage <- rnaNormalize(coverage, df, dfr, tx[names(region)])
        title <- paste0(title, " RNA-normalized")
      }
      for(s in scores) {
        a <- windowCoveragePlot(coverage, scoring = s, colors = colors, title = title)
        ggsave(paste0(outdir, paste0(df@experiment,"_cr_all_", s, "_", title,".png"))
               , a, height = 10)
      }
    }
  }
  if (!exists("a")) a <- NULL
  return(a)
}

#' Meta coverage of region using two groups
#' @param dfr RNA seq ORFik experiment to normalize positions against
regionWindowAll <- function(one, two, tx, df = getTCPdf(), outdir,
                            upstream = 75, downstream = 74, scores = c("sum", "zscore"),
                            seperations = c("lowly", "highly"), title = "cds metacoverage",
                            colors = c("skyblue4", "orange"), allTogether = TRUE, dfr = NULL) {
  # find most translated genes
  windowsOne <- startRegion(one, tx, TRUE, upstream = upstream, downstream = downstream)
  windowTwo <- startRegion(two, tx, TRUE, upstream = upstream, downstream = downstream)

  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, windowsOne)
    coverage <- data.table()
    if (!allTogether) {
      stop("fix!")
    } else { # all combined
      for (i in varNames) { # For each stage
        print(i)
        hitMapOne <- metaWindow(get(i), windowsOne, scoring = "meanPos", withFrames = FALSE,
                                fraction = paste0(seperations[1], "_", i), feature = "")
        hitMapTwo <- metaWindow(get(i), windowTwo, scoring = "meanPos", withFrames = FALSE,
                                fraction = paste0(seperations[2], "_", i), feature = "")

        coverage <- rbindlist(list(coverage, hitMapOne, hitMapTwo))
      }


      for(s in scores) {
        a <- windowCoveragePlot(coverage, scoring = s, colors = colors, title = title)
        ggsave(ORFik:::pasteDir(outdir, paste0(df@experiment,"_cr_all_", s, "_", title,".png"))
               , a, height = 10)
      }
    }
  }
  if (!exists("a")) a <- NULL
  return(a)
}

regionBarPlot <- function(region, tx, reads,
                          outdir, upstream = 75, downstream = 74,
                          scores = "sum", proportion = TRUE) {
  windows <- startRegion(region, tx, TRUE,
                         upstream = upstream, downstream = downstream)
  hitMap <- metaWindow(reads, windows, zeroPosition = upstream,
                       feature = "1", fraction = "123",
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

regionBarPlotAll <- function(region, tx, df, outdir, format = ".png",
                          upstream = 75, downstream = 74, scores = "sum",
                          shifting = c("5prime", "3prime"),
                          proportion = TRUE, location = "TIS") {
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
          out <- ifelse(!is.null(shifting), paste0(outdir, df@experiment,"_bp_", location, "_",i ,
                                                   "_", shift, "_", score, format),
                        paste0(outdir, df@experiment,"_bp_", location, "_", i,"_", score, format))

          reads <- convertToOneBasedRanges(get(i), method = shift, addSizeColumn = TRUE)
          regionBarPlot(region, tx, reads,
                        outdir = out,
                        upstream = up, downstream = down, scores = score, proportion = TRUE)
        }
      }
    }
  }
}

#' Find ratio of iMet to met amino acids
#' Also get basic tRNA statistics
#' @param fa location of fasta genome file which tRNA aligned to, or a biostring object of preloaded
#' @param bamFiles paths to bamFiles matched to fa. A character vector of path to bam files aligned to fa.
#' @param trna
tRNA_stats <- function(fa, bamFiles,
                       trna = "/export/valenfs/data/references/Zv10_zebrafish/tRNAscan-SE-2.0/GRcZ10_tRNA_scan_output.fa") {
  if (is(fa, "FaFile")) {
    getFasta(trna)
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

#' @param min_cutoff either a quantile if input is string[0-1], like "0.99", or numeric value if input is numeric.
removeBadTxByRegion <- function(tx, reads, upstream, downstream, median_multiplier = "0.99",
                                extension = max(0, upstream + 1), min_cutoff = "0.999") {
  region <- ORFik:::startRegion(tx, extendLeaders(tx, extension = extension), upstream = upstream, downstream = downstream)
  coverage <- coveragePerTiling(region, reads, T, as.data.table = T)
  #counts <- countOverlaps(region, reads);print(summary(counts)); print(sum(counts > value)); print(which(counts > value))
  coverage[, median_per_gene := median(count), by = genes]
  if (is(min_cutoff, "character")) {
    if (is.na(as.numeric(min_cutoff))) stop("min_cutoff must be numeric or character of numeric value")
    min_cutoff <- quantile(coverage$count, as.numeric(min_cutoff))
  }
  if (is(median_multiplier, "character")) {
    if (is.na(as.numeric(median_multiplier))) stop("median_multiplier must be numeric or character of numeric value")
    median_multiplier <- quantile(coverage$count, as.numeric(median_multiplier))
  }
  toFilter <- coverage[count > median_per_gene*median_multiplier + min_cutoff]
  print("Number of 0-score positions and ratio:")
  print(coverage[, .(zeropos_num = sum(count == 0), ratio = sum(count == 0) / .N)])
  print("Using median_multiplier and minimum cutoff of:")
  print(c(median_multiplier, min_cutoff))
  print("Read count distribution:")
  print(summary(coverage$count))
  print("Hit positions in genes that will be filtered:")
  print(toFilter)
  print("Transcripts filtered out:")
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

#' Helper function for coverage plots
#' @inheritParams transcriptWindow
#' @param returnCoverage (defualt: FALSE), return the ggplot object (TRUE)
#'  or NULL (FALSE).
#' @param title Title to give plot
#' @param plotFunction Which plot function, default: windowCoveragePlot
#' @return NULL (or ggplot object if returnCoverage is TRUE)
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
        plot <- coverageHeatMap(coverage = coverage, scoring = s)
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
#'
#'@param splitList character vector where each unique group will be an output
#'
splitMetacoverage <- function(region, splitList, df, outdir, scores = ifelse(!is.null(dfr), c("sum", "zscore"), "sum"),
                              colors = c("skyblue4", "orange"), allTogether = TRUE,
                              title = "cds metacoverage", dfr = NULL) {
  # find most translated genes
  if (!is.character(splitList) | (length(region) != length(splitList)))
    stop("splitList must be character and equal length of |region|")
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    print(paste("Experiment:", df@experiment))
    coverage <- data.table()
    varNames <- bamVarName(df)
    outputLibs(df, region)

    for (i in varNames) { # For each stage
      print(i)
      hitMapOne <- metaWindow(get(i), region, scoring = "meanPos", withFrames = FALSE,
                              fraction = i, forceUniqueEven = FALSE)
      coverage <- rbindlist(list(coverage, hitMapOne))
    }
    coverage[, feature := splitList[genes]]

    if (!is.null(dfr)) {
      coverage <- rnaNormalize(coverage, df, dfr, region,
                               normalizeMode = "feature")
      title <- paste0(title, " RNA-normalized")
    }


    for(s in scores) {
      a <- windowCoveragePlot(coverage, scoring = s, colors = colors, title = title)
      width <- 7 + length(unique(splitList))*0.7
      ggsave(paste0(outdir, paste0(df@experiment,"_smc_all_", s, "_", title,".png"))
             , a, height = 10, width = width)
    }
  }
  if (!exists("a")) a <- NULL
  return(a)
}

#' Differential group expression
#' @param cov a data.table with score and feature column
#' @param type what is the type of score ? RNA-seq, Translational efficiency etc
#' @param grouping what are the category of groups ? (stop codons, start codon, motif groups etc)
#' @param test (defualt: nonparametric)
difGroupAnalysis <- function(cov, dir, type, grouping, test = "nonparametric") {
  library(rstatix)
  ecdfplotName <- paste0(dir, type, "_", grouping, "_ecdf.png")
  boxplotName <- paste0(dir, type, "_", grouping, "_boxplot.png")
  # ecdf
  dot <- ggplot(data = cov, aes(x = score, color = feature)) +
    stat_ecdf() +
    scale_x_log10() +
    xlab(type) +
    ylab("Cumulative frequency")
  plot(dot)
  ggsave(ecdfplotName, dot, height=100, width=250, units = 'mm', limitsize = F, dpi = 300)

  dot <- ggplot(data = cov, aes(x = feature, y = log2(score))) +
    geom_violin() +
    geom_boxplot(alpha = 0.8) +
    xlab(grouping) +
    ylab(paste("Log2", type))
  plot(dot)
  ggsave(boxplotName, dot, height=100, width=250, units = 'mm', limitsize = F, dpi = 300)


  if (test == "parametric") {  # Is it normally distributed
    message("Parametric data")
    res.aov <- aov(score ~ feature, data = cov)
    print(summary(res.aov))
    print(TukeyHSD(res.aov))
    print(plot(res.aov, 2)) # NO

  }else if (test == "nonparametric") {# Non-parametric
    message("Non-parametric data")
    print(kruskal.test(score ~ feature, data = cov))
    message("fdr p-adjust")
    wc <- wilcox_test(score ~ feature, data = cov, p.adjust.method = "fdr")
    print(data.table(group1 = wc$group1, group2 = wc$group2, p.adj = wc$p.adj,
                     is.significant = wc$p.adj.signif))
    message("holm p-adjust")
    wc <- wilcox_test(score ~ feature, data = cov, p.adjust.method = "holm")
    print(data.table(group1 = wc$group1, group2 = wc$group2, p.adj = wc$p.adj,
                     is.significant = wc$p.adj.signif))
    message("bonferroni p-adjust")
    wc <- wilcox_test(score ~ feature, data = cov, p.adjust.method = "bonferroni")
    print(data.table(group1 = wc$group1, group2 = wc$group2, p.adj = wc$p.adj,
                     is.significant = wc$p.adj.signif))
    print(cov[, .(mean = mean(score), median = median(score), max = max(score),
                  ratio = .N/sum(nrow(scanning))), by = feature])
  }

  return(NULL)
}

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
