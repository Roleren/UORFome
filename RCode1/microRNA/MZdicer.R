getMir430 <- function(path = "/export/valenfs/projects/adam/TCP_seq/data_files/targetscanFish_180829_tidy.csv", 
                      filter = -0.20) {
  targets <- read.csv(path)
  targets <- targets[!(is.na(targets$miR430_score) | is.nan(targets$miR430_score)),]
  targets <- targets[targets$miR430_score <= filter,]
  genes <- targets$gene_id
  return(genes)
}

getMir430Tx <- function(path = "/export/valenfs/projects/adam/TCP_seq/data_files/targetscanFish_180829_tidy.csv", 
                      filter = -0.20, txdb, tx) {
  # miR430 target genes
  genes <- getMir430(path, filter)
  
  # Get miR430 target transcripts
  len <- GenomicFeatures::transcriptLengths(txdb)
  match <- len[len$gene_id %in% genes,]
  match <- match[match$tx_name %in% names(tx),]
  txNames <- match$tx_name
  print(paste("mir targets length: ", length(txNames)))
  return(txNames)
}

getBazziniList <- function(file = "/export/valenfs/projects/Hakon/mir430/data/Bazzini_Supplementary_DataSet1.txt", 
                           type = "tx") {
  x <- fread(cmd = paste0("grep '^>' ",file), header = FALSE)
  if (type == "tx") {
    x <- x$V4 # Transcripts
    x <- sub(pattern = ")", x = x, "")
  } else { # gene
    x <- x$V1 
    x <- sub(pattern = ">", x = x, "")
  }
  return(x)
}

difPlot <- function(dd, filename = "/export/valenfs/projects/Hakon/mir430/fold_changes/2dFigure_bazzini.png",
                    type = "RFP", lim = 6, logIt = FALSE, lineScoring = "median") {
  if (logIt) lim <- 15
  d <- copy(dd)
  
  d$cols <- c("gray", "red")[as.integer(d$ismiRNA) + 1]
  d <- rbind(d[ismiRNA == F, ], d[ismiRNA == T, ])
  if (logIt) d$RNA <- convertLog2(d$RNA)
  
  if (type == "RFP") { # Ribo-seq
    if (logIt) d$RFP <- convertLog2(d$RFP)
    if (lineScoring == "median") {
      means <- d[,.(RNAs = median(RNA), RFPs = median(RFP)), by = .(ismiRNA, stage)]
    } else {
      means <- d[,.(RNAs = mean(RNA), RFPs = mean(RFP)), by = .(ismiRNA, stage)]
    }
    
    means$ismiRNA <- factor(means$ismiRNA)
    means$col <- c("black", "red")[means$ismiRNA]
    
    p <- ggplot() +
      geom_point(data = d,aes(x = RNA, y = RFP, alpha = 0.1), color = factor(d$cols)) + 
      ylab(expression("WT vs MZ dicer "*Delta*" RFP")) + 
      geom_hline(data = means, aes(yintercept =  RFPs, color = col))
    
  } else { # RCP
    if (logIt) d$LSU <- convertLog2(d$LSU)
    if (logIt) d$SSU <- convertLog2(d$SSU)
    if (lineScoring == "median") {
      means <- d[,.(RNAs = median(RNA), SSUs = median(SSU), LSUs = median(LSU)), by = .(ismiRNA, stage)]
    } else {
      means <- d[,.(RNAs = mean(RNA), SSUs = mean(SSU), LSUs = mean(LSU)), by = .(ismiRNA, stage)]
    }
    means$ismiRNA <- factor(means$ismiRNA)
    means$col <- c("black", "red")[means$ismiRNA]
    
    p <- ggplot() +
      geom_point(data = d,aes(x = RNA, y = LSU, alpha = 0.1), color = factor(d$cols)) + 
      ylab(expression("WT vs MZ dicer "*Delta*" LSU")) + 
      geom_hline(data = means, aes(yintercept =  LSUs, color = col))
    p  
    
    pp <- ggplot() +
      geom_point(data = d,aes(x = RNA, y = SSU, alpha = 0.1), color = factor(d$cols)) + 
      ylab(expression("WT vs MZ dicer "*Delta*" SSU")) +
      xlim(-lim, lim) +
      ylim(-lim, lim) + 
      xlab(expression("WT vs MZ dicer "*Delta*" mRNA")) + 
      facet_grid( ~ stage, scales = "free") + 
      geom_vline(data = means, aes(xintercept = RNAs, color = col)) + 
      geom_hline(data = means, aes(yintercept =  SSUs, color = col)) +
      scale_fill_manual(values = c("gray", "red"),
                        name = "", aesthetics = "cols") + 
      scale_fill_manual(values = c("black", "red"),
                        name = "", aesthetics = "col") + 
      theme(legend.position = "none")
      
    print(pp)
    ggsave(pp, filename = filename[2], width = 10, height = 4)
  }
  p <- p + 
    xlim(-lim, lim) +
    ylim(-lim, lim) + 
    xlab(expression("WT vs MZ dicer "*Delta*" mRNA")) + 
    facet_grid( ~ stage, scales = "free") + 
    geom_vline(data = means, aes(xintercept = RNAs, color = col)) + 
    scale_fill_manual(values = c("gray", "red"),
                      name = "", aesthetics = "cols") + 
    scale_fill_manual(values = c("black", "red"),
                      name = "", aesthetics = "col") + 
    theme(legend.position = "none")
  
  ggsave(p, filename = filename[1], width = 10, height = 4)
  return(p)  
}

convertLog2 <- function(col) {
  coll <- col
  col[(col > 0)] <- log2(col[(col > 0)])
  col[(coll < 0)] <- - log2(abs(col[(coll < 0)]))
  return(col)
}
