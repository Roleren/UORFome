#clusterHelper(table = cageTissuesPrediction, saveLocation = "/clustering/cageByTissuesWithPrediction.pdf", method = "binary", mainText = "Dendrogram of cage tissues clustering with uORF predictions")
# database plots

# plot info :
# 1. general stat, median, mean, count, max, mean width, max width, min width, etc.
# 2. Cage per tissue, difference (counts per tissue), unique per tissue, variance per uorf on cage
# 3. ribo seq, how does the cage change with fpkm support ? 
# 4. Expl.
plotGeneralStatistics <- function(){
  rm(list=ls())
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./DataBaseSetup.R")
  
  
  #1. general statistics:
  
  #check max/min set of uorfs, check that everything looks good
  
  
  #find how few orfs are minimum
  rm(fiveUTRs)
  getLeaders()
  numberOfTx <- length(fiveUTRs) 
  
  widthsMeanOriginal <- mean(ORFik:::widthPerGroup(fiveUTRs,F))
  originalUorfsByTx <- getUnfilteredUORFs(fiveUTRs, assignRanges = F)
  
  gr <- unlist(originalUorfsByTx, use.names = F)
  originalUorfs <- groupGRangesBy(gr, gr$names)
  originalUorfs <- removeORFsWithinCDS(originalUorfs)
  gr <- unlist(originalUorfs, use.names = F)
  originalUorfsByTx <- groupGRangesBy(gr, names(gr))
  nTxWithUorfs <- length(originalUorfsByTx)
  uorfTxRatio <- nTxWithUorfs/numberOfTx
  nOrfsOriginal <- length(originalUorfs)
  uorfsPerTx <- nOrfsOriginal/numberOfTx
  uorfsPerTxNormalized <- uorfsPerTx/uorfTxRatio
  
  
  #find how many orfs are maximum
  rm(fiveUTRs)
  fiveUTRsMax <- leaderAllSpanning()
  
  widthsMeanMax <- mean(ORFik:::widthPerGroup(fiveUTRsMax,F))
  maxUorfsByTx  <- getUnfilteredUORFs(fiveUTRsMax, assignRanges = F)
  nTxWithUorfsMax <- length(maxUorfsByTx)
  uorfTxRatioMax <- nTxWithUorfsMax/numberOfTx
  gr <- unlist(maxUorfsByTx, use.names = F)
  maxUorf <- groupGRangesBy(gr, gr$names)
  nOrfsMax <- length(maxUorf)
  uorfsPerTxMax <- nOrfsMax/numberOfTx
  uorfsMaxPerTxNormalized <- uorfsPerTxMax/uorfTxRatioMax
  
  
    
  #now data-base results:
  grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
  
  #first some sanity checks
  # ORfs found can not be more than maxOrfs
  if (length(grl) > nOrfsMax) stop("number of orfs in atlas > than maxOrfs possible!")
  # ORfs found must all be in uidMax
  uidGRL <- toUniqueIDFromGR(grl)
  if (length(uidGRL) != length(unique(uidGRL))) stop("uids as gr was not unique, something is wrong!")
  uidMax <- toUniqueIDFromGR(maxUorf)
  uidMax <- unique(uidMax)
  overlap <- uidGRL %in% uidMax
  uidOriginal <- toUniqueIDFromGR(originalUorfs)
  nuorfsOriginalUnique <- length(unique(uidOriginal)) 
  if (sum(overlap) +1  < min(length(uidGRL),length(uidMax))) stop("Not all orfs in atlas are in maxOrfs")
  #whats up with that one that is not there ????
  
      
  #width per uorf
  widths <- ORFik:::widthPerGroup(grl,F)
  medianWidths <- median(widths)
  meanWidths <- mean(widths)
  maxWidth <- max(widths)
  minWidth <- min(widths)
  
  #width per exon
  widths <- width(unlist(grl, use.names = F))
  
  exonmedianWidths <- median(widths)
  exonmeanWidths <- mean(widths)
  maxWidth <- max(widths)
  minWidth <- min(widths)
  
  #exon stats
  numberOfexons <- length(widths)
  exonsPerGroup <- GroupGRangesExonsPerGroup(grl)
  medianExons <- median(exonsPerGroup)
  meanExons <- mean(exonsPerGroup)
  
  #some max checks
  a <- findOverlaps(grl,fiveUTRsMax, type = "within")
  index <- which.max(runLength(Rle(to(a))))
  leaderWithMostUORFs <- fiveUTRsMax[index]
  
  #check width changes from cage, the effect
  
  
  df <- data.frame()
  
  library(ggplot2)
  
  #make table of stats
  
  
  
  #plots: 
  
  #1. width histogram
  
  widthsGRL <- ORFik:::widthPerGroup(grl,F)
  df <- data.frame(widths = widthsGRL)
  lengthPlot <- ggplot() +
    geom_histogram(data = df, aes(x = log10(widths)),color="orange") +
    xlab("length of uorfs (log)") + ylab("number of counts")  +
    labs(title="Length of UORFs")
  lengthPlot
}

#2. cage average hits per tissue
plotCageTissues <- function(){
  tissueAtlas <- readTable("tissueAtlasByCage", with.IDs = F)
  
  colsums <- colSums(tissueAtlas)
  badIds <- which(colsums == 0) # <-no matches found for tissue in cage
  names(badIds) <- NULL
  tissueAtlas[,badIds] <- NULL
  colsums <- colsums[colsums > 0]
  
  df <- data.frame(colsums = colsums)
  line <- data.frame(y = c(nrow(tissueAtlas), nrow(tissueAtlas)), x = c(0, nrow(df)) )
  distributionPlot <- ggplot()  +
    geom_histogram(data = df, aes(x = 1:length(colsums), y = colsums),color="orange",stat="identity") +
    xlab("tissue") + ylab("uorfs in tissue")  +
    labs(title="cage hits of UORFs")  +  geom_line(aes(x = line$x, y = line$y), color = "blue")
  distributionPlot
  
  
  
  
  uorfsInAllTissues <- sum(rowSums(tissueAtlas) == (length(tissueAtlas)))
  uorfsRatioTissue <- uorfsInAllTissues/nrow(tissueAtlas)
  
  uorfsToUniqueTissue <- rowSums(tissueAtlas) == 1
  
  uorfsToUniqueTissueSum <- sum(uorfsToUniqueTissue)
  uorfsAtleastTwoTissues <- nrow(tissueAtlas) - uorfsInAllTissues - uorfsToUniqueTissueSum
  # pie chart of types
  df <- data.frame(value = c(uorfsInAllTissues, uorfsToUniqueTissueSum
                              ,uorfsAtleastTwoTissues),
                   group = c("in all tissues", "in 1 tissue",
                             "at least 2 tissues, not all"))
  piePlot <-  ggplot(df, aes(x="", y=value, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + 
    geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                  label = scales::percent(value/sum(value))), size=4)
    
  piePlot
  
  # variation <- colsums - mean(colsums)
  # goodIDs <- which(variation > quantile(variation, 0.90))
  # colsums <- colsums[goodIDs]
  
  # 5 and 95% quantile plots
  q95 <- which(colsums > quantile(colsums, 0.95))
  q05 <- which(colsums < quantile(colsums, 0.05))
  c95 <- colsums[q95]
  c05 <- colsums[q05]
  cs <- c(c95, c05)
  
  df <- data.frame(cs = cs, names = names(cs))
  line <- data.frame(y = c(mean(cs), mean(cs)), x = c(0, nrow(df)) )
  line2 <- data.frame(y = c(uorfsInAllTissues, uorfsInAllTissues),
                      x = c(0, nrow(df)) )
  distributionPlot <- ggplot()  +
    geom_histogram(data = df, aes(x = factor(names), y = cs),color="orange",stat="identity") +
    xlab("tissue") + ylab("uorfs in tissue")  +
    labs(title="cage hits of UORFs")  +  geom_line(aes(x = line$x, y = line$y), color = "blue") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      geom_line(aes(x = line2$x, y = line2$y), color = "red")
  distributionPlot
  
  
  # Which tissues does have the most unique uORFs
  uniqueTissueUorfs <- tissueAtlas[uorfsToUniqueTissue,]
  uniqueTissueUorfsColSums <- colSums(uniqueTissueUorfs)
}

#3. Show that fiveUTRs acctually gets reduced in many samples
# stacked bar for 3 points, unchanged, down, up
PlotLeaderWidthChanges <- function() {
  
  # load("leaderWidthChanges.rdata")
  # df <- as.data.frame(dt)
  # widthChangeCagePlotSame <- ggplot()  +
  #   geom_histogram(data = df, aes(x = 1:nrow(df), y = same),color="orange",stat="identity") +
  #   xlab("cage experiment") + ylab("width")  +
  #   labs(title="width changes per cage experiment: same")  
  # widthChangeCagePlotSame
  # 
  # widthChangeCagePlotBigger <- ggplot()  +
  #   geom_histogram(data = df, aes(x = 1:nrow(df), y = bigger),color="orange",stat="identity") +
  #   xlab("cage experiment") + ylab("width")  +
  #   labs(title="width changes per cage experiment: bigger") 
  # widthChangeCagePlotBigger
  # 
  # widthChangeCagePlotSmaller <- ggplot()  +
  #   geom_histogram(data = df, aes(x = 1:nrow(df), y = smaller),color="orange",stat="identity") +
  #   xlab("cage experiment") + ylab("width")  +
  #   labs(title="width changes per cage experiment: smaller") 
  # widthChangeCagePlotSmaller
  # 
  # widthChangeCagePlotMeanBigger <- ggplot()  +
  #   geom_histogram(data = df, aes(x = 1:nrow(df), y = meanBigger),color="orange",stat="identity") +
  #   xlab("cage experiment") + ylab("width")  +
  #   labs(title="width changes per cage experiment: mean of bigger") 
  # widthChangeCagePlotMeanBigger
  
  load("CageFiveUTRs.rdata")
  getLeaders()
  
  widthCage <- widthPerGroup(CageFiveUTRs)
  widthNormal <- widthPerGroup(fiveUTRs)
  bigger <- sum(widthCage > widthNormal)
  equal <- sum(widthCage == widthNormal)
  smaller <- sum(widthCage < widthNormal)
  
  df <- data.frame(value = c(bigger, equal
                             ,smaller),
                   group = c("bigger", "equal",
                             "smaller"))
  piePlot <-  ggplot(df, aes(x="", y=value, fill=group)) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + 
    geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                  label = scales::percent(value/sum(value))), size=4)
  
  piePlot
  
  # total changes
  
  load("leaderWidthChanges.rdata")
  
  df <- data.frame(values = c(dt$bigger, dt$smaller, dt$same),
                   names = c(rep("bigger", nrow(dt)),rep("smaller", nrow(dt)),rep("same", nrow(dt)) ))
  colnames(df) <- c("names", "values")
  g <- ggplot(df, aes(x = names, y = log10(values),  fill = names)) + 
    geom_violin() + 
    theme(axis.text.x = element_text(angle=0, size = 14), axis.text.y = element_text(size = 12), 
          plot.title = element_text(size = 15)) + 
    labs(title="Leader size changes by CAGE", 
         subtitle="For all 1829 CAGE-seq experiments", x = "", y = "Amount")
  plot(g)
}

#' plot ribo support differences
plotRiboDifferences <- function(){
  ribo <- readTable("RiboByTissueTF", with.IDs = F)
  riboColSums <- colSums(ribo)
  df.ribo <- as.data.frame(riboColSums)
  
  riboSupportTissue <- ggplot()  +
    geom_histogram(data = df.ribo, aes(x = rownames(df.ribo), y = riboColSums, fill=rownames(df.ribo)), width = 0.5, stat="identity") +
    xlab("ribo tissue") + ylab("uorfs supported")  +
    labs(title="Ribo-seq support, number of uorfs per tissue", 
         subtitle="FPKM > 1 in 2 experiments for inclusion")  
  riboSupportTissue
  
  # case of brain, which of the ribo-seq supported reads
  # are not in standard transcriptome ?
  
  
  brainRibo <- ribo$brain
  brainCage <- tissueAtlas$brain
  supportedUorfs <- sum(brainRibo)
  indeces <- which(brainRibo == 1)
  idsSupported <- ribo$uorfID[indeces]
  uniqueSupportedForBrain <- idsSupported %in% uidOriginal
  idsResult <- sum(uniqueSupportedForBrain)
  
  
  
  # get unique brain grl
  # find genes to check for interesting features
  uidsBrain <- idsSupported[!uniqueSupportedForBrain]
  grlBrain <- toGRFromUniqueID(uidsBrain)
  overlaps <-  findOverlaps(grlBrain, fiveUTRsMax, type = "within")
  txname <- names(fiveUTRsMax[to(overlaps)])
  
  getAllTranscriptLengths()
  geneID <- allLengths$gene_id[allLengths$tx_name == txname[2]]
  
  
  # now find a brain example to use:
  # using HTR3A gene, this made no hit, so not valid
  
  allLengths[allLengths$gene_id == "ENSG00000166736",]
  findOverlaps(grl, fiveUTRsMax["ENST00000504030"], type = "within")
  grl[27154] # <- verified as uorf
  findOverlaps(grlBrain, fiveUTRsMax["ENST00000504030"], type = "within")
  
  # lets try huntington gene HTT
  # ENSG00000197386
  allLengths[allLengths$gene_id == "ENSG00000197386",]
  findOverlaps(grlBrain, fiveUTRsMax["ENST00000355072"], type = "within")
  # no hit, lets try with all allowed
  grlBrainAll <- toGRFromUniqueID(idsSupported)
  findOverlaps(grlBrainAll, fiveUTRsMax["ENST00000355072"], type = "within")
  # found a match, but is not unique to brain
  indecesKid <- which(ribo$kidney == 1)
  idsSupportedKid <- ribo$uorfID[indecesKid]
  grlKid<- toGRFromUniqueID(idsSupportedKid)
  findOverlaps(grlKid, fiveUTRsMax["ENST00000355072"], type = "within")
  
}

plotCageStatistics <- function(){
  cage <- fread.bed(standardCage)
  tab <- table(cage$score)[-1][-1]
  
  a <- lapply(cageFiles[1:50], function(x){
    return(sum(table(fread.bed(p(cageFolder, x))$score)[-1]))
  })
  b <- unlist(a)
  mean(b)
}


plotBasicFeatures <- function(){
  library(RColorBrewer)
  starts <- readTable("StartCodons")
  
  r <- starts$startCodon
  tab <- table(r)
  
  df <- data.frame(tab)
  colnames(df) <- c("names", "values")
  g <- ggplot(df, aes(x = names, y = values)) + 
    geom_bar(aes(fill=1:length(values)), width = 0.3, stat = "identity") + 
    theme(axis.text.x = element_text(angle=65, size = 14), axis.text.y = element_text(size = 12), 
          plot.title = element_text(size = 15)) + 
    labs(title="Start codon usage", 
         subtitle="All 1 base variations from ATG", x = "Start codon", y = "Amount")
  plot(g)
  
  # for riboseq
  riboTab <- readTable("Ribofpkm", with.IDs = F)
  rATG <- mean(rowMeans(riboTab[r == "ATG",]))
  rCTG <- mean(rowMeans(riboTab[r == "CTG",]))
  rTTG <- mean(rowMeans(riboTab[r == "TTG",]))
  rGTG <- mean(rowMeans(riboTab[r == "GTG",]))
  rAAG <- mean(rowMeans(riboTab[r == "AAG",]))
  rAGG <- mean(rowMeans(riboTab[r == "AGG",]))
  rACG <- mean(rowMeans(riboTab[r == "ACG",]))
  rATC <- mean(rowMeans(riboTab[r == "ATC",]))
  rATA <- mean(rowMeans(riboTab[r == "ATA",]))
  rATT <- mean(rowMeans(riboTab[r == "ATT",]))
  
  df <- data.frame(values <- as.numeric(c(rATG,rCTG,rTTG,rGTG,rAAG,rAGG,rACG,rATC,rATA,rATT)), 
                   names = c("ATG","CTG","TTG","GTG","AAG","AGG",
                             "ACG","ATC","ATA","ATT"))
  rg <-  ggplot(df, aes(x = names, y = values)) + 
    geom_bar(aes(fill=1:length(values)), width = 0.5, stat = "identity") + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
    labs(title="Riboseq FPKM of ORFs by Start codon", 
         subtitle="All 1 base variations from ATG")
  plot(rg)
  #for te
  
  te <- readTable("teFiltered", with.IDs = F)
  teATG <- mean(rowSums(te[r == "ATG",])/ncol(te))
  teCTG <-  mean(rowSums(te[r == "CTG",])/ncol(te))
  teTTG <-  mean(rowSums(te[r == "TTG",])/ncol(te))
  teGTG <-  mean(rowSums(te[r == "GTG",])/ncol(te))
  teAAG <-  mean(rowSums(te[r == "AAG",])/ncol(te))
  teAGG <-  mean(rowSums(te[r == "AGG",])/ncol(te))
  teACG <-  mean(rowSums(te[r == "ACG",])/ncol(te))
  teATC <-  mean(rowSums(te[r == "ATC",])/ncol(te))
  teATA <-  mean(rowSums(te[r == "ATA",])/ncol(te))
  teATT <-  mean(rowSums(te[r == "ATT",])/ncol(te))
  
  df <- data.frame(values <- as.numeric(c(teATG,teCTG,teTTG,teGTG,
                                          teAAG,teAGG,teACG,teATC,teATA,teATT)), 
                   names = c("ATG","CTG","TTG","GTG","AAG","AGG",
                             "ACG","ATC","ATA","ATT"))
  gg <- ggplot(df, aes(x = names, y = values)) + 
    geom_bar(aes(fill=1:length(values)), width = 0.5, stat = "identity") + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
    labs(title="Start codon TE", 
         subtitle="All 1 base variations from ATG")
  plot(gg)
  
  kozak <- readTable("kozak", with.IDs = F)
  
  kATG <- mean(kozak$kozak[r == "ATG"])
  kCTG <- mean(kozak$kozak[r == "CTG"])
  kTTG <- mean(kozak$kozak[r == "TTG"])
  kGTG <- mean(kozak$kozak[r == "GTG"])
  kAAG <- mean(kozak$kozak[r == "AAG"])
  kAGG <- mean(kozak$kozak[r == "AGG"])
  kACG <- mean(kozak$kozak[r == "ACG"])
  kATC <- mean(kozak$kozak[r == "ATC"])
  kATA <- mean(kozak$kozak[r == "ATA"])
  kATT <- mean(kozak$kozak[r == "ATT"])
  
  df <- data.frame(values <- as.numeric(c(kATG,kCTG,kTTG,kGTG,kAAG,kAGG,kACG,kATC,kATA,kATT)), 
                   names = c("ATG","CTG","TTG","GTG","AAG","AGG",
                             "ACG","ATC","ATA","ATT"))
  kg <- ggplot(df, aes(x = names, y = values)) + 
    geom_bar(aes(fill=1:length(values)), width = 0.5, stat = "identity") + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
    labs(title="Start codon mean Kozak scores", 
         subtitle="All 1 base variations from ATG")
  plot(kg)
  
  # stop codons all
  stops <- readTable("StopCodons")
  s <- stops$stopCodon
  TAG <- sum(s == "TAG")
  TGA <- sum(s == "TGA")
  TAA <- sum(s == "TAA")
  
  
  df <- data.frame(values <- as.numeric(c(TAG, TGA, TAA)), 
                   names = c("TAG", "TGA", "TAA"))
  stg <- ggplot(df, aes(x = names, y = values)) + 
    geom_bar(aes(fill=1:length(values)), width = 0.5, stat = "identity") + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
    labs(title="Stop codon usage")
  plot(stg)
  
  # stop codons predicted
  s <- s[uniqueOrder]
  s <- s[finalCagePred]
  tab <- table(s)
  
  df <- data.frame(tab)
  colnames(df) <- c("names", "values")
  
  g <- ggplot(df, aes(x = names, y = values)) + 
    geom_bar( width = 0.3, stat = "identity") + 
    theme(axis.text.x = element_text(angle=30, size = 14), axis.text.y = element_text(size = 12), 
          plot.title = element_text(size = 15)) +
    scale_color_brewer(palette="Dark2") +
    labs(title="Stop codon usage", 
         subtitle="", x = "Stop codon codon", y = "Amount")
  plot(g)
}

#1. Check cds te of predicted vs not predicted
validatePredictionPlot <- function(){
  grl <- getUorfsInDb()
  getCDS()
  cdsTE <- readTable("cdsTETissueMean", with.IDs = T)
  link <- readTable("linkORFsToTx")
  order <- ORFik:::uniqueOrder(grl)
  finalCagePred <- readTable("allUorfsByCageAndPred")
  f <- finalCagePred[order]
  positiveCDS <- link$txNames[f]
  negCDS <- link$txNames[!f]
  
  #positiveCDS <- positiveCDS[!(positiveCDS %in% unique(negCDS))]
  negCDS <- negCDS[!(negCDS %in% unique(positiveCDS))]
  
  posTE <- cdsTE[cdsTE$txNames %in% positiveCDS, ]
  negTE <- cdsTE[cdsTE$txNames %in% negCDS, ]
  
  t <- data.frame(counts = c(rowMeans(posTE[,2:ncol(posTE)]),
                             rowMeans(cdsTE[,2:ncol(cdsTE)]),
                             rowMeans(negTE[,2:ncol(negTE)])), 
                  variable = c(rep("Predicted translated", nrow(posTE)),
                               rep("Average", nrow(cdsTE)),
                               rep("Predicted non-translated", nrow(negTE))))
  #t <- t[t$counts > 1.1, ]
  values <- 1:length(unique(t$variable))
  hp <- ggplot(t, aes(variable, counts, fill = variable)) +    
    geom_boxplot() +
    theme(axis.text.x=element_text(angle=0)) + 
    labs(title="Te of CDS by uORF prediction", 
         subtitle="", y = "CDS TE", x = "CDS uORF grouping") + 
    ylim(0, 5) 
  
  plot(hp)
  
  p <- ggplot(t, aes(x=variable, y=log10(counts), fill = variable)) + 
    geom_violin() + 
    labs(title="Te of CDS by uORF prediction", 
         subtitle="", y = "CDS TE (log10)", x = "CDS uORF grouping") + 
    geom_boxplot(width=0.1) + 
    scale_color_brewer(type='qual', palette="Blues") + 
    geom_hline(yintercept = log10(median(t$counts[t$variable == "Predicted translated"]))) 
  p
  
  
  # start codon usage
  tab <- table(StartCodons$startCodon[uniqueOrder][finalCagePred])
  df <- data.frame(counts = tab)
  colnames(df) <- c("variable", "counts")
  values <- 1:nrow(df)
  hp <- ggplot(df) +
    geom_bar(aes(y = counts, x = variable, fill=values), width = 0.3, stat = "identity") +
    theme(axis.text.x = element_text(angle=45, vjust=0.6)) + 
    labs(title="Start codon usage of uORF prediction with Cage filter", 
         subtitle="All 1 base variations from ATG", x = "start codons", y = "counts")
  plot(hp)
  
  #kozak
  kozak <- readTable("kozak", with.IDs = F)
  kozak <- kozak$kozak[uniqueOrder]
  cdsKozak <- readTable("cdsKozak", with.IDs = F)
  counts <- c(mean(kozak[finalCagePred]), mean(cdsKozak$cdsKozak), mean(kozak), mean(kozak[!finalCagePred]))
  variable <- c("predicted true uORFs", "CDS", "all uORFs", "all predicted false uORFs")
  df <- data.frame(counts, variable)
  values <- 1:nrow(df)
  hp <- ggplot(df) +
    geom_bar(aes(y = counts, x = variable, fill=values), width = 0.2, stat = "identity") +
    theme(axis.text.x = element_text(angle=45)) + 
    labs(title="Average kozak-score for different ORF types", 
         subtitle="", x = "ORF type", y = "Kozak") + 
    guides(fill = F)
  plot(hp)
}


