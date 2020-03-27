#' Comparison with uORFome of McGilivray et al.
#'
#' Uses the supplemental data from article
#'
#' Remap tool:
#' https://www.ncbi.nlm.nih.gov/genome/tools/remap
#' @importFrom openxlsx read.xlsx
verifyOthersWorks <- function(work = "/export/valenfs/projects/uORFome/Mcgilivray_Supplemental_Data_Tables.xlsx") {
  dataBaseFolder <- p(mainFolder,"/dataBase")
  setwd(dataBaseFolder)
  uorfDB <- createDataBase("uorfCatalogue.sqlite")
  # old grch37 annotation, top 10% uORFs
  d <- openxlsx::read.xlsx(work, sheet = 6, colNames = T, startRow = 3)
  
  strand <- d$strand
  startSite <- as.integer(d$start_coordinate)
  stopSite <- as.integer(d$end_coordinate)
  transcript <- sub("\\..*","",d$uORF_ID)
  chromosome <- d$chromosome
  
  sta <- startSite
  sta[strand == "-"] <- stopSite[strand == "-"]
  sto <- stopSite
  sto[strand == "-"] <- startSite[strand == "-"]
  getLeaders()
  goodHits <- which(transcript %in% names(fiveUTRs))
  bed <- GRanges(seqnames = chromosome, ranges = IRanges(sta, sto), strand = strand, name = transcript, score = seq.int(1, length(chromosome)))
  
  rtracklayer::export.bed(object = bed, con = "to_be_mapped_mcgilivray_uorfs.bed") # now do conversion
  print("@ NCBI Genome Remapping Service GrCH37 p.13 -> 38 no patch")
  print("File: to_be_mapped_mcgilivray_uorfs.bed")
  print("Go to: https://www.ncbi.nlm.nih.gov/genome/tools/remap")
  my.name <- readline(prompt="Press enter when you finished remapping and downloaded annotation data: ")
  if (my.name == "") {
    message("input empty, loading standard file!")
    g <- rtracklayer::import.bed(con = "/export/valenfs/projects/uORFome/dataBase/remapped_to_be_mapped_mcgilivray_uorfs.bed")
  } else {
    g <- rtracklayer::import.bed(con = my.name)
  }
  g <- g[order(g$score)]
  g$original_id <- g$score
  g$peptide_score <- as.numeric(d$peptide_score[g$original_id])
  g$score <- NULL
  gr <- g
  g <- groupGRangesBy(g, seq(length(g)))
  if( length(strand) != length(startSite) | length(strand) != length(transcript)){
    stop("wrong input readings for others")
  }
  # Test their start codons:
  getFasta()
  codons <- ORFik:::startCodons(g, is.sorted = T)
  hit.seqnames <- seqnamesPerGroup(codons, F) %in% seqlevels(fa)
  startcods <- ORFik:::txSeqsFromFa(codons[hit.seqnames],
                                    fa, is.sorted = TRUE, keep.names = FALSE)
  original.startcods <- d$start_codon[gr$original_id][hit.seqnames]
  length(startcods); length(original.startcods)
  table(startcods);table(original.startcods)
  bad.mappings <- which(startcods != original.startcods)
  g <- g[hit.seqnames[-bad.mappings]]; gr <- gr[hit.seqnames[-bad.mappings]]
  length(g); length(gr)
  
  grl <- getUorfsInDb()
  
  #how many did we find ?
  startsOur <- startSites(grl, is.sorted = T, keep.names = FALSE)
  stopsOur <- stopSites(grl, is.sorted = T, keep.names = FALSE)
  strandOur <- strandPerGroup(grl, keep.names = F)
  chromosomeOurs <- seqnamesPerGroup(grl, F)
  
  ours <- paste(startsOur, stopsOur, chromosomeOurs, strandOur)
  theirs <- paste(startSites(g, is.sorted = T, keep.names = FALSE),
                  stopSites(g, is.sorted = T, keep.names = FALSE),
                  seqnamesPerGroup(g, F), strandPerGroup(g, keep.names = F))
  hitsOurs <- which(ours %in% theirs)
  hitsTheirs <- which(theirs %in% ours)
  
  print(paste("Mapping step found from theirs, total uORFs of:", length(hitsOurs), "that is:", length(hitsOurs) /length(g), "of full set"))
  #how many match our prediction
  
  # Using only predicted uORFs
  all.scores <- readTable("finalCAGEuORFPrediction")$Matrix
  predicted <- which(all.scores == 1)
  ours.pred <- ours[predicted]
  
  hitsOurs.pred <- which(ours.pred %in% theirs)
  hitsTheirs.pred <- which(theirs %in% ours.pred)
  print(paste("Mapping step found from theirs, total uORFs of:", length(hitsOurs.pred), "that is:", length(hitsOurs.pred) /length(g), "of full set"))
  # What is uORFome prediction score of subset of ours in theirs vs other part
  summary(all.scores[hitsOurs]);summary(all.scores[-hitsOurs])
  plot_comparison()
}

#' Comparison with uORFome of McGilivray et al. all positive!
#'
#' Uses the supplemental data from article, sheet 5
#'
#' Remap tool:
#' https://www.ncbi.nlm.nih.gov/genome/tools/remap
#' @importFrom openxlsx read.xlsx
verifyOthersWorks <- function(work = "/export/valenfs/projects/uORFome/Mcgilivray_Supplemental_Data_Tables.xlsx") {
  dataBaseFolder <- p(mainFolder,"/dataBase")
  setwd(dataBaseFolder)
  uorfDB <- createDataBase("uorfCatalogue.sqlite")
  # old grch37 annotation, top 10% uORFs
  d <- openxlsx::read.xlsx(work, sheet = 5, colNames = T, startRow = 3)
  
  strand <- d$strand
  startSite <- as.integer(d$start_coordinate)
  stopSite <- as.integer(d$end_coordinate)
  transcript <- sub("\\..*","",d$uORF_ID)
  chromosome <- d$chromosome
  
  sta <- startSite
  sta[strand == "-"] <- stopSite[strand == "-"]
  sto <- stopSite
  sto[strand == "-"] <- startSite[strand == "-"]
  getLeaders()
  goodHits <- which(transcript %in% names(fiveUTRs))
  bed <- GRanges(seqnames = chromosome, ranges = IRanges(sta, sto), strand = strand, name = transcript, score = seq.int(1, length(chromosome)))
  
  rtracklayer::export.bed(object = bed, con = "all_to_be_mapped_mcgilivray_uorfs.bed") # now do conversion
  print("@ NCBI Genome Remapping Service GrCH37 p.13 -> 38 no patch")
  print("File: to_be_mapped_mcgilivray_uorfs.bed")
  print("Go to: https://www.ncbi.nlm.nih.gov/genome/tools/remap")
  my.name <- readline(prompt="Press enter when you finished remapping and downloaded annotation data: ")
  if (my.name == "") {
    message("input empty, loading standard file!")
    g <- rtracklayer::import.bed(con = "/export/valenfs/projects/uORFome/dataBase/remapped_all_to_be_mapped_mcgilivray_uorfs.bed")
  } else {
    g <- rtracklayer::import.bed(con = my.name)
  }
  g <- g[order(g$score)]
  g$original_id <- g$score
  g$peptide_score <- as.numeric(d$peptide_score[g$original_id])
  g$score <- NULL
  gr <- g
  g <- groupGRangesBy(g, seq(length(g)))
  if( length(strand) != length(startSite) | length(strand) != length(transcript)){
    stop("wrong input readings for others")
  }
  # Test their start codons:
  getFasta()
  codons <- ORFik:::startCodons(g, is.sorted = T)
  hit.seqnames <- seqnamesPerGroup(codons, F) %in% seqlevels(fa)
  startcods <- ORFik:::txSeqsFromFa(codons[hit.seqnames],
                                    fa, is.sorted = TRUE, keep.names = FALSE)
  original.startcods <- d$start_codon[gr$original_id][hit.seqnames]
  length(startcods); length(original.startcods)
  table(startcods);table(original.startcods)
  bad.mappings <- which(startcods != original.startcods)
  g <- g[hit.seqnames[-bad.mappings]]; gr <- gr[hit.seqnames[-bad.mappings]]
  length(g); length(gr)
  
  grl <- getUorfsInDb()
  
  #how many did we find ?
  startsOur <- startSites(grl, is.sorted = T, keep.names = FALSE)
  stopsOur <- stopSites(grl, is.sorted = T, keep.names = FALSE)
  strandOur <- strandPerGroup(grl, keep.names = F)
  chromosomeOurs <- seqnamesPerGroup(grl, F)
  
  ours <- paste(startsOur, stopsOur, chromosomeOurs, strandOur)
  theirs <- paste(startSites(g, is.sorted = T, keep.names = FALSE),
                  stopSites(g, is.sorted = T, keep.names = FALSE),
                  seqnamesPerGroup(g, F), strandPerGroup(g, keep.names = F))
  hitsOurs <- which(ours %in% theirs)
  hitsTheirs <- which(theirs %in% ours)
  
  print(paste("Mapping step found from theirs, total uORFs of:", length(hitsOurs), "that is:", length(hitsOurs) /length(g), "of full set"))
  #how many match our prediction
  
  # Using only predicted uORFs
  all.scores <- readTable("finalCAGEuORFPrediction")$Matrix
  predicted <- which(all.scores == 1)
  ours.pred <- ours[predicted]
  
  hitsOurs.pred <- which(ours.pred %in% theirs)
  hitsTheirs.pred <- which(theirs %in% ours.pred)
  print(paste("Mapping step found from theirs, total uORFs of:", length(hitsOurs.pred), "that is:", length(hitsOurs.pred) /length(ours.pred), "of full set"))
  # What is uORFome prediction score of subset of ours in theirs vs other part
  summary(all.scores[hitsOurs]);summary(all.scores[-hitsOurs])
  summary(gr$peptide_score[hitsTheirs]);summary(gr$peptide_score[-hitsTheirs])
  plot_comparison()
}

plot_comparison <- function() {
  ggplot(data.frame(all.scores = all.scores), aes(x = all.scores, fill = seq(length(grl)) %in% hitsOurs)) + geom_bar()
  #ratio from ours to theirs
  
  ratio <- 1675/645 # ours by theirs to get real number of hits
  ratio <- length(hitsOurs)
  theirSearchSpace <- length(theirs)
  ourSearchSpace <- length(ours)
  library(grid)
  grid.newpage()
  library(VennDiagram)
  ggplot <- draw.pairwise.venn(ourSearchSpace,
                               theirSearchSpace,
                               ratio, category = c("uORFome prediction", "McGillivray et al. prediction"),
                               lty = rep("blank", 2), fill = c("cyan", "red"),
                               alpha = rep(0.5, 2), cat.pos = c(0, 0),
                               cat.dist = rep(0.025, 2), cat.cex = c(1,1), cex = c(1,1,1))
  library(gridExtra)
  vennPred <- grid.arrange(gTree(children=ggplot), top=textGrob("Overlap between prediction pipelines", gp=gpar(fontsize=20,font=8)),
                           bottom="")
  
  # Start codon metrics
  uorfData <- getAllSequenceFeaturesTable()
  StartCodons <- uorfData$StartCodons
  table(StartCodons[hitsOurs])
  #table(StartCodons[uniqueOrder][finalCagePred])
  table(StartCodons[predicted])
  tab1 <- table(StartCodons[hitsOurs])/sum(table(StartCodons[hitsOurs]))
  #tab2 <- table(StartCodons$startCodon[uniqueOrder][finalCagePred])/sum(table(StartCodons$startCodon[uniqueOrder][finalCagePred]))
  tab2 <- table(StartCodons[predicted])/sum(table(StartCodons[predicted]))
  
  df <- data.frame(value = c(tab1, tab2), variable =c(names(tab1), names(tab2)),
                   pred  = c(rep("McGillivray et al. prediction", length(tab1)), rep("uORFome prediction", length(tab2))))
  
  cstarts <- ggplot(df, aes(x=variable,y=value,fill=factor(pred)))+
    geom_bar(stat="identity",position="dodge")+
    scale_fill_discrete(name="Prediction pipeline")+
    xlab("Start codon")+ylab("percentage")
  
  library(cowplot)
  plot_grid(vennPred,cstarts, align='hv',nrow=2,labels=c('A','B'))
}