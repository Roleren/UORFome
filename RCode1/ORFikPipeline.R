rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/") 
source("./HelperLibraries.R")
source("./SummarizedExperimentHelpers.R")
source("./tcp_pipeline.R")
source("./edfClass.R")
source("./edfLists.R")
source("./transcriptHelpers.R")
source("./slackr_support.R")

p <- paste0 # just for parsing relative paths together

#' Save a reduced version of bam
#' 
# save.quickBam <- function(bam = "/export/valenfs/data/processed_data/TCP-seq/valen_all_withrRNA/aligned/17-Shield13_S3_R1_001_Aligned.sortedByCoord.out.bam", 
#                       outdir = "/export/valenfs/data/processed_data/") {
#   system.time(a <- readBam(bam))
#   nameShort <- basename(bam)
#   nameShort <- gsub(pattern = ".bam", replacement = ".quickBam", nameShort)
#   longName <- ORFik:::pasteDir(outdir, nameShort)
#   fwrite(x = as.data.table(a)[,.(seqnames, cigar, start, strand)], file = longName, compress = "gzip")
#   return(NULL)
# }

#' Load a reduced version of bam
#' 
#' Must be of format .quickBam
#' Gzipped compressed table of 4 columns from bam file:
#' seqnames, cigar, start and strand
# read.quickBam <- function(bam) {
#   if (tools::file_ext(bam) != "quickBam") stop("Format must be .quickBam")
#   
#   system.time(A <- fread(cmd =  paste("gunzip -c", longName), header = TRUE, nThread = 4, 
#                          colClasses = c("character", "character", "integer", "character")))
#   system.time(B <- GAlignments(seqnames = Rle(as.factor(A$seqnames)), cigar = A$cigar,
#                                pos = A$start, strand = factor(A$strand, levels = c("+","-","*"))))
# }