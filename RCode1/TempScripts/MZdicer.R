getMir430 <- function(path = "/export/valenfs/projects/adam/TCP_seq/data_files/targetscanFish_180829_tidy.csv", 
                      filter = -0.20) {
  targets <- read.csv(path)
  targets <- targets[!(is.na(targets$miR430_score) | is.nan(targets$miR430_score)),]
  targets <- targets[targets$miR430_score <= filter,]
  genes <- targets$gene_id
  return(genes)
}


getmzDicerDf <- function() {
  RFPPath <- "/export/valenfs/data/processed_data/Ribo-seq/bazzini_2012_zebrafish/final_results/aligned_Zv9/merged/"
  RFP <- c(p(RFPPath, "ribo-Seq_MZdicer_2hpf.bam"), 
             p(RFPPath, "ribo-Seq_MZdicer_4hpf.bam"),
             p(RFPPath, "ribo-Seq_MZdicer_6hpf.bam"),
             p(RFPPath, "ribo-Seq_Wild_type_2hpf.bam"),
             p(RFPPath, "ribo-Seq_Wild_type_4hpf.bam"),
             p(RFPPath, "ribo-Seq_Wild_type_6hpf.bam"))
  RNAPath <- "/export/valenfs/data/processed_data/RNA-seq/bazzini_2012_zebrafish/final_results/aligned_Zv9/merged/"
  RNA <- c(p(RNAPath, "mRNA-Seq_MZdicer_2hpf.bam"), 
           p(RNAPath, "mRNA-Seq_MZdicer_4hpf.bam"),
           p(RNAPath, "mRNA-Seq_MZdicer_6hpf.bam"),
           p(RNAPath, "mRNA-Seq_Wild_type_2hpf.bam"),
           p(RNAPath, "mRNA-Seq_Wild_type_4hpf.bam"),
           p(RNAPath, "mRNA-Seq_Wild_type_6hpf.bam"))
  
  stage <- c(2, 4, 6, 2, 4, 6)
  type <- c(rep("MZ", 3), rep("WT", 3))
  df <- data.frame(RNA, RFP, stage, type)
  return(df)
}

bamVarName <- function(df, skip.replicate = FALSE) {
  libTypes <- libraryTypes(df)
  varName <- c()
  for (i in 1:nrow(df)) {
    for (lib in libTypes) {
      varName <- c(varName, ifelse(is.null(df$rep) | skip.replicate, 
                                   paste(lib, df$type[i], df$stage[i], sep = "_"),
                                   paste(lib, df$type[i], df$stage[i], paste0("r",df$rep[i]), sep = "_")))
    }
  }
  return(varName)
}

outputBams <- function(df, chrStyle = NULL) {
  validateExperiments(df)

  
  libTypes <- libraryTypes(df)
  varNames <- bamVarName(df)
  j <- 0
  for (i in 1:nrow(df)) { # For each stage
    print(i)
    for (lib in libTypes) { # For each library of that stage (SSU, LSU, RNA-seq, RIBO-seq)
      j <- j + 1
      if (exists(varNames[j])) next
      
      bam <- ORFik:::readBam(df[i, lib], chrStyle)
      assign(varNames[j], bam, envir = .GlobalEnv)
    }
  }
}

#' Try to optimze bam reader
#' 
#' Extract fields 3,4,6 and 7
#' 3: chromosome
#' 4: left most position, 
#' @param path a path to bam file
#' @param samtoolsPath If not samtools in global scope, give full path here to samtools binary
#' @return a GAlignment object of bamfile
readBamNew <- function(path, samtoolsPath = "/Home/ii/hakontj/bin/samtools-1.9/samtools") {
  
  a <- fread(cmd = paste(samtoolsPath, "view -@ 30", path, "| cut -f 3-4,6-7"), 
             strip.white = F, colClasses = c("character", "integer", "character", "character"))
  a[data.table::`%chin%`(V4, "="), V4 := "*"]
  return(GAlignments(Rle(factor(a$V1)), pos = as.integer(a$V2), cigar = a$V3,
                     strand = Rle(factor(a$V4, levels = c("+", "-", "*")))))
}

#' Try to optimze bam reader
#' 
#' Extract fields 3,4,6 and 7
#' 3: chromosome
#' 4: left most position, 
#' @param path a path to bam file
#' @param samtoolsPath If not samtools in global scope, give full path here to samtools binary
#' @return a GAlignment object of bamfile
# readBam <- function(path, chrStyle = NULL, sambambaPath = "~/bin/sambamba-0.7.0-linux-static") {
#   
#   a <- fread(cmd = paste(sambambaPath, "view -t 30", path, "| cut -f 3-4,6-7"), 
#              strip.white = F, colClasses = c("character", "integer", "character", "character"), 
#              nThread = 5)
#   a[data.table::`%chin%`(V4, "="), V4 := "*"]
#   return(GAlignments(Rle(factor(a$V1)), pos = as.integer(a$V2), cigar = a$V3,
#                      strand = Rle(factor(a$V4, levels = c("+", "-", "*")))))
# }

