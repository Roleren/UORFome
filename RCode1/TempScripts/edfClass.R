experiment <- setClass("experiment", 
                       slots=list(experiment = "character",txdb = "character",
                                  fafile = "character", expInVarName = "logical"), 
                       contains = "DataFrame")

setMethod("show",
          "experiment",
          function(object) {
            cat("experiment:", object@experiment, "with", 
                length(unique(object@listData$libtype)), "library types and",
                length(object@listData$libtype), "runs","\n")
            show(as.data.table(as(object@listData, Class = "DataFrame"))
                 [,-"filepath"])
          }
)

setMethod("nrow", 
          "experiment",
          function(x) {
  nrow(as.data.table(as(x@listData, Class = "DataFrame")))
  }
)

#' Read ORFik experiment
#' 
#' @param file a .csv file following ORFik experiment style.
read.experiment <-  function(file) {
  info <- read.table(file, sep = ",", nrows = 3, stringsAsFactors = FALSE)
  exper <- info[1,2]
  txdb <- ifelse(is.na(info[2,2]),  "", info[2,2])
  fa <- ifelse(is.na(info[3,2]),  "", info[3,2])
  listData <- read.csv2(file, skip = 3, header = T, sep = ",",
                        stringsAsFactors = FALSE)
  df <- experiment(experiment = exper, txdb = txdb, fafile = fa,
                   listData = listData, expInVarName = TRUE)
  
  validateExperiments(df)
  return(df)
}

#' Which type of experiments?
#' @param df an ORFik experiment data.frame
#' @return NULL
libraryTypes <- function(df){
  if (is(df, "experiment")) {
    return(unique(df$libtype))
  } else if (is(df, "character") | is(df, "factor")) {
    return(gsub("_.*", x = df, replacement = ""))
  } else stop("library types must be data.frame or character vector!")
}

#' Validate ORFik experiment
#' @param df an ORFik experiment data.frame
#' @return NULL
validateExperiments <- function(df) {
  libTypes <- libraryTypes(df)
  if (!is(df, "experiment")) stop("df must be experiment!")
  if (!all((c("stage", "libtype") %in% colnames(df))))
    stop("stage, libtype and experiment must be colnames in df!")
  if (length(libTypes) == 0) stop("df have no valid sequencing libraries!")
  if (nrow(df) == 0) stop("df must have at least 1 row!")
  
  emptyFiles <- c()
  for (i in df$filepath) {
    emptyFiles <- c(emptyFiles, as.numeric(sapply(as.character(i),
                                                  file.size)) == 0)
  }
  if (any(is.na(emptyFiles)))
    stop(paste("File is not existing:\n",df$filepath[is.na(emptyFiles)]))
  if (any(emptyFiles)) {
    print(cbind(df[which(emptyFiles),]))
    stop("Empty files in list, see above for which")
  }
  if (length(bamVarName(df)) != length(unique(bamVarName(df))))
    stop("experiment table has non-unique rows!")
}

#' Get variable names from experiment
#' @param df an ORFik experiment data.frame
#' @param skip.replicate a logical (FALSE), don't include replicate
#' in variable name.
#' @param skip.condition a logical (FALSE), don't include condition
#' in variable name.
#' @param skip.stage a logical (FALSE), don't include stage
#' in variable name.
#' @return NULL
bamVarName <- function(df, skip.replicate = length(unique(df$rep)) == 1,
                       skip.condition = length(unique(df$condition)) == 1,
                       skip.stage = length(unique(df$stage)) == 1, 
                       skip.experiment = !df@expInVarName) {
  libTypes <- libraryTypes(df)
  varName <- c()
  for (i in 1:nrow(df)) {
    varName <- c(varName, bamVarNamePicker(df[i,], skip.replicate, 
                                           skip.condition, skip.stage, 
                                           skip.experiment))
  }
  return(varName)
}

#' Get variable names from experiment
#' @param df an ORFik experiment data.frame
#' @param skip.replicate a logical (FALSE), don't include replicate
#' in variable name.
#' @param skip.condition a logical (FALSE), don't include condition
#' in variable name.
#' @param skip.stage a logical (FALSE), don't include stage
#' in variable name.
#' @return NULL
bamVarNamePicker <- function(df, skip.replicate = FALSE, skip.condition = FALSE,
                             skip.stage = FALSE, skip.experiment = FALSE) {
  if(nrow(df) != 1) stop("experiment must only input 1 row")
  lib <- df$libtype
  stage <- df$stage
  cond <- df$condition
  rep <- df$rep
  current <- lib
  if(!skip.condition)
    current <- paste(current, cond, sep = "_")
  if (!skip.stage)
    current <- paste(current, stage, sep = "_")
  if (!(skip.replicate | is.null(rep)))
    current <- paste(current, paste0("r", rep), sep = "_")
  if (! (skip.experiment | is.null(df@experiment)))  
    current <- paste(df@experiment, current, sep = "_")
  return(current)
}

#' Output bam/bed/wig files to R as variables
#' 
#' Variable names defined by df
#' @param df an ORFik experiment data.frame
#' @param chrStyle the sequencelevels style (GRanges object or chr)
#' @return NULL
outputLibs <- function(df, chrStyle = NULL) {
  validateExperiments(df)
  
  
  libTypes <- libraryTypes(df)
  varNames <- bamVarName(df)
  for (i in 1:nrow(df)) { # For each stage
    print(i)
    if (exists(varNames[i])) next
    reads <- ORFik:::fimport(df[i,]$filepath, chrStyle)
    assign(varNames[i], reads, envir = .GlobalEnv)
  }
}

#' Get all bam files in folder
#' @param dir The directory to find .bam files.
#' @return (character vector) all .bam files
allBamFilesInFolder <- function(dir) {
  files <- grep(pattern = ".bam", x = list.files(dir, full.names = T), value = T)
  bai <- -grep(pattern = ".bai", x = files)
  if (length(bai)) {
    return(files[bai])
  }
  return(files)
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