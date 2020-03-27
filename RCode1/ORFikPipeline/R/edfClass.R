#' Read standard experiment on server
#' @param relPath name of experiment, run list.experiments() for candidates
#' @inheritParams list.experiments
#' @param expInVarName keep experiment name in variable names? defualt (FALSE)
#' @export
read.experimentl <- function(relPath, dir = "/export/valenfs/data/processed_data/experiment_tables_for_R/",
                             expInVarName = FALSE) {
  if (tools::file_ext(relPath) == "") {
    df <- read.experiment(p(dir, relPath, ".csv"))
  } else if (tools::file_ext(relPath) == "csv"){
    df <- read.experiment(p(dir, relPath))
  } else stop("invalid file type of file")
  df@expInVarName <- expInVarName
  return(df)
}

#' Create standard experiment on server
#' Default zebrafish annotation
#' @inheritParams ORFik::create.experiment
#' @param txdb ("/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db")
#' @param fa a fasta index path ("/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db")
#' @param saveDir directory for ORFik experiments: default ("/export/valenfs/data/processed_data/experiment_tables_for_R/")
#' @export
create.experimentl <- function(dir, exper,
                               saveDir = "/export/valenfs/data/processed_data/experiment_tables_for_R/", types = c("bam", "bed", "wig"),
                               txdb = "/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db",
                               fa = "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa",
                               viewTemplate = TRUE) {
  create.experiment(dir, exper, saveDir, types, txdb, fa, viewTemplate)
}

#' Save template experiment to disc
#' @param temp a temperary ORFik experiment template
#' @param saveDir directory for ORFik experiments: default ("/export/valenfs/data/processed_data/experiment_tables_for_R/")
#' @export
save.experimentl <- function(temp, saveDir = "/export/valenfs/data/processed_data/experiment_tables_for_R/") {

  save.experiment(temp, file = p(saveDir, temp[1,2], ".csv"))
}

#' List current experiment available
#'
#' Will only search .csv extension, also exclude any experiment with the word template.
#' @param dir directory for ORFik experiments: default ("/export/valenfs/data/processed_data/experiment_tables_for_R/")
#' @param pattern allowed patterns in experiment file name: default ("*", all)
#' @param libtypeExclusive search for experiments with exclusivly this libtype, default (NULL, all)
#' @export
list.experiments <- function(dir = "/export/valenfs/data/processed_data/experiment_tables_for_R/",
                             pattern = "*", libtypeExclusive = NULL) {
  experiments <- list.files(path = dir, pattern = "\\.csv")
  experiments <- grep(experiments, pattern = pattern, value = TRUE)
  experiments <- experiments[grep(experiments, pattern = "template", value = FALSE, invert = TRUE)]
  es <- lapply(experiments, function(x) {
    e <- read.experimentl(x)
    return(e)
  })
  libtypes <- lapply(es, function(e) {
    return(unique(e$libtype))
    })
  runs <- lapply(es, function(e) {
    return(length(e$libtype))
  })
  dt <- data.table(name = gsub(".csv", "", experiments), libtypes, samples = runs)
  if (!is.null(libtypeExclusive)) {
    message(paste("subset on libtype:", libtypeExclusive))
    dt <- dt[libtypes %in% libtypeExclusive,]
  }
  return(dt)
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
