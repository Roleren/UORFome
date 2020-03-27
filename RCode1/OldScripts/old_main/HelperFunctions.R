
#as number coverter function
an = function(fac){
  return(as.numeric(as.character(fac)))
}

#' remove non finite in data.table
#'
#' Removes inf, NA, and NaN values
#' @param DT a data.table
#' @param replacement (0) what should the non finite values be changed to ?
removeNonFinite <- function(DT, replacement = 0){
  invisible(lapply(names(DT),function(.name)
    set(DT, which(!is.finite(DT[[.name]])), j = .name, value = replacement)))
  return(DT)
}

getRelativePathName = function(name){
  return (gsub(".*/", "", name))
}

#' Make directoy structure for orf finding
#'
#' The main Path is ./.. relative to RCode1/ location
orfikDirs <- function(mainPath, makeDatabase = F){
  setwd(mainPath)
  print(paste("main path for project will be: ", mainPath))
  resultsLoc <- resultsFolder
  if (!dir.exists(resultsLoc)) dir.create(resultsLoc)

  dir.create(p(resultsLoc,"/New_Cage_Leaders"))
  dir.create(p(resultsLoc,"/regionUORFs"))
  dir.create(p(resultsLoc,"/rangesOfUORFs"))
  dir.create(p(resultsLoc,"/fasta"))
  dir.create(p(resultsLoc,"/uorfIDs"))

  if (makeDatabase) {
    dir.create("dataBase")
    dir.create("dataBase/forests/")
    dir.create("dataBase/forests/predicateTables")
  }

  print("directories created successfully")
}

#Check if uorfRanges exist already, or if must be created.
###########Should make this more failsafe!!!!!!!!!! add possibility to give ranges!!!!!!!!
UorfRangesNotExists <- function(assignUorf = F, givenCage = NULL){
  if(!exists("rangesOfuORFs")){
    if(file.exists(getUORFRDataName(givenCage))){#!!!Will not work for single run now!!!
      if(assignUorf){
        cat("loading rangesOfuorf from folder\n",uorfFolder,"\n")
        load(getUORFRDataName(givenCage),envir = .GlobalEnv)
      }
      return(F)
    }else{return(T)}
  }
  return(F)
}

#' set up cluster for pipeline
#'
#' Cluster object saved as cl
#' @param maxCores what is max cores, if Null set to half of available cores
pipelineCluster <- function(maxCores = NULL, reset = FALSE, outfile = NULL){
  if(exists("cl") && is(cl, "cluster") && !reset){
    message("cluster already exists, abort if not correct")
  } else {
    if(exists("cl") && is(cl, "cluster") && reset) {
      setDefaultCluster(NULL)
      rm(cl)
    }
    library(doParallel)
    if(is.null(maxCores)){
      maxCores = as.integer(detectCores()-(detectCores()/2)) # using half
    }
    if (is.null(outfile)) {
      cl <- makeCluster(maxCores, type = "PSOCK")
    } else {
      cl <- makeCluster(maxCores, outfile = outfile, type = "PSOCK")
    }

    registerDoParallel(cl)
    assign("cl",cl,envir = .GlobalEnv)
  }

  message("running with number of threads: ", maxCores)
}

writeFasta <- function(input, file = p(mainFolder,"/fasta.fasta")){
  if(is(input, "GRangesList") | is(input, "GRanges")){ getSequencesFromFasta(input)
  } else {
    seqs <- input
  }
  writeXStringSet(seqs, filepath = file)
}

updateORFik <- function(branch = "master", user = "Roleren")  {
  devtools::install_github(paste0(user, "/ORFik"), ref = branch)
}

GRToWig <- function(gr, outputPrefix = p(mainFolder,"/test_")) {

  if(!(is(gr, "GRanges") | is(gr, "GAlignments"))){
    if(is.character(gr)){
      if (tools::file_ext(gr) == "bam") {
        gr <- GenomicAlignments::readGAlignments(gr)
      } else if(tools::file_ext(gr) == "bed") {
        gr <- fread.bed(gr)
      } else gr <- import(gr)
    } else stop("character or GR, nothing else supported!")
  }
  gr <- GRanges(gr)
  gr <- resize(gr, width = 1, fix = "start")
  neg <- strand(gr) == "-"
  forward <- gr[!neg]
  reverse <- gr[neg]
  forward <- GRanges(coverage(forward))
  forward <- resize(forward[forward$score > 0], width = 1, fix = "start")
  reverse <- GRanges(coverage(reverse))
  reverse <- resize(reverse[reverse$score > 0], width = 1, fix = "start")

  export.wig(forward, p(outputPrefix, "forward.wig"))
  export.wig(reverse, p(outputPrefix, "reverse.wig"))

}

# GRToWig_new <- function(gr, outputPrefix, shift = NULL) {
#
#   if(!(is(gr, "GRanges"))){
#     if(is.character(gr)){
#       gr <- ORFik:::fimport(gr)
#     } else stop("gr must either be character or GRanges!")
#   }
#   if (!is.null(shift)) {
#     if ("5prime" %in% shift) {
#
#     }
#   }
#
#   neg <- strand(gr) == "-"
#   export.wig(gr[!neg], p(outputPrefix, "_forward.wig"))
#   export.wig(gr[neg], p(outputPrefix, "_reverse.wig"))
#
# }



#' Extension of countOverlaps that can use $score column
countOverlapsScore <- function(query, subject, maxgap = -1L, minoverlap = 0L, type = "any") {
  if (is.null(subject$score)) {
    return(countOverlaps(query, subject, maxgap, minoverlap, type))
  }

  a <- findOverlaps(query, subject, maxgap, minoverlap, type)
  dt <- data.table(from = from(a), to = to(a), score = subject$score[to(a)])
  dt <- dt[, .(counts = sum(score)), by = from]
  res <- rep(0, length(query))
  res[dt$from] <- dt$counts
  names(res) <- names(query)
  return(res)
}

fpkmScore <- function(grl, reads, pseudoCount = 0)
{
  grl_len <- widthPerGroup(grl, FALSE)
  overlaps <- countOverlapsScore(grl, reads)
  librarySize <- length(reads)
  return(ORFik:::fpkm_calc(overlaps, grl_len, librarySize) + pseudoCount)
}

p <- paste0 # just for parsing relative paths together
