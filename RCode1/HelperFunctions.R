
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
    set(DT, which(!is.finite(DT[[.name]])), j = .name,value = replacement)))
  return(DT)
}

getRelativePathName = function(name){
  return (gsub(".*/", "", name))
}

#' fast import
#' 
#' versatile loader function, supports many formats
#' @param fileName the file path
fimport <- function(fileName){
  library(tools)
  if(!file.exists(fileName)) stop(paste(fileName, "does not name a file"))
  
  format <- file_ext(fileName)
  if(format == "") stop(paste("file", fileName,"have no format"))
  
  if (format == "csv") {
    return(fread(fileName))
  } else if (format == "tab") {
    return(fread(fileName, sep = "\t"))
  } else if (format == "bed") {
      return(ORFik:::cageFromFile(fileName))
  } else if (format == "bam") {
      return(readGAlignments(fileName))
  } else if (format == "rds") {
      return(readRDS(fileName))
  } else if (format == "rdata") {
    rdataname <- load(fileName)
    message(paste("file is loaded to global environment as:", rdataname))
    return(NULL)
  } else if (format == ".gtf") {
    return(makeTxDbFromGFF(fileName))
  } else {
    stop(paste("loading function for", fileName,
               "is not supported, do it manually."))
  }
}

#' Make directoy structure for orf finding
#' 
#' The main Path is ./.. relative to RCode1/ location
orfikDirs <- function(mainPath, makeDatabase = F){
  setwd(mainPath)
  print(paste("main path for project will be: ", mainPath))
  resultsLoc <- "test_results/"
  dir.create(resultsLoc)
  
  dir.create(p(resultsLoc,"New_Cage_Leaders"))
  dir.create(p(resultsLoc,"New_Cage_bedLeaders"))
  dir.create(p(resultsLoc,"rangesOfUORFs"))
  dir.create(p(resultsLoc,"bedUORFS"))
  dir.create(p(resultsLoc,"fasta"))
  
  if (makeDatabase) {
    dir.create(p(resultsLoc,"uorfIDs"))
    dir.create("dataBase")
  }
  
  print("directories created successfully")
}

#Check if uorfRanges exist already, or if must be created.
###########Should make this more failsafe!!!!!!!!!! add possibility to give ranges!!!!!!!!
UorfRangesNotExists <- function(assignUorf = F, givenCage = NULL){
  if(exists("rangesOfuORFs") == F){
    if(!assignUorf){ #if not loading or assigning to global, we need cage name
      thisCage <- givenCage
      assign("thisCage",thisCage,envir = .GlobalEnv)
    }
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
pipelineCluster <- function(maxCores = NULL, outfile = NULL){
  if(exists("cl") && any(class(cl) == "cluster")){
    message("cluster already exists, abort if not correct")
  } else {
    library(doParallel)
    if(is.null(maxCores)){
      maxCores = as.integer(detectCores()-(detectCores()/2)) # using half
    }
    if (is.null(outfile)) {
      cl <- makeCluster(maxCores)
    } else {
      cl <- makeCluster(maxCores, outfile = outfile)
    }
    
    registerDoParallel(cl)
    assign("cl",cl,envir = .GlobalEnv)
  }
  
  message("running with number of threads: ", maxCores)
}


loadMatrix = function(mname,toGlobalEnv = T){
  return(fread(input = p(matrixFolder,mname),header = T))
}

#' save image of current session, makes name automaticly
saveRData <- function(){
  if(exists("generalName") == F){
    save.image(p(RdataFolder,"results.Rdata"))
  }else{
    save.image(paste0(RdataFolder,detailedFullName,".results.Rdata",sep=""))
  }
}

loadRData = function(rdataname, toGlobalEnv = T){
  if(toGlobalEnv)
    load(p(RdataFolder,rdataname),envir = .GlobalEnv)
  else
    load(p(RdataFolder,rdataname))
}

loadUorfID = function(rdataname){
  load(paste0(resultsFolder,"/uorfIDs/",rdataname), envir = .GlobalEnv)
}

updateORFik <- function(branch = "master") {
  devtools::install_github("JokingHero/ORFik", ref = branch)
}