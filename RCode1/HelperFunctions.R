
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
pipelineCluster <- function(maxCores = NULL, outfile = NULL){
  if(exists("cl") && any(class(cl) == "cluster")){
    message("cluster already exists, abort if not correct")
  } else {
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

updateORFik <- function(branch = "master", user = "JokingHero")  {
  devtools::install_github(paste0(user, "/ORFik"), ref = branch)
}