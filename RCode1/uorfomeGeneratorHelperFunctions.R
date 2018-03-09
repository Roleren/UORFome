source("./HelperVariables.R")
source("./NameCreator.R")
source("./GenomicGetters.R")
source("./HelperLibraries.R")
source("./CageDataIntegration.R")
source("./scanUORFs.R")
source("./Plotting&Statistics/PlotUORFome.R")
source("./HelperFunctions.R")
source("./uOrfFeatures.R")
source("./GRangesHelpers.R")
source("./bedMaker.R")


#' check if uorf ranges already exists,
#' 
#' if not make load them or make them from scratch
decideHowToGetUORFRanges <- function(assignUorf = F,givenCage = NULL){ ###Watch out for assign uorf, might be buggy
  cat("started finding UORFS\n")
  if(UorfRangesNotExists(assignUorf,givenCage))
    rangesOfuORFs <- scanUORFs(fiveUTRs,saveToFile = T, assignUorf = assignUorf)
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
