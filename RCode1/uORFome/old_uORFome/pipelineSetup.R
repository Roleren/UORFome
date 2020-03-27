library(ORFikPipeline)
#source("./HelperLibraries.R")
source("./uORFome/helperScripts.R")
idFiles <- list.files(idFolder)
cageFiles <- grep(list.files(cageFolder), pattern = ".bed", value = TRUE)

source("./uORFome/PipelineParts.R")
