source("./HelperLibraries.R")
source("./uORFome/helperScripts.R")
idFiles = list.files(idFolder)
cageFiles = list.files(cageFolder)
cageFiles <- cageFiles[grep(cageFiles, pattern = ".bed")]

source("./uORFome/PipelineParts.R")