library(data.table)

source("./HelperLibraries.R")

idFiles = list.files(idFolder)
cageFiles = list.files(cageFolder)
cageFiles <- cageFiles[grep(cageFiles, pattern = ".bed")]

source("./PipelineParts.R")