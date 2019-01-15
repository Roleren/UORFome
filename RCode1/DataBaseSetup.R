library(DBI)

setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./pipelineSetup.R")

source("./databaseHelpers.R")
source("./DataBaseGetters.R")
source("./DataBaseInfo.R")
source("./DataBaseValidation.R")
source("./DataBaseGroupers.R")
source("./DataBaseAtlasFunctions.R")
source("./DataBaseCreator.R")

source("./features.R")
source("./Classifier.R")
source("./ClassifierHelpers.R")
source("./TissueTables.R")
source("./experiment.R")
source("./bedMaker.R")
source("./Plotting&Statistics/PlotUORFome.R")

dataBaseFolder <- p(mainFolder,"/dataBase")
if(!dir.exists(dataBaseFolder)){
  stop("data base folder not found")
}

setwd(dataBaseFolder)
#databaseName = "final_results/uorfCatalogue"
databaseName = "uorfCatalogue"
databaseName = p(databaseName,".sqlite")

uorfDB <- createDataBase(databaseName)
