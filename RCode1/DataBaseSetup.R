library(DBI)
library(data.table)


setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./HelperLibraries.R")

source("./databaseHelpers.R")
source("./DataBaseGetters.R")
source("./DataBaseInfo.R")
source("./DataBaseValidation.R")
source("./TissueTables.R")
source("./PipelineParts.R")
source("./DataBaseGroupers.R")
source("./DataBaseAtlasFunctions.R")
source("./DataBaseCreator.R")
source("./features.R")
source("./Classifier.R")
source("./ClassifierHelpers.R")
source("./experiment.R")

#goals:
#1. get some way to load matrices, they are big and many!
#2. how should we filter here, only passed uorfs ? yes, i think so
#3. combine the equal uorf ids from samme tissue, so tissue matrices
#4. check for something interesting, what could it be ?
#5. I have one idea, test te uorf healthy vs sick for each matrix

uorfFiles = list.files(uorfFolder)
idFiles = list.files(idFolder)
cageFiles = list.files(cageFolder)
cageFiles <- cageFiles[grep(cageFiles, pattern = ".bed")]


dataBaseFolder <- p(mainFolder,"/dataBase")
if(!dir.exists(dataBaseFolder)){
  stop("data base folder not found")
}

setwd(dataBaseFolder)
#databaseName = "final_results/uorfCatalogue"
databaseName = "uorfCatalogue"
databaseName = p(databaseName,".sqlite")

uorfDB <- createDataBase(databaseName)
