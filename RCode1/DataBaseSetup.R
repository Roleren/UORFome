library(DBI)
library(data.table)


setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./uorfomeGeneratorHelperFunctions.R")
source("./databaseHelpers.R")
source("./DataBaseGetters.R")
source("./DataBaseInfo.R")
source("./DataBaseValidation.R")
source("./TissueTables.R")
source("./PipelineParts.R")

source("./HelperFunctions.R")
source("./HelperVariables.R")
source("./GRangesHelpers.R")
source("./DataBaseAtlasFunctions.R")
source("./DataBaseCreator.R")
source("./Classifier.R")
source("./ClassifierHelpers.R")

#goals:
#1. get some way to load matrices, they are big and many!
#2. how should we filter here, only passed uorfs ? yes, i think so
#3. combine the equal uorf ids from samme tissue, so tissue matrices
#4. check for something interesting, what could it be ?
#5. I have one idea, test te uorf healthy vs sick for each matrix

uorfFiles = list.files(uorfFolder)
idFiles = list.files(idFolder)
cageFiles = list.files(cageFolder)
cageFiles <- cageFiles[grep(cageFiles, pattern = "bed")]


dataBaseFolder <- "/export/valenfs/projects/uORFome/dataBase"
if(!dir.exists(dataBaseFolder)){
  stop("ribo-seq folder not found")
}

setwd(dataBaseFolder)
databaseName = "uorfCatalogue"
databaseName = p(databaseName,".sqlite")

uorfDB <- createDataBase(databaseName)
