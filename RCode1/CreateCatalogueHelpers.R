library(DBI)
library(data.table)

source("./HelperFunctions.R")
source("./HelperVariables.R")
source("./GRangesHelpers.R")
source("./DataBaseMatchExperiments.R")
#goals:
#1. get some way to load matrices, they are big and many!
#2. how should we filter here, only passed uorfs ? yes, i think so
#3. combine the equal uorf ids from samme tissue, so tissue matrices
#4. check for something interesting, what could it be ?
#5. I have one idea, test te uorf healthy vs sick for each matrix


matrixFiles = list.files(matrixFolder)
uorfFiles = list.files(uorfFolder)
idFiles = list.files(idFolder)
cageFiles = list.files(cageFolder)



setwd("/export/valenfs/projects/uORFome/dataBase/")
databaseName = "uorfCatalogue"
databaseName = paste0(databaseName,".sqlite")
name = databaseName

# for(i in list.files(matrixFolder)){
#   load(p(matrixFolder,i))
#   #find a way to combine the matrices into a bigger format
# }
# 
# plottingFolder = p(resultsFolder,"/Plotting/Comparisons_plots/")
# #plotting(matrix,paste0(plottingFolder,gsub("%","_",detailedFullName),".pdf")) #plot results

#'Group the matricies by tissue
groupMatriciesByTissue = function(){
  
}

isHealthy = function(matrixName){
  if(length(grep(pattern = "cancer|PC3|HEK|tumor|sick",x = matrixName,ignore.case = T)) > 0)
    return(F)
  
  return(T) 
}
#'split specific tissue into healthy and sick
splitIntoHealthyAndSick = function(){
  healthMatrix = Bigmatrix[healthy == T]
  assign("healthMatrix",healthMatrix,envir = .GlobalEnv)
  sickMatrix = Bigmatrix[healthy == F]
  assign("sickMatrix",sickMatrix,envir = .GlobalEnv)
}

getExperimentName = function(name){
  i = gsub(".csv","",name)
  i = gsub("%",".",i)
  i = gsub("\\.","",i)
  return(gsub("-","_",i))
}

listAllExperiments = function(){
  l = rep(NA,length(matrixFiles))
  j = 1
  for(i in matrixFiles){
    l[j] = getExperimentName(i)
    j = j+1
  }
  return(l)
}
