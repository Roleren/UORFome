setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./CreateCatalogueHelpers.R")




#'Load all matricies and store them some how
loadAllMatrices = function(){
  Bigmatrix = data.table()
  for(i in matrixFiles){
    
    matrix = loadMatrix(i)
    matrix = matrix[matrix$pass_filter == T,]
    if(isHealthy(i)){
      healthy = rep(T,nrow(matrix))
    }else{
      healthy = rep(F,nrow(matrix))
    }
    matrix$healthy = healthy
    Bigmatrix = rbindlist(l = list(Bigmatrix,matrix))
    rm(matrix)
    #find a way to combine the matrices into a bigger format
  }
  assign("Bigmatrix",Bigmatrix,envir = .GlobalEnv)
}


loadAllMatrices()
splitIntoHealthyAndSick()
