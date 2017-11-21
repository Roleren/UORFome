library(DBI)
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./CreateCatalogueHelpers.R")
setwd("/export/valenfs/projects/uORFome/dataBase/")
#databaseFolder = "../dataBase/"
databaseName = "uorfCatalogue"
databaseName = paste0(databaseName,".sqlite")
name = databaseName


createDataBase = function(name){
  return (dbConnect(RSQLite::SQLite(), name))
}

deleteDataBase = function(name){
  dbDisconnect(uorfDB)
  unlink(name)
}
deleteTable = function(tableName){
  dbRemoveTable(uorfDB,tableName)
}

insertTable= function(Matrix, tableName, appends = F){
  dbWriteTable(uorfDB, tableName, Matrix,append = appends)
}
innerJoinInsertTableTE = function(Matrix, table1, table2){
  dbSendQuery(uorfDB,paste('SELECT',paste0(table1,'.uorfID,'),paste0(table1,'teUORF,'),paste0(table2,'teUORF'), 'FROM', table1,'INNER JOIN', table2,'ON', "teUORF" ,'=', "teUORF"))
}


readTable = function(tableName){
  return(as.data.table(dbReadTable(uorfDB,tableName)))
}

createCatalogueDB = function(name, bigMatrix, tableName){
  uorfDB = createDataBase(name)
  insertTable(bigMatrix,tableName)
  
  dbListTables(uorfDB)
  #test, get first row
  dbGetQuery(uorfDB, paste('SELECT * FROM', tableName ,'ORDER BY ROWID ASC LIMIT 1') )
  
  
  deleteDataBase(name)
  rm(uorfDB)
}
uorfDB = createDataBase(databaseName)
experiments = listAllExperiments()

# insertTable(as.data.table(experiments),"experiments")

# i = matrixFiles[1]
# i = gsub(".csv","",i)
# i = gsub("%",".",i)
# i = gsub("\\.","",i)
# tableName = gsub("-","_",i)
#nameOfBigMatrix = i

createUORFAtlas = function(){
  j = 1
  for(i in idFiles[1735:length(idFiles)]){
    load(p(idFolder, i))
    
    tes = uorfID[!duplicated(uorfID)]
    if(j == 1) {
      tesTotal = data.table(cbind(tes,rep(T,length(tes))))
      colnames(tesTotal) = c("uorfID", "1")
    }else{
      tes = data.table(cbind(tes,rep(T,length(tes))))
      colnames(tes) = c("uorfID", toString(j))
      tesTotal = merge(tesTotal,tes,by = "uorfID", all = T)
      colnames(tesTotal) = c("uorfID", as.character(1:(ncol(tesTotal)-1)))
    }
    j = j+1
    
  }
  #save(tesTotal,file = "tempUORFAtlas.rdata")
  load("./tempUORFAtlas.rdata")
  ncol(tesTotal)
  nrow(tesTotal)
  insertTable(Matrix = tesTotal,tableName = "uorfAtlas")
  rm(list=ls())
}

createTeTable = function(){
  teTable = "TeValues"
  j = 1
  for(i in matrixFiles){
    matrix = loadMatrix(i)
    tes = matrix[ , .(uorfID,teUORF)]
    tes = tes[!duplicated(tes$uorfID)]
    if(j == 1) {
      tesTotal = tes
      colnames(tesTotal) = c("uorfID", "1")
    }else{
      tesTotal = merge(tesTotal,tes,by = "uorfID")
      colnames(tesTotal) = c("uorfID", as.character(1:(ncol(tesTotal)-1)))
    }
    j = j+1
    
  }
  
  insertTable(Matrix = tesTotal,tableName = teTable)
  
  teFilteredTable = "filteredTeValues"
  filtered = tesTotal
  for(i in 2:ncol(tesTotal)){
    tempAnswer = tesTotal[,i, with =F]
    tempIndexes = is.na(tempAnswer) | is.infinite(unlist(tempAnswer))
    tempAnswer[tempIndexes] = 0
    filtered[,i] = tempAnswer
  }
  insertTable(Matrix = filtered,tableName = teFilteredTable)
}

teSumPerUORF = rowSums(filtered[,2:ncol(filtered), with = F])
uorfsWithValidTE = sum(teSumPerUORF > 0)
ratioOfValidTEs =  uorfsWithValidTE/nrow(filtered)

createPassFilter = function(){
  pfTable = "pfTable"
  j = 1
  for(i in matrixFiles){
    matrix = loadMatrix(i)
    tes = matrix[ , .(uorfID,pass_filter)]
    tes = tes[!duplicated(tes$uorfID)]
    if(j == 1) {
      tesTotal = tes
      colnames(tesTotal) = c("uorfID", "1")
    }else{
      tesTotal = merge(tesTotal,tes,by = "uorfID")
      colnames(tesTotal) = c("uorfID", as.character(1:(ncol(tesTotal)-1)))
    }
    j = j+1
    
  }
  insertTable(Matrix = tesTotal,tableName = pfTable)
}

doubleFiltered = filtered
for(i in 2:ncol(tesTotal)){
  tempAnswer = doubleFiltered[,i, with =F]
  tempIndexes = !tesTotal[,i,with=F]
  tempAnswer[tempIndexes] = 0
  doubleFiltered[,i] = tempAnswer
}

fteSumPerUORF = rowSums(doubleFiltered[,2:ncol(doubleFiltered), with = F])
fuorfsWithValidTE = sum(teSumPerUORF > 0)
fratioOfValidTEs =  uorfsWithValidTE/nrow(doubleFiltered)

l = data.table(1:nrow(doubleFiltered))
for(i in 1:nrow(doubleFiltered)){
  l[i,]= min(doubleFiltered[i,2:ncol(doubleFiltered),with=F])
  
}
insertTable(Matrix = l,tableName = "teRowSums")
#createCatalogueDB(databaseName,matrix,tableName)
#0st name of experiments
###1st cage supprt per uorf, bool, 
###2nd rna ribo seq supprt and te per uorf, matrix
###3rd 

###example: find uorf that are only in certain tissue