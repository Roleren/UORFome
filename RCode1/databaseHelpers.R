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

insertTable= function(Matrix, tableName, appends = F, rmOld = F){
  if (rmOld){
    deleteTable(tableName)
  }
  dbWriteTable(uorfDB, tableName, as.data.table(Matrix),append = appends)
}

readTable = function(tableName, asGR = F){
  if (asGR){
    grl <- as.data.table(dbReadTable(uorfDB,tableName))
    return(makeGRangesListFromDataFrame(grl, split.field = "group"))
    
  } else{
    return(as.data.table(dbReadTable(uorfDB,tableName)))
  }
}

listTables = function(){
  dbListTables(uorfDB)
}

tableNotExists <- function(name){
  sum(grep(pattern = name, x = listTables())) == 0
}