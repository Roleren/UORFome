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
    if(!tableNotExists(tableName))
      deleteTable(tableName)
  }
  dbWriteTable(uorfDB, tableName, as.data.table(Matrix),append = appends)
}

readTable = function(tableName, asGR = FALSE, with.IDs = TRUE){
  if (asGR){
    grl <- as.data.table(dbReadTable(uorfDB,tableName))
    return(makeGRangesListFromDataFrame(grl, split.field = "group",
                                        names.field = "group_name",
                                        keep.extra.columns = TRUE))
    
  } else{
    if (!with.IDs) {
      dt <- as.data.table(dbReadTable(uorfDB,tableName))
      if (!is.numeric(dt[1,1][[1]])) {
        dt <- dt[, -1]
        if (!is.numeric(dt[1,1][[1]])) {
          dt <- dt[, -1]
        }
      }
      return(dt)
    }
    return(as.data.table(dbReadTable(uorfDB,tableName)))
  }
}

listTables = function(){
  sort(dbListTables(uorfDB))
}

tableNotExists <- function(name){
  sum(grep(pattern = name, x = listTables())) == 0
}