createDataBase = function(name){
  return (dbConnect(RSQLite::SQLite(), name))
}

deleteDataBase = function(name){
  dbDisconnect(uorfDB)
  unlink(name)
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

      return(removeIDColumns(dt))
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

deleteTable = function(tableName){
  dbRemoveTable(uorfDB,tableName)
}

#' Delete uorf tables
#' 
#' For rerunning
deleteUorfTables <- function() {
  deleteTable("uorfsAsGRWithTx")
  deleteTable("uniqueIDs")
  deleteTable("SplittedByExonsuniqueUORFs")
  
  deleteTable("RSS")
  deleteTable("RRS")
  deleteTable("Ribofpkm")
  deleteTable("ORFScores")
  deleteTable("RiboByTissueMean")
  deleteTable("RiboByTissueTF")
  deleteTable("TEByTissueMeanWithInf")
  deleteTable("TEByTissueMeanWithoutInf")
  deleteTable("TEByTissueTF")
  deleteTable("teFiltered")
  deleteTable("teUnfiltered")
  deleteTable("tissueAtlasByCage")
  
  deleteTable("numberOfUorfsPerTx")
  deleteTable("rankInTx")
  deleteTable("disengagementScores")
  
  deleteTable("entropyRFP")
  deleteTable("fractionLengths")
  deleteTable("ioScore")
  deleteTable("kozak")
  deleteTable("inFrameCDS")
  deleteTable("isOverlappingCds")
  
  deleteTable("distORFCDS")
  deleteTable("linkORFsToTx")
  deleteTable("isOverlappingCds")
  deleteTable("floss")
  
  
  if(!tableNotExists("biggestTEVariance")) {
    deleteTable("biggestTEVariance")
    deleteTable("smallestTEVariance")
  }
}