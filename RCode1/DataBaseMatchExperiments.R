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
  
  teSumPerUORF = rowSums(filtered[,2:ncol(filtered), with = F])
  uorfsWithValidTE = sum(teSumPerUORF > 0)
  ratioOfValidTEs =  uorfsWithValidTE/nrow(filtered)
}



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
}

createExperimentsTable = function(){
  experiments = listAllExperiments()
  
  insertTable(as.data.table(experiments),"experiments")
}