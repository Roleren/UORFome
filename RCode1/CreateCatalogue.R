p = function(nameA,nameB){
  paste0(nameA,nameB)
}

mainFolder = "./.."
resultsFolder = p(mainFolder,"/test_results") #output folder
matrixFolder = p(resultsFolder,"/Matrices/")
matrixFiles = list.files(matrixFolder)

# for(i in list.files(matrixFolder)){
#   load(p(matrixFolder,i))
#   #find a way to combine the matrices into a bigger format
# }
# 
# plottingFolder = p(resultsFolder,"/Plotting/Comparisons_plots/")
# #plotting(matrix,paste0(plottingFolder,gsub("%","_",detailedFullName),".pdf")) #plot results

exists(p(matrixFolder,matrixFiles[2]))

matrix = read.csv(p(matrixFolder,matrixFiles[5]))



