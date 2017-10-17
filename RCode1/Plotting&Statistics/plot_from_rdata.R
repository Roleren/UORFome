p = function(nameA,nameB){
  paste0(nameA,nameB)
}

mainFolder = "./.."
resultsFolder = p(mainFolder,"/test_results") #output folder
RdataFolder = p(resultsFolder,"/Rdata/")
for(i in list.files(RdataFolder)){
  load(p(RdataFolder,i))
  plottingFolder = p(resultsFolder,"/Plotting/Single_result_Plots/")
  plotting(matrix,paste0(plottingFolder,gsub("%","_",detailedFullName),".pdf")) #plot results
}
