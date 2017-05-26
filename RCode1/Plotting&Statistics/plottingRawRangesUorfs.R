library(ggplot2)
library(reshape2)

###Output plot to inputfolder
###plot contains te values against eachother, normalizations against eachother
###Lengths of uorfs and the 
# plottingRawRangesUorfs = function(matrix,plottingName = "resultPlots.pdf"){
#   
#   
# }

plottingName = "resultPlots.pdf"
cat("plotting results\n")

#create list of counts
setwd("/export/valenfs/projects/uORFome/test_results/rangesOfUORFs/")
filenames <- list.files(pattern="*.rdata", full.names=TRUE)

ORFcounts = sapply(filenames, function(x){ 
  load(x,envir = .GlobalEnv)
  return(cbind(length(rangesOfuORFs),x)) 
  })
ORFcounts = t(ORFcounts)

ORFcounts1 = as.data.frame(matrix(ORFcounts[,1],ncol = 1))
colnames(ORFcounts1) = "y"
setwd("/export/valenfs/projects/uORFome/test_results/Plotting/")

pdf(plottingName)

####Histogram of UORF counts
# lengthPlot = ggplot(data = ORFcounts1,aes(x = 1:nrow(ORFcounts1))) + 
#   labs("counts of orfs") + xlab("sample") + ylab("counts") + geom_histogram(fill="orange")
lengthPlot = hist(x = unlist(ORFcounts1))
print(lengthPlot)

dev.off()