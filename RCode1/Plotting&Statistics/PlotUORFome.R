###Check what filter cutoff is a good one, linear or exponential
library(ggplot2)
library(reshape2)

###Output plot to inputfolder
###plot contains te values against eachother, normalizations against eachother
###Lengths of uorfs and the 
plotting = function(matrix,plottingName = "resultPlots.pdf"){
  
  cat("plotting results\n")
  #old filter
#   index2 = matrix[,"normUORFRNA"] < 0.1 | is.na(matrix[,"normUORFRNA"] ) | is.infinite(matrix[,"normUORFRNA"] ) | is.nan(matrix[,"normUORFRNA"] )
#   matrix = matrix[!index2,]
#   index3 = matrix[,"normUORFRFP"] < 0.1 | is.na(matrix[,"normUORFRFP"] ) | is.infinite(matrix[,"normUORFRFP"] ) | is.nan(matrix[,"normUORFRFP"] )
#   matrix = matrix[!index3,]
  ###Remove all bad values from teCDS
  #   index = matrix[,"teCDS"] == 0 | is.na(matrix[,"teCDS"] ) | is.infinite(matrix[,"teCDS"] ) | is.nan(matrix[,"teCDS"] )
  #   matrix = matrix[!index,]
  pdf(plottingName)
  maxUORFTEList = aggregate(as.numeric(as.matrix(matrix$teUORF)),by = list(matrix$tx_name),max)
  maxCDSTEList = aggregate(as.numeric(as.matrix(matrix$teCDS)),by = list(matrix$tx_name),max)
  ###plot of teuorfMax vs teCDS with 0 values for uorf removed
  
  maxTEUORFPLOT = log(maxUORFTEList$x)
  maxTEUORFPLOT[which(maxTEUORFPLOT == -Inf)] = NA
  maxTECDSPLOT = log(maxCDSTEList$x)
  plot(maxTEUORFPLOT,maxTECDSPLOT, main = "teCDS vs maxTEUORF, 0 values for uorf removed")
  #abline(lm(maxTEUORFPLOT ~ maxTECDSPLOT))
  
  
  ###plot of teuorfmax vs teCDS with 0 values increased to min value
  minValueTEUORF = min(maxTEUORFPLOT, na.rm = T)
  maxTEUORFPLOT[which(is.na(maxTEUORFPLOT))] = minValueTEUORF
  plot(maxTEUORFPLOT,maxTECDSPLOT, main = "teCDS vs maxTEUORF, 0 values increased to min value")
  #abline(lm(maxTEUORFPLOT ~ maxTECDSPLOT))
  
 
  
  ###plot normalizations against eachother
  CDSPLOT <- ggplot(matrix, aes(x = log(an(normCDSRNA)), y = log(an(normCDSRFP)))) + geom_point()
  print(CDSPLOT)
  fiveUTRPLOT <- ggplot(matrix, aes(x = log(an(norm5UTRRNA)), y = log(an(norm5UTRRFP)))) + geom_point()
  print(fiveUTRPLOT)
  threeUTRPLOT <- ggplot(matrix, aes(x = log(an(norm3UTRRNA)), y = log(an(norm3UTRRFP)))) + geom_point()
  print(fiveUTRPLOT)
  UORFPLOT <- ggplot(matrix, aes(x = log(an(normUORFRNA)), y = log(an(normUORFRFP)))) + geom_point()
  print(UORFPLOT)
  
  ####Histogram of frame proportions
  framePlot = ggplot(matrix,aes(x = frame,fill = pass_filter) ) + geom_bar(position = "dodge") +
    title(main = "Reading frame proportions uorf to cds start") + xlab("Frame 1,2,3") + ylab("counts of orfs per frame") 
  print(framePlot)
  
  #CHECK OVERLAPING DIFFERENT FRAME WHAT IS DIFFERENCE TE ?
#   > mean(na.exclude(matrix$teUORF[matrix$frame == 0 & !is.infinite(matrix$teUORF)]))
#   [1] 50.81061
#   > mean(na.exclude(matrix$teUORF[matrix$frame == 1 & !is.infinite(matrix$teUORF)]))
#   [1] 50.594
#   > mean(na.exclude(matrix$teUORF[matrix$frame == 2 & !is.infinite(matrix$teUORF)]))
#   [1] 55.56346
  
  ###plot width vs te value
  widthTEuorfPlot <- ggplot(matrix, aes(x = log(an(matrix$width)), y = log(an(teUORF)))) + geom_point()
  print(widthTEuorfPlot)
  ###plot quantiles of width, vs te use "cut" with quantiles
  
  ####Histogram of UORF counts
  
  lengthPlot = ggplot(matrix,aes(x = log10(width)) ) + 
      labs("counts of orfs") + xlab("log10 length of uorf") + ylab("number of counts") + geom_histogram(fill="orange")
  print(lengthPlot)
  
  #####Pie chart of amount of transcripts with uorfs
  if(exists("numberOfOverlaps") == F){ #This should be changed, so it uses the matrix!!!
    numberOfOverlaps = getUOrfOverlaps()
    assign("numberOfOverlaps",numberOfOverlaps,envir = .GlobalEnv)
  } 
  pie(c(numberOfOverlaps,length(cds)),labels = c("with uorfs","without")
      ,main = paste0("Ratio of transcripts with uorfs\n Total:",length(cds)))
  ####BOXPLOT of uorf te's
  matrix$teCDS <- as.numeric(as.matrix(matrix$teCDS)) ###TODO: fix this, put class def into helperfunc
  matrix$te5UTR <- as.numeric(as.matrix(matrix$te5UTR))
  matrix$te3UTR <- as.numeric(as.matrix(matrix$te3UTR))
  matrix$teUORF <- as.numeric(as.matrix(matrix$teUORF))
  reshapedMat = melt(matrix,id.vars = "uorfName",measure.vars = c("teCDS","te5UTR","te3UTR","teUORF"),variable.name = "te")
  
  teUORFBOXplot <- ggplot(reshapedMat,aes(x = te ,y =log(value), fill = te)) + 
    xlab("Feature Type") + ylab("log values of TE") + geom_boxplot()
  print(teUORFBOXplot)
  dev.off()
}
