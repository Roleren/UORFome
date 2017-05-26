library(ggplot2)
library(reshape2)

###Output plot to inputfolder
###plot contains te values against eachother, normalizations against eachother
###Lengths of uorfs and the 
plotting = function(matrix,plottingName = "resultPlots.pdf"){
  
  cat("plotting results\n")
  index2 = matrix[,"normUORFRNA"] < 0.1 | is.na(matrix[,"normUORFRNA"] ) | is.infinite(matrix[,"normUORFRNA"] ) | is.nan(matrix[,"normUORFRNA"] )
  matrix = matrix[!index2,]
  index3 = matrix[,"normUORFRFP"] < 0.1 | is.na(matrix[,"normUORFRFP"] ) | is.infinite(matrix[,"normUORFRFP"] ) | is.nan(matrix[,"normUORFRFP"] )
  matrix = matrix[!index3,]
  
  pdf(plottingName)
  maxUORFTEList = aggregate(as.numeric(as.matrix(matrix$teUORF)),by = list(matrix$tx_name),max)
  maxCDSTEList = aggregate(as.numeric(as.matrix(matrix$teCDS)),by = list(matrix$tx_name),max)
  ###plot of teuorfMax vs teCDS with 0 values for uorf removed
  
  maxTEUORFPLOT = log(maxUORFTEList$x)
  maxTEUORFPLOT[which(maxTEUORFPLOT == -Inf)] = NA
  maxTECDSPLOT = log(maxCDSTEList$x)
  plot(maxTEUORFPLOT,maxTECDSPLOT, main = "teCDS vs maxTEUORF, 0 values for uorf removed")
  abline(lm(maxTEUORFPLOT ~ maxTECDSPLOT))
  
  
  ###plot of teuorfmax vs teCDS with 0 values increased to min value
  minValueTEUORF = min(maxTEUORFPLOT, na.rm = T)
  maxTEUORFPLOT[which(is.na(maxTEUORFPLOT))] = minValueTEUORF
  plot(maxTEUORFPLOT,maxTECDSPLOT, main = "teCDS vs maxTEUORF, 0 values increased to min value")
  abline(lm(maxTEUORFPLOT ~ maxTECDSPLOT))
  
 
  
  ###plot normalizations against eachother
  CDSPLOT <- ggplot(matrix, aes(x = log(an(normCDSRNA)), y = log(an(normCDSRFP)))) + geom_point()
  print(CDSPLOT)
  fiveUTRPLOT <- ggplot(matrix, aes(x = log(an(norm5UTRRNA)), y = log(an(norm5UTRRFP)))) + geom_point()
  print(fiveUTRPLOT)
  threeUTRPLOT <- ggplot(matrix, aes(x = log(an(norm3UTRRNA)), y = log(an(norm3UTRRFP)))) + geom_point()
  print(fiveUTRPLOT)
  UORFPLOT <- ggplot(matrix, aes(x = log(an(normUORFRNA)), y = log(an(normUORFRFP)))) + geom_point()
  print(UORFPLOT)
  
  ####Histogram of UORF counts
  lengthPlot = ggplot(ORFLengths,aes(x = log10(x)) ) + 
      labs("counts of orfs") + xlab("log10 length of uorf") + ylab("number of counts") + geom_histogram(fill="orange")
  print(lengthPlot)
  
  #####Pie chart of amount of transcripts with uorfs
  pie(c(numberOfOverlaps,length(cds)),labels = c("with uorfs","without")
      ,main = paste0("Ratio of transcripts with uorfs\n Total:",length(cds)))
  ####BOXPLOT of uorf te's
  matrix$teCDS <- as.numeric(as.matrix(matrix$teCDS))
  matrix$te5UTR <- as.numeric(as.matrix(matrix$te5UTR))
  matrix$te3UTR <- as.numeric(as.matrix(matrix$te3UTR))
  matrix$teUORF <- as.numeric(as.matrix(matrix$teUORF))
  reshapedMat = melt(matrix,id.vars = "namesOFuorfs",measure.vars = c("teCDS","te5UTR","te3UTR","teUORF"),variable.name = "te")
  
  teUORFBOXplot <- ggplot(reshapedMat,aes(x = te ,y =log(value), fill = te)) + 
    xlab("Feature Type") + ylab("log values of TE") + geom_boxplot()
  print(teUORFBOXplot)
  dev.off()
}
