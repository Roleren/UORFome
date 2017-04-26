#########################PLOTSCRIPT FOR UORFOME########################
require(ggplot2)
matrix = read.csv("matrix.csv")

createLogPlot(matrix[,"teRNA"],matrix[,"teUORF"],index4)
summary(predRNA)
createLogPlot(matrix[,"teUTR"],matrix[,"UORF"],index4)
summary(predUTR)

####################FUNCTIONS#####################

createLogPlot = function(xData,yData, removedXValues = NULL){
  
  if(!is.null(removedXValues)){
    xData = xData[!removedXValues]
  }
  xData = log(xData)
  yData = log(yData)
  
  df = data.frame(xData,yData)
  pred = predict(lm(yData ~ xData), df)
  
  ggplot(df,
         aes(x = xData, y = yData )
  ) + geom_point() +
    geom_line(aes(y = pred))
  
}
