#Col1 is healthy, Col2 is unhealthy, if there is such a difference
#Else the difference is specified
library(ggplot2)

data1 = "/export/valenfs/projects/uORFome/Gonzales brain-tissue/test_data2/matrix.csv"
data2 = "/export/valenfs/projects/uORFome/Gonzales brain-tissue/test_data4/matrix.csv"
lab1 = "healthy"
lab2 = "tumor"
matrix1 = read.csv(data1)
matrix2 = read.csv(data2)

#conbine them
teuORFs = cbind(matrix1$teUORF,matrix2$teUORF)
#fix bad values
index2 = teuORFs < 0.1 | is.na(teuORFs) | is.infinite(teuORFs ) | is.nan(teuORFs)
teuORFs[index2] = 0 

teuORFs = as.data.frame(teuORFs)


#finding quantiles
quant = quantile(teuORFs$V1, .95)

# teuORFs$quant  <- with(teuORFs, factor(ifelse(value < quants[1], 0, 
#                                       ifelse(value < quants[2], 1, 2))))

teuORFs$quant  <- with(teuORFs, factor(ifelse(V1 < quant[1],0  ,1)))

ggplot(teuORFs,aes(x = log(V1),y = log(V2) )) + geom_point(aes(color = quant)) + xlab(lab1) + ylab(lab2) + scale_colour_manual(values =  c("black", "red"))
