#Col1 is healthy, Col2 is unhealthy, if there is such a difference
#Else the difference is specified
library(ggplot2)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(Biostrings)
setwd("/export/valenfs/projects/uORFome/test_results/comparisons/GonzalesHBrainAdult1/")
data1 = "/export/valenfs/projects/uORFome/test_results/comparisons/GonzalesHBrainAdult1/testWithoutCage/matrix.csv"
data2 = "/export/valenfs/projects/uORFome/test_results/comparisons/GonzalesHBrainAdult1/testWithCage/matrix.csv"
load("/export/valenfs/projects/uORFome/test_results/comparisons/GonzalesHBrainAdult1/testWithoutCage/rangesOfuorfs.rdata")
rangesOfuORFs1 = rangesOfuORFs
rm(rangesOfuORFs)
load("/export/valenfs/projects/uORFome/test_results/comparisons/GonzalesHBrainAdult1/testWithCage/rangesOfuorfs.rdata")
rangesOfuORFs2 = rangesOfuORFs
rm(rangesOfuORFs)

lab1 = "teCDS D/H"
lab2 = "teUORF D/H"

matrix1 = read.csv(data1)
matrix2 = read.csv(data2)

#combine them
# teuORFs = cbind(matrix1$normUORFRFP,matrix2$normUORFRFP)
# teCDSs = cbind(matrix1$normCDSRFP,matrix2$normCDSRFP)

ag1 = aggregate( .~tx_name  ,data = matrix1,FUN = sum)
ag2 = aggregate( .~tx_name  ,data = matrix2,FUN = sum)


totalMat = merge(ag1, ag2, by.x = "Group.1", by.y = "Group.1", all = F)

grl1 = as.data.frame(rangesOfuORFs1)
grl2 = as.data.frame(rangesOfuORFs2)
ORFLengths1 = aggregate(grl1[,"width"],by=list(grl1[,"names"]),sum)
ORFLengths2 = aggregate(grl1[,"width"],by=list(grl2[,"names"]),sum)
#In orflength, find overlaps by transcripname.

rename1 = gsub("_[0-9]*","", ORFLengths1$Group.1)

vector = totalMat$tx_len.x
names(vector) = totalMat$tx_name.x
index = match(rename1,names(vector))
vector1 = vector[index]

rename2 = gsub("_[0-9]*","", ORFLengths2$Group.1)
vector = totalMat$tx_len.y
names(vector) = totalMat$tx_name.y
index = match(rename2,names(vector))
vector2 = vector[index]
numberOfUTRSCage = length(vector2)

index2 = totalMat[,"teUORF.y"] < 1 | is.na(totalMat[,"teUORF.y"]) | is.infinite(totalMat[,"teUORF.y"]) | is.nan(totalMat[,"teUORF.y"])
totalMat = totalMat[!(index2[,1] | index2[,2]),]

##ggplot(totalMat,aes(x = log(/V1),y = log(V4/V3) )) + geom_point() + xlab(lab1) + ylab(lab2)
pdf("comparison.pdf")
lengthPlot = ggplot() +
   geom_density(data = ORFLengths1,aes(x = log10(x)) ,color="orange") +
   geom_density(data = ORFLengths2,aes(x = log10(x)) ,color="red") +
   xlab("log10 length of uorf") + ylab("number of counts") 
print(lengthPlot)
tePlot = ggplot() + 
  geom_point(data = totalMat,aes(x = log(teUORF.x),color="red")) + 
  geom_point(data = totalMat,aes(x = log(teUORF.y,color="blue"))) + 
  xlab("teUORF non cage") + ylab("teUORF  cage") 
print(tePlot)
teCDSPlot = ggplot() + 
  geom_point(data = totalMat,aes(x = log(teCDS.x),color="red")) + 
  geom_point(data = totalMat,aes(x = log(teCDS.y,color="blue"))) + 
  xlab("teUORF non cage") + ylab("teUORF  cage") 
print(teCDSPlot)

###Show summary
###Per transcript
dev.off()

###Compare length histogram
###te uorfs 

