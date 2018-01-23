#Col1 is healthy, Col2 is unhealthy, if there is such a difference
#Else the difference is specified
library(ggplot2)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(Biostrings)
library(reshape2)
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
ORFLengths2 = aggregate(grl2[,"width"],by=list(grl2[,"names"]),sum)
#In orflength, find overlaps by transcripname.

rename1 = gsub("_[0-9]*","", ORFLengths1$Group.1)
ORFLengths1$Group.1 = rename1 

rename2 = gsub("_[0-9]*","", ORFLengths2$Group.1)
ORFLengths2$Group.1 = rename2 

ORFLengthsMean1 = aggregate(ORFLengths1[,"x"],by=list(ORFLengths1[,"Group.1"]),mean)
ORFLengthsMean2 = aggregate(ORFLengths2[,"x"],by=list(ORFLengths2[,"Group.1"]),mean)

totalORF = merge(ORFLengthsMean1,ORFLengthsMean2,by.x = "Group.1",by.y = "Group.1",all= F)
colnames(totalORF) = c("transcriptName","Width_nonCage","Width_Cage")
totalORFLong = melt(totalORF)
numberOfUTRSCage = length(vector2)

total2Plot = subset(totalMat,select=c("tx_name.x","teUORF.x","teCDS.x","tx_name.y","teUORF.y","teCDS.y"))
total2Plot2 = merge(total2Plot[,1:3],total2Plot[,4:6],by.x = "tx_name.x",by.y = "tx_name.y",all= F)
totalPlot3a = total2Plot2[,"teUORF.x"]/total2Plot2[,"teCDS.x"]
totalPlot3b = total2Plot2[,"teUORF.y"]/total2Plot2[,"teCDS.y"]

totalPlot
###Make means
print("summary of nonCage")
s1 = summary(totalORF$Width_nonCage)
print("summary of Cage")
s2 = summary(totalORF$Width_Cage)
means = paste("\n mean of nonCage is:",s1[4])
means = paste(means,"\n mean of Cage is",s2[4])
#index2 = totalMat[,"teUORF.y"] < 1 | is.na(totalMat[,"teUORF.y"]) | is.infinite(totalMat[,"teUORF.y"]) | is.nan(totalMat[,"teUORF.y"])
#totalMat = totalMat[!(index2[,1] | index2[,2]),]

##ggplot(totalMat,aes(x = log(/V1),y = log(V4/V3) )) + geom_point() + xlab(lab1) + ylab(lab2)
#how many
#lengths
#lengths vs te
#teUORF/teCDS vs teUORF/teCDS 
#how many datasets 
#teuorf sort
#teUORF
pdf("comparison.pdf")
theme_set(theme_gray(base_size = 18))
lengthPlot = ggplot() +
   geom_density(data = totalORFLong,aes(x = log10(value),group = variable,color = variable)) +
   xlab(paste("length of UORFs(log)",means)) + ylab("counts(log) of UORFs at width")
print(lengthPlot)
tePlot = ggplot() + 
  geom_point(data = totalMat,aes(x = log(teUORF.x),y = log(teCDS.x)),color="red") + 
  geom_point(data = totalMat,aes(x = log(teUORF.y),y = log(teCDS.y)),color="blue") + 
  xlab("teUORF value") + ylab("teCDS value") 
print(tePlot)
###Show summary
###Per transcript
dev.off()

###Compare length histogram
###te uorfs 

