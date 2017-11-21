
#PLAN: Plot the healthy vs unhealthy, first start with just total global difference
#Then we need to aggregate by uorf ID
#use same cds matrix for all same uorfs
#the have 1 column each for each match ? bad idea I think

#!! how did I get max ucdist of 600k????

#Col1 is healthy, Col2 is unhealthy, if there is such a difference
#Else the difference is specified
library(ggplot2)

library(rtracklayer)


lab1 = "teCDS D/H"
lab2 = "teUORF D/H"

plottingFolder = p(resultsFolder,"/Plotting/Comparisons_plots/")



# ag1 = aggregate( .~tx_name  ,data = healthyMatrix,FUN = sum) #fix this to DT
# ag2 = aggregate( .~tx_name  ,data = sickMatrix,FUN = sum)


# totalMat = merge(ag1, ag2, by.x = "Group.1", by.y = "Group.1", all = F)


##ggplot(totalMat,aes(x = log(/V1),y = log(V4/V3) )) + geom_point() + xlab(lab1) + ylab(lab2)
pdf(p(plottingFolder,"comparison.pdf"))

###Global differences plots###
lengthPlot = ggplot() + 
  geom_density(data = healthMatrix,aes(x = log10(width)) ,color="orange") +
  geom_density(data = sickMatrix,aes(x = log10(width)) ,color="red") +
  xlab("length of uorfs (log)") + ylab("number of counts") +
  labs(title="Length of UORFs")
print(lengthPlot)
tePlot = ggplot() + 
  geom_density(data = healthMatrix,aes(x = log10(teUORF),color="red")) + 
  geom_density(data = sickMatrix,aes(x = log10(teUORF),color="blue"))  + 
  xlab("teUORF (log)") + ylab("ratio of counts") 
print(tePlot)
teCDSPlot = ggplot() + 
  geom_density(data = healthMatrix,aes(x = log10(teCDS),color="red")) + 
  geom_density(data = sickMatrix,aes(x = log10(teCDS),color="blue")) + 
  xlab("teCDS (log)") + ylab("ratio of counts") 
print(teCDSPlot)
rankPlot = ggplot() +
  geom_density(data = healthMatrix,aes(x = rank) ,color="orange") +
  geom_density(data = sickMatrix,aes(x = rank) ,color="red") +
  xlab("exon rank(from left) of uorfs ") + ylab("ratio of counts") +
  xlim(1,100)
print(rankPlot)
ucDistPlot = ggplot() +
  geom_density(data = healthMatrix,aes(x = UCdist) ,color="orange") +
  geom_density(data = sickMatrix,aes(x = UCdist) ,color="red") +
  xlab("distance between uorf end and cds start") + ylab("ratio of counts") + 
  xlim(0,300)
print(ucDistPlot)
####Histogram of frame proportions
framePlot = ggplot(healthMatrix,aes(x = frame,) ) + geom_bar(position = "dodge") +
  labs(title="Healthy: Reading frame proportions uorf to cds start") + xlab("Frame 1,2,3") + ylab("counts of orfs per frame") 
print(framePlot)
framePlot = ggplot(sickMatrix,aes(x = frame) ) + geom_bar(position = "dodge") +
  labs(title="UnHealthy: Reading frame proportions uorf to cds start") + xlab("Frame 1,2,3") + ylab("counts of orfs per frame") 
print(framePlot)
#then plots for average of all unique uorfs summed by aggregation
#filter only the ones that are in both healthy and unhealthy ? 
#which are in many of the healthy but not unhealthy
#what way do they differ in features ?
agHealth = healthMatrix[]
#plot for te difference between each 
###Show summary
###Per transcript
dev.off()

###Compare length histogram
###te uorfs 

