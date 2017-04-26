#4 Create statistics
par(mai=c(0.5,0.5,0.5,0.5))
#Correlation 5'utr and cds and plot
#correlation = cor(normUTR,normUTRBED)
#print(correlation)
#plot(normUTR,normUTRBED)

#Correlation te's
plot(log(as.numeric(matrix$teRNA)),log(as.numeric(matrix$V2)))

correlation = cor(matrix[,"teRNA"],matrix[,"teUTR"])
print(correlation)
plot(x = matrix[,"teRNA"],y = matrix[,"teUTR"])
cor.test(matrix[,"teRNA"],matrix[,"teUTR"])

### TODO LIST ###
#1. statistics:
# Test if the utrs have a high translational efficiancy (t test, wilcoxen)
hist(fiveUTRs) #get index
index = 1

t-test = t.test(cds[0:index],cds[index+1:end])

#a. See if they have homogeneous variance
if(var.test(fiveUTRs,cds)$p.value < 0.05){
   T = t.test(x = fiveUTRs,y = cds,alternative = "two.sided",var.equal = T)
   W = wilcox.test(x = fiveUTRs,y = cds,alternative = "two.sided", var.equal = T)
}

 
 
 #create pdf of plot
 pdf("pdfOfPlot.pdf")
