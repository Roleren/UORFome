#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(gridExtra)

##
# 1 unshifted heatmaps
##

start_m5 <- read.table(args[1])

mono.start5.raw <-

start_m5_50 <- subset(start_m5, V1<=50 & V1 >= -100);
start_m5_50 <- start_m5_50 %>% group_by(V2) %>% mutate(LenSum=sum(V3), LenSD=sd(V3, na.rm = T), LenMean=mean(V3)) %>% mutate(LenZ=(V3-LenMean)/LenSD)
start_m5_50 <- start_m5_50 %>% group_by(V1) %>% mutate(PosSum=sum(V3), PosSD=sd(V3, na.rm = T), PosMean=mean(V3)) %>% mutate(PosZ=(V3-PosMean)/PosSD)


#hack to fix plotting problem with NaN values and the heatmap legend
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

start_m5_50[is.nan(start_m5_50)] <- NA

#mLen5<-ggplot(subset(start_m5_50, V2 <= 100 & V2 >= 15 & V1 <= 100 & V1 >= -50) , aes(x=V1, y=V2, fill=LenZ)) + geom_tile()  +
mLen5<-ggplot(subset(start_m5_50, V1 <= 50 & V1 >= -100) , aes(x=V1, y=V2, fill=LenZ)) + geom_tile()  +
  scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(-2,0,2,4,6,8,10,12), rescaler = function(x,...) x, oob = identity, name="Z-score") +
  xlab("Position relative to start codon") + ylab("Protected fragment length") +
  scale_x_continuous(limits = c(-101,51), expand = c(0, 0), breaks=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  #               scale_y_continuous(limits = c(15,100), expand = c(0, 0), breaks=c(20, 30, 40, 50, 60, 70, 80 ,90, 100), labels=c("0.20","0.30","0.40","0.50","0.60","0.70","0.80","0.90","1.00")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 12))


mPos5<-ggplot(subset(start_m5_50, V1 <= 50 & V1 >= -100) , aes(x=V1, y=V2, fill=PosZ)) + geom_tile()  +
  scale_fill_gradientn(colours=c("yellow","lightblue","blue", "navy"), values=c(-2,0,2,4,6,8,10,12), rescaler = function(x,...) x, oob = identity, name="Z-score") +
  xlab("Position relative to start codon") + ylab("Protected fragment length") +
  scale_x_continuous(limits = c(-101,51), expand = c(0, 0), breaks=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 12))


##
# 2 unshifted barplots
##


#scale the data by the sum of the window
subset_mono5.raw<-filter(mono.start5.raw, V1 >=-100, V1 <= 50)
subset_mono5.raw <- mutate(subset_mono5.raw, proportion = V2/sum(V2))

max_y <- max(subset_mono3.raw$proportion,subset_mono5.raw$proportion)

mb5 <- ggplot(data=subset_mono5.raw, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
  geom_bar(stat="identity",width=1, colour="black", size=0.3) +
  theme_bw() +
  xlab("Position relative to start codon") +
  ylab("Proportion of reads") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-101,51), expand = c(0, 0), breaks=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50)) +
  scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
  scale_fill_manual(values=c("#db811a", "#9a9a99", "#6d6e63"), name = "Frame") +
  theme(text = element_text(size = 12))

mb3 <- ggplot(data=subset_mono3.raw, aes(x=V1, y=proportion, fill=as.factor((V1 %% 3)+1))) +
  geom_bar(stat="identity",width=1, colour="black", size=0.3) +
  theme_bw() +
  xlab("Position relative to start codon") +
  ylab("Proportion of reads") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-50,100), expand = c(0, 0), breaks=c(-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100)) +
  scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
  scale_fill_manual(values=c("#db811a", "#9a9a99", "#6d6e63"), name = "Frame") +
  theme(text = element_text(size = 12))

#grid.arrange(mb5,mb3)

###
# Pretify plots
###

# Function to save legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Save the legend
#+++++++++++++++++++++++
legend_heatmapP5<- get_legend(mPos5)    
legend_heatmapL5<- get_legend(mLen5)
legend_barchart5 <- get_legend(mb5)

mPos5 <- mPos5 + theme(legend.position="none")
mLen5 <- mLen5 + theme(legend.position="none")
mb5 <- mb5 + theme(legend.position="none")

##
# Output
##

lay <- rbind(c(1,1,1,1,1,1,1,1,2,3,3,3,3,3,3,3,3,4),
             c(5,5,5,5,5,5,5,5,6,7,7,7,7,7,7,7,7,8),
             c(9,9,9,9,9,9,9,9,10,11,11,11,11,11,11,11,11,12))

grid.arrange(mb5, legend_barchart5, mLen5, legend_heatmapL5, mPos5, legend_heatmapP5, layout_matrix = lay)

out1<-arrangeGrob(mb5, legend_barchart5, mLen5, legend_heatmapL5, mPos5, legend_heatmapP5, layout_matrix = lay)
ggsave(file = outName, out1, width=500, height=400, unit="mm", dpi=300) 












