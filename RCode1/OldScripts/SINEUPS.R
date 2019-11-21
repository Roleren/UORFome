# SINES
rm(list=ls())
setwd("/export/valenfs/projects/uORFome/RCode1/") 
source("./pipelineSetup.R")
source("./TempScripts/tcp_pipeline.R")
getFasta()
sineFolder <- "/export/valenfs/projects/SINEUPS/SINEUPS_tcp_seq"
mainFolder <- sineFolder

# !!!!!!!!!!!!!!OLD WAY !!!!!!!!!!!!!!!!!!!!!!!!!!
# pfam <- data.table::fread(file  = "/export/valenfs/projects/uORFome/Annotations/Repeats/hg38_dfam.hits", sep = "\t", header = TRUE)
# 

# 
# seqlevels(fa)
# 
# a <- pfam[pfam$`#seq_name` %in% seqlevels(fa), ]
# a[,`sq-len` := NULL]; a[,kimura_div := NULL]; a[,family_acc := NULL]
# 
# a <- a[a$`e-value` < 1e-50,] # sum(a$`e-value` < 1e-50)/length(a$`e-value`)
# 
# strand <- (a$strand == "+")
# starts <- rep.int(0, length(strand))
# ends <- rep.int(0, length(strand))
# 
# starts[strand] <- a$`ali-st`[strand]
# starts[!strand] <- a$`ali-en`[!strand]
# ends[strand] <- a$`ali-en`[strand]
# ends[!strand] <- a$`ali-st`[!strand]
# 
# gr <- GRanges(seqnames = a$`#seq_name`, IRanges(starts, ends),  strand = a$strand)
# gr <- gr[(width(gr) > 100 & width(gr) < 400)] # SINES definition
# samples <- sample(x = seq.int(length(gr)), size = 1000000)
# gr <- gr[samples]

# !!!!!!!!!!!!!!!!!!! NEW WAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# load sines from UCSC database
gr <- fread.bed(filePath = paste0(sineFolder, "/Repeats/SINEs.bed.gz"))
names(gr) <- paste(seqnames(gr) ,start(gr), end(gr), strand(gr))

countsPerLibraryOverTranscriptPerSubunit(gr)

b <- fread(file =  p(sineFolder,"/Repeats/countsPerRepeat.csv"))
b[, txNames := NULL]; b[, V1 := NULL]

c <- melt(b)
colnames(c)[2] <- "counts"
library(ggplot2)
p1 <- ggplot(c, aes(x=variable, y=counts)) + 
  geom_violin() + theme(axis.text.x = element_text(face="bold", color="#993333", 
                                                   size=14, angle=45)) + 
  ggtitle("TCP-seq counts", subtitle = "distribution per library") + 
  ylab(label = "distribution of counts per transcript")

ggsave(plot = p1, filename = p(mainFolder, "/Repeats/distributions/distribution_violin.pdf"))

d <- c[, .(counts = sum(counts)), by = variable]

p2 <- ggplot(d, aes(x=variable, y = counts)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(face="bold", color="#993333",  size=14, angle=45)) + 
  ggtitle("TCP-seq counts", subtitle = "sum per library")

ggsave(plot = p2, filename = p(mainFolder, "/Repeats/sum_plots/sum_per_library.pdf"))

e <- readLengthsPerLibrary()
e <- melt(e)

f <- copy(d)
f$counts <- f$counts / (e$value / 1000000)

p3 <- ggplot(f, aes(x=variable, y = counts)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(face="bold", color="#993333",  size=14, angle=45)) + 
  ggtitle("TCP-seq counts", subtitle = "sum per library normalized") + 
  ylab(label = "TPM")
ggsave(plot = p3, filename = p(mainFolder, "/Repeats/sum_plots/sum_per_library_normalized.pdf"))

# top SINEs
g <- rowSums(b)
best <- order(g, decreasing = TRUE)
top100 <- gr[best[1:100]]
export.bed(top100, con = p(mainFolder, "/Repeats/SINEs_top_100_counts.bed"))

gtfPath <- p(dataFolder, "/Zebrafish/zebrafish_GRCh10_81.gtf.db")
txdb <- loadDb(gtfPath); seqlevelsStyle(txdb) <- seqlevelsStyle(gr)
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
outSide <- countOverlaps(gr, tx) == 0


topSinesOutside <- gr

# dot plots

h <- b[g > 1,]
gr <- gr[g > 1]

dt <- data.table(id = names(gr))
for(i in 1:ncol(h)) {
  dt <- cbind(dt, ORFik:::fpkm_calc(h[,get(colnames(h)[i])], lengthSize = width(gr), librarySize = e$value[i]))
}
dt[, id := NULL]
colnames(dt) <- colnames(h)

df = getTCPdf()
for (i in 1:nrow(df)) {
  print(i)
  stage <- df$stage[i]
  type <- df$type[i]
  max <- max(quantile(dt[,get(colnames(h)[(i*2)-1])], 0.95), quantile(dt[,get(colnames(h)[(i*2)])], 0.95))
  p4 <- ggplot(dt, aes(x=dt[,get(colnames(h)[(i*2)-1])], y=dt[,get(colnames(h)[(i*2)])])) + 
    geom_jitter(width = .3, size=0.5) + ylim(0, max) + xlim(0, max) + 
    ggtitle(label = "SINEs Distribution TCP-seq") + 
    xlab("SSU") + ylab("LSU")
  
  ggsave(plot = p4, filename = paste0(mainFolder, "/Repeats/distributions/", stage,"_", type, "_dotplot_distribution.pdf"))
}

f <- melt(dt)[, .(counts = sum(value)), by = variable]
p3 <- ggplot(f, aes(x=variable, y = counts)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(face="bold", color="#993333",  size=14, angle=45)) + 
  ggtitle("TCP-seq counts", subtitle = "sum per library normalized") + 
  ylab(label = "FPKM counts")
ggsave(plot = p3, filename = p(mainFolder, "/Repeats/sum_plots/sum_per_library_normalized_fpkm.pdf"))

# Write out wig files
df = getTCPdf()

for (i in 1:nrow(df)){
  stage <- df$stage[i]
  type <- df$type[i]
  GRToWig(gr = df$SSU, outputPrefix = paste0(mainFolder, "/tcp_wig/SSU_",stage,"_",type,"_"))
  GRToWig(gr = df$LSU, outputPrefix = paste0(mainFolder, "/tcp_wig/LSU_",stage,"_",type,"_"))
}


df <- c("/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_2hpf_Tota_mRNA.bam",
        "/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_6hpf_Total_mRNA_merged.bam")
for (i in 1){
  stage <- c("64", "shield")
  type <- "1"
  GRToWig(gr = df[i], outputPrefix = paste0(mainFolder, "/RNAseq/RNA_",stage[i],"_",type,"_"))
  GRToWig(gr = df[i+1], paste0(mainFolder, "/RNAseq/RNA_",stage[i+1],"_",type,"_"))
}


# count
og <- g[best]
valid <- og[og > 10]
validOnes <- data.frame(values = valid)
pa <- ggplot(data = validOnes, aes(values)) + geom_density() + xlim(10, 690) + xlab("counts of tcp seq") +
  ggtitle(label = "SINE count distribution", subtitle = p("All SINEs with > 10 reads, \nnumber of SINEs: ", nrow(validOnes)))
ggsave("../test.pdf", plot = pa)
