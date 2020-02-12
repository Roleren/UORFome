# Go analysis
# Start or stop ? contexts
# This first part depends on a lot of packages, if you only want data before ploting, go step 3
# 1: get adams data (only run if you have all packages)
source("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline.R")
library(ggpubr)
setwd("/export/valenfs/projects/Hakon/AdamVienna/")

# df <- combined.dat.ranked.filter # <- this is output from adams file, so can no run this from here.
# Location: 
#dt <- setDT(df)
# txNames <- dt$transcript_id
#saveRDS(df, file = "kozakTx.rds") 


# 2: Add go terms column
dt <- setDT(readRDS("kozakTx.rds")) # <- instead you can load this file
orfGo <- getORFsGoTerms(dt$X.gene_id, organism = "Danio rerio")

dt$go <- orfGo
res <- list()
for (i in unique(dt$initiation_sequence_sub)[1:2]) {
  res <- list(res, table(dt$go[dt$initiation_sequence_sub == i]))
}

subs <- dt$initiation_sequence_sub
res <- lapply(unique(subs), function(i) {
  table(dt$go[subs == i])
})
names(res) <- unique(dt$initiation_sequence_sub)

dat <- data.table(g0 = unique(dt$go), matrix(0, nrow = length(unique(dt$go)), ncol = length(unique(subs))))
colnames(dat)[2:ncol(dat)] <- unique(dt$initiation_sequence_sub)
for (i in 2:ncol(dat)) {
  dat[(unique(dt$go) %in% names(res[[i-1]])), i] <- as.numeric(res[[i-1]])
}

# 3: melt and filter data
dat <- dat[rowSums2(as.matrix(dat[, 2:ncol(dat)])) > 130,] # filter on significanse of GO terms, remove most of them
saveRDS(dat, file = "goData.rds")

dat <- setDT(readRDS("goData.rds")) # <----------------- Start point, load go matrix
mt <- melt(data = dat)

mt <- mt[, .(g0, value = value/sum(value)), by = variable]


# 4: plot data
plot <- ggplot(mt, aes(x = g0, y = variable, fill = value)) +
  geom_tile() + scale_fill_gradientn(colours = c("white", "yellow2", "yellow3", "lightblue", "blue", "navy")) + 
  xlab("TIS context") + 
  ylab("Protected fragment length") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + theme(text = element_text(size = 12)) + 
  theme(axis.text.x = element_text(angle=65, vjust=0.5, size = 12))

plot
ggsave("goheatmap.pdf", plot, width = 25, height = 30, dpi = 300)

# new version
library(clusterProfiler)
library(org.Dr.eg.db)

dt <- setDT(readRDS("kozakTx.rds")) # <- instead you can load this file
subs <- dt$initiation_sequence_sub
# valid <- names(table(subs))[table(subs) > 20]
valid <- valid[valid %in% c("TGGAATG", "AAACATG", "AAGCATG")]

# Do gorilla
for (v in valid) {
  df <- data.frame(genes = unique(dt$X.gene_id[dt$initiation_sequence_sub %in% v]))
  write.table(df, file = paste0(v, ".csv"))
  ddf <- data.frame(genes = unique(dt$X.gene_id[!(dt$X.gene_id %in% df$genes)]))
  write.table(ddf, file = paste0(v, "_back.csv"))
}
ggg <- unique(dt$X.gene_id[dt$go =="ion binding"])

# AAGCATG hits http://cbl-gorilla.cs.technion.ac.il/GOrilla/l6at4vmc/GOResults.html
# AAACATG hits http://cbl-gorilla.cs.technion.ac.il/GOrilla/ty5183wa/GOResults.html
# TGGAATG no hits


ego <- lapply(valid, function(i) {
  a <- clusterProfiler::enrichGO(gene = dt$X.gene_id[subs == i],
                                 OrgDb = org.Dr.eg.db,
                                 keyType = 'ENSEMBL',
                                 ont = "ALL")
  # if(is.null(a)) return(NULL)
  # p <- dotplot(a, 10)
  # p + ggtitle(label = i)
  return(a)
})

a <- sapply(unique(subs), function(i) {sum(subs == i)})
p <- dotplot(ego, showCategory=30)
library(gridExtra)
library(grid)
pp <- do.call("grid.arrange", rectGrob(ego))

pp <- lapply(ego, function(ii) grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii))) 
b <- grid.arrange(grobs=pp, ncol=2, 
                  top="top label", bottom="bottom\nlabel", 
                  left="left label", right="right label")
bb <- marrangeGrob(ego, nrow = ceiling(102 / 2) , ncol=2)

b <- unlist(ego)


ggsave("multipage.pdf", ml, width = 23, height = 49)


pdf("abc.pdf")
for(i in ego) {
  result = tryCatch({
    print(i)
  }, warning = function(w) {
    0
  }, error = function(e) {
    0
  }, finally = {
    0
  })
} 
dev.off()

# Analysis

dt <- setDT(readRDS("kozakTx.rds")) # <- instead you can load this file
subs <- dt$initiation_sequence_sub
valid <- names(table(subs))
# valid <- names(table(subs))[table(subs) > 20]
v <- valid[valid %in% c("TGGAATG", "AAACATG", "AAGCATG")]

d <- data.table()
for (i in v) {
  a <- clusterProfiler::enrichGO(gene = dt$X.gene_id[subs == i],
                                 OrgDb = org.Dr.eg.db,
                                 keyType = 'ENSEMBL',
                                 ont = "ALL", 
                                 readable = TRUE)
  if(nrow(a@result) == 0) next
  
  d <- rbind(d, cbind(a@result, context = i))
}

write.csv2(d, file = "enrichedTerms.csv")

e <-  data.table()
for (i in v) {
  a <- as.character(dt$X.gene_id[subs == i])
  if(length(a) == 0) next
  
  e <- rbind(e, cbind(gene = a, context = i, symbol = geneSymbolsTo(a, org.db = org.Dr.eg.db)))
}
write.csv2(e, file = "genes.csv")


v <- valid[valid %in% c("AAAAATG")]
e <-  data.table()
for (i in v) {
  a <- as.character(dt$X.gene_id[subs == i])
  if(length(a) == 0) next
  
  e <- rbind(e, cbind(gene = a, context = i, symbol = geneSymbolsTo(a, org.db = org.Dr.eg.db)))
}
write.csv2(e, file = "genesInThe7Group.csv")