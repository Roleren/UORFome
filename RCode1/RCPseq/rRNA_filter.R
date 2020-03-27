#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# STAR aligner, run statistics
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Used to get some transcript read overlap statistics
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO on columns:
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# name = name of file,
# ,raw_reads = number of raw reads (fastq library size)
# ,raw_reads_size = percentage size compared to max set of reads (one is always 100%)
# ,not_trimmed = What percentage of reads was not trimmed away (fastq score and minimum length)
# ,rRNA_content = What percentage of reads was trimmed away by silva rRNA database
# ,aligned_to_genome = What percentage of unique reads aligned to genome
# ,aligned_to_genome_reads = How many unique reads aligned to genome
# ,ratio_aligned_raw = What percentage of unique reads aligned to genome,
# ratio_mrna_aligned = How many of the aligned reads overlap mrna (multimapped + unique reads, so number can be > 1 if overlapping genes)


library(jsonlite)
library(slackr)
folder <- "/export/valenfs/data/processed_data/TCP-seq/valen_2020_human_C11/"
subset_of_libs <- 1:13 # Which libraries
# dirname(df$filepath[1])

# Raw reads
files <- list.files(paste0(folder, "/trim"), full.names = T)
files <- grep("\\.json", files, value = T)[subset_of_libs]
bas <- gsub("_S.*", "", basename(files))
bas <- gsub("report_", x = bas, "")
dt <- data.table(name = bas)

values <- c()
values_after <- c()
for (f in files) {
  js <- jsonlite::read_json(f)
  values <- c(values, js$read1_before_filtering$total_reads)
  values_after <- c(values_after, js$read1_after_filtering$total_reads)
}
dt$raw_reads <- paste0(values)
dt$raw_reads_size <- paste0(round(values / max(values) * 100, digits = 2), "%")
dt$not_trimmed <- paste0(round((values_after / values) * 100 , digits = 2), "%")

# rRNA
files <- list.files(paste0(folder, "/rRNA_depletion/"), full.names = T)
files <- grep("Log.final.out", files, value = T)[subset_of_libs]

values <- c()
for (f in files) {
  uniqes <- fread(f, sep = "\t", skip = 8)[2,]$V2
  multi <- fread(f, sep = "\t", skip = 23)[2,]$V2
  values <- c(values, sum(as.numeric(gsub("%", "", c(uniqes,multi)))))
}
dt$rRNA_content <- paste0(values, "%")

# Aligned to genome
files <- list.files(paste0(folder, "/aligned/LOGS/"), full.names = T, include.dirs = F)
files <- grep("Log.final.out", files, value = T)
values <- c()
values_num <- c()
for (f in files) {
  uniques_num <- fread(f, sep = "\t", skip = 8)[1,]$V2
  uniqes <- fread(f, sep = "\t", skip = 8)[2,]$V2
  multi_num <- fread(f, sep = "\t", skip = 23)[1,]$V2
  multi <- fread(f, sep = "\t", skip = 23)[2,]$V2
  values <- c(values, sum(as.numeric(gsub("%", "", c(uniqes,multi)))))
  values_num <- c(values_num, sum(as.numeric(uniques_num), as.numeric(multi_num)))
}
dt$aligned_to_genome <- paste0(values, "%")
dt$aligned_to_genome_reads <- values_num
dt$ratio_aligned_raw <- paste0(round((dt$aligned_to_genome_reads / as.numeric(dt$raw_reads)) * 100, digits = 2), "%")

# To mrna lvl, from ORFikQC:
mrnastats <- paste0(folder, "/aligned/QC_STATS/STATS.csv")
dtmrna <- fread(mrnastats)[subset_of_libs,]
dt$ratio_mrna_aligned <-  round(dtmrna$ratio_mrna_aligned, 4)
dt$mrna_reads<-  dtmrna$mRNA
# Order and clean up
dt <- dt[order(raw_reads, decreasing = TRUE), ]
fwrite(dt, paste0(folder, "alignment_stats.csv"))
slackr_upload(paste0(folder, "alignment_stats.csv"), channels = "#visualizations")

# Statistics
dt$raw_reads <- as.numeric(dt$raw_reads)
dt$libtype <- c("WT", "WT", "ribo", rep("WT", 3), rep("ribo", 2), "RNA", "ribo", "WT", rep("ribo", 2))
dt[, mean(raw_reads), by = libtype]
dt[, mean(ratio_mrna_aligned), by = libtype]
dt[, median(ratio_mrna_aligned), by = libtype]
dt[, mean(mrna_reads), by = libtype]
dt[, median(mrna_reads), by = libtype]
t(dt[, .(summary = summary(mrna_reads), type = names(summary(mrna_reads))), by = libtype][1:12,])
