# Libraries
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("ORFik")
}
if (!requireNamespace("slackr", quietly=TRUE)) {
  install.packages("slackr")
}


library(data.table)
library(ggplot2)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(Biostrings)
library(Rsamtools)
library(reshape2)
library(slackr)
library(ORFik)
library(slackr)
