# test 
fastq <- readDNAStringSet("/export/valenfs/projects/uORFome/WT_SSU.fastq")

pos <- startRegion(l, extendLeaders(l, 76), T, 1, -1)
match <- countOverlaps(A16_SSU_5p, pos)
what <- which((match > 0))
test <- A16_SSU_5p[-what]
regionBarPlot(l, extendLeaders(l, 76), outdir = NULL, test, upstream = upshort[1], downstream = downshort[1])

fasts <- fastq[what]
subs <- substring(fasts, 1, 1)
table(subs)
gs <- subs == "G"
total_filter <- what[gs]
res <- A16_SSU_5p[total_filter]
start(res[strandBool(res)]) <- start(res[strandBool(res)]) + 1
end(res[strandBool(res)]) <- end(res[strandBool(res)]) + 1
start(res[!strandBool(res)]) <- start(res[!strandBool(res)]) - 1
end(res[!strandBool(res)]) <- end(res[!strandBool(res)]) - 1

results <- A16_SSU_5p
results[total_filter] <- res
regionBarPlot(l, extendLeaders(l, 76), outdir = NULL, results, upstream = upshort[1], downstream = downshort[1])
