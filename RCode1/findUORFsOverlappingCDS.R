require(data.table)
###Find all uorf into the cds, check if reading frame changes. %3 = 0
#rm(list = ls())
#load("/export/valenfs/projects/uORFome/test_results/rangesOfUORFs/CD133%2b%20stem%20cells%20-%20adult%20bone%20marrow%20derived%2c%20pool1.CNhs12552.12224-129F1.RangesUorf.rdata",envir = .GlobalEnv)
#Gtf = makeTxDbFromGFF("/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf")
#cds = cdsBy(Gtf,"tx",use.names = T)

overlap1 = findOverlaps(cds,rangesOfuORFs)
overlapCount = countOverlaps(cds,rangesOfuORFs)
numberOfOverlaps = sum(overlapCount == 1)
overlapHitsIndex = overlapCount[overlapCount == 1]


#get start cds, end uorf, find difference, check mod3
#check start upstream, ending downstream
df.rangesOfuORFs = as.data.frame(rangesOfuORFs)
df.rangesOfuORFsOnlyFirstExon = setDT(df.rangesOfuORFs)[, if(.N>1) if(strand=="+") head(.SD, 1) else tail(.SD,1) else .SD , by = names]
df.rangesOfuORFsOnlyFirstExon = as.data.frame(df.rangesOfuORFsOnlyFirstExon)
df.pos = df.rangesOfuORFsOnlyFirstExon[df.rangesOfuORFsOnlyFirstExon$strand == "+",]
df.neg = df.rangesOfuORFsOnlyFirstExon[df.rangesOfuORFsOnlyFirstExon$strand == "-",]
df.pos.finished = df.pos[,cbind("names","start","strand","seqnames")]
df.neg.finished = df.neg[,cbind("names","end","strand","seqnames")]

transcriptNames.pos = gsub("_[0-9]*","", df.pos.finished$names)
transcriptNames.neg = gsub("_[0-9]*","", df.neg.finished$names)

df.pos.finished$names = transcriptNames.pos
df.neg.finished$names = transcriptNames.neg

df.cds = as.data.frame(cds)
df.cds.onlyFirstExon = setDT(df.cds)[, if(.N>1) if(strand=="+") head(.SD, 1) else tail(.SD,1) else .SD , by = group_name]
df.cds.onlyFirstExon = as.data.frame(df.cds.onlyFirstExon)

df.cds.pos = df.cds.onlyFirstExon[df.cds.onlyFirstExon$strand == "+",]
df.cds.neg = df.cds.onlyFirstExon[df.cds.onlyFirstExon$strand == "-",]



posMerge = merge(df.cds.pos,df.pos.finished[,c("names","start","strand","seqnames")],by.x="group_name",by.y="names")
posMergesresult = posMerge[which(posMerge$start.x == posMerge$start.y),]

a = df.pos.finished[which(df.pos.finished$start != df.cds.pos$start),]

b = merge(df.cds.neg,df.neg.finished[,c("names","start","strand","seqnames")],by.x="group_name",by.y="names")

start(rangesOfuORFs)






