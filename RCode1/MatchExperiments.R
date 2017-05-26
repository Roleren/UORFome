library(data.table)
rnaPath = "/export/valenfs/data/processed_data/RNA-seq/"
rfpPath = "/export/valenfs/data/processed_data/Ribo-seq/"

dat = read.csv(file = "/export/valenfs/projects/uORFome/DATA/fantom6_Hakon_new.csv",sep = "\t")

gonzales = dat[dat$Study == "Gonzalez C,2014",]

humanGonzales = gonzales[gonzales$Species == "Human",]

linkerFile = fread("/export/valenfs/projects/uORFome/SRA_Accessions.tab",sep = "\t")

linkerFileSmall = linkerFile[grep("SRR",linkerFile$Accession)]
linkerFileSmall = linkerFileSmall[grep("GSM",linkerFileSmall$Alias)]
save(linkerFileSmall, file = "linkerFileSmall.rdata")
load("linkerFileSmall.rdata")

linkerFileSmall[grep(humanGonzales$Sample_ID[1],linkerFileSmall$Alias)]$Accession

sorted = humanGonzales[order(humanGonzales$Sample_description),]

if(nrow(sorted) %% 2 != 0){
  cat("warning, not equal number of rna and rfp files!")
  stop()
}

SRR = lapply(sorted$Sample_ID, function(x) getSRRs(x))

getSRRs = function(x){
  linkerFileSmall[grep(x,linkerFileSmall$Alias)]$Accession
}

study = strsplit(as.character(sorted$Study[1])," ")[[1]][1]
sorted$SRR = SRR
currentRnaPath = paste0(rnaPath,grep(study,list.files(rnaPath),ignore.case = T,value = T),"/final_results/aligned_GRCh38/")
currentRfpPath = paste0(rfpPath,grep(study,list.files(rfpPath),ignore.case = T,value = T),"/final_results/aligned_GRCh38/")

folders = data.frame(matrix("A",  nrow = 1,ncol = 1))
for(i in 1:nrow(sorted)){
  if(sorted[i,]$Sample_Type == "RPF"){
    folders[i] = paste0(currentRfpPath,grep(sorted[i,]$SRR,list.files(currentRfpPath),ignore.case = T,value = T))
  }else{
    folders[i] = paste0(currentRnaPath,grep(sorted[i,]$SRR,list.files(currentRnaPath),ignore.case = T,value = T))
  }
}
folders = t(folders)
sorted$RnaRfpFolders = folders

cageFiles = "/export/valenfs/projects/uORFome/DATA/CAGE/human/"
currentCageFiles = grep(sorted$Tissue_or_CellLine[1],list.files(cageFiles),ignore.case = T,value = T)
currentCageFiles = grep(".bed",currentCageFiles,ignore.case = T,value = T)

if(length(currentCageFiles) != length(unique(currentCageFiles))){
  print("Warning duplicate cage files!")
}

getMatrix(usingNewCage = T,cageName = currentCageFiles[1],rnaSeq = sorted[2,]$RnaRfpFolders,rfpSeq = sorted[1,]$RnaRfpFolders)
