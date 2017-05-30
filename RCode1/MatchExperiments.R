library(data.table)
rnaPath = "/export/valenfs/data/processed_data/RNA-seq/"
rfpPath = "/export/valenfs/data/processed_data/Ribo-seq/"
linkerPath = "/export/valenfs/projects/uORFome/SRA_Accessions.tab"
linkerSmallPath = "/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data/linkerFileSmall.rdata"
experimentsIDPath = "/export/valenfs/projects/uORFome/DATA/fantom6_Hakon_new.csv"
studyName = "Gonzalez C,2014"
speciesName = "Human"
rnaRFPMiddlePathName = "/final_results/aligned_GRCh38/" #Change this if mouse is needed!!!
bashScriptLocation = "/export/valenfs/projects/uORFome/prepareDataWithoutPredefinedLeaders.sh"

if(!file.exists(linkerSmallPath)){
  dat = read.csv(file = experimentsIDPath,sep = "\t")
  
  linkerFile = fread(linkerPath,sep = "\t")
  
  linkerFileSmall = linkerFile[grep("SRR",linkerFile$Accession)]
  linkerFileSmall = linkerFileSmall[grep("GSM",linkerFileSmall$Alias)]
  save(linkerFileSmall, file = linkerSmallPath)
}else
  load(linkerSmallPath)
#Choose study and species of study
studyGroup = dat[dat$Study == studyName,]

SpeciesStudyGroup = studyGroup[studyGroup$Species == speciesName,]

linkerFileSmall[grep(SpeciesStudyGroup$Sample_ID[1],linkerFileSmall$Alias)]$Accession

sorted = SpeciesStudyGroup[order(SpeciesStudyGroup$Sample_description),]

if(nrow(sorted) %% 2 != 0){ # bad test!!!!!!!!
  cat("warning, not equal number of rna and rfp files!")
  stop()
}
#Make srr
SRR = lapply(sorted$Sample_ID, function(x) getSRRs(x))

getSRRs = function(x){
  linkerFileSmall[grep(x,linkerFileSmall$Alias)]$Accession
}

study = strsplit(as.character(sorted$Study[1])," ")[[1]][1]
sorted$SRR = SRR

####Make rna seq and rfp paths
currentRnaPath = paste0(rnaPath,grep(study,list.files(rnaPath),ignore.case = T,value = T),rnaRFPMiddlePathName)
currentRfpPath = paste0(rfpPath,grep(study,list.files(rfpPath),ignore.case = T,value = T),rnaRFPMiddlePathName)

folders = data.frame(matrix("A",  nrow = 1,ncol = 1))
for(i in 1:nrow(sorted)){ #Scary way to do it!!!!!!
  if(sorted[i,]$Sample_Type == "RPF"){
    folders[i] = paste0(currentRfpPath,grep(sorted[i,]$SRR,list.files(currentRfpPath),ignore.case = T,value = T))
  }else{
    folders[i] = paste0(currentRnaPath,grep(sorted[i,]$SRR,list.files(currentRnaPath),ignore.case = T,value = T))
  }
}
folders = t(folders)
sorted$RnaRfpFolders = folders

###Make cage files specific for rna seq and rfp
cageFiles = "/export/valenfs/projects/uORFome/DATA/CAGE/human/"
currentCageFiles = grep(sorted$Tissue_or_CellLine[1],list.files(cageFiles),ignore.case = T,value = T)
currentCageFiles = grep(".bed",currentCageFiles,ignore.case = T,value = T)
currentCageFiles = paste0(cageFiles,currentCageFiles)
if(length(currentCageFiles) != length(unique(currentCageFiles))){
  print("Warning duplicate cage files!")
}
rnaSeq = as.character(sorted[2,]$RnaRfpFolders)
rfpSeq = as.character(sorted[1,]$RnaRfpFolders)

###run test experiment
#getMatrix(usingNewCage = T,cageName = currentCageFiles[1],rnaSeq = rnaSeq,rfpSeq = rfpSeq)

#Do all combinations of cage, rna-seq and rfp and call bash script to run parallel
setwd("/export/valenfs/projects/uORFome/")
a = 0
tooManyRunning = F
for(i in 1:nrow(sorted)/2){ #NB watch the indeces for rna and rfp, might be wrong!!!!
  rnaIndex = i*2-1
  rfpIndex = i*2
  rnaName = as.character(sorted[rnaIndex,]$RnaRfpFolders)
  rfpName = as.character(sorted[rfpIndex,]$RnaRfpFolders)
  for(n in currentCageFiles){
    a = a+1
    if(a > 60)
      tooManyRunning = T
    system(paste(bashScriptLocation,T,n,rnaName,rfpName),wait = tooManyRunning) 
    }
}
cat("Running ",a," cores, dont do too many!!")