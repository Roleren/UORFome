####LIBS####
library(parallel)
library(data.table)
source("./HelperVariables.R")


####VARIABLES####
cat("wd is: ",getwd())
rnaPath = "/export/valenfs/data/processed_data/RNA-seq/"
rfpPath = "/export/valenfs/data/processed_data/Ribo-seq/"
linkerPath = "./../SRA_Accessions.tab"
linkerSmallPath = p(helperMainFolder,"/linkerFileSmall.rdata")
experimentsIDPath = p(dataMainFolder,"/fantom6_Hakon_new.csv")

studyName = "Gonzalez C,2014"
speciesName = "Human"
rnaRFPMiddlePathName = "/final_results/aligned_GRCh38/" #Change this if mouse is needed!!!
bashScriptLocation = "./../prepareDataWithoutPredefinedLeaders.sh"

maxCores = as.integer(detectCores()-(detectCores()/2)) #dont use too many, 60 on kjempetuja


####FUNCTIONS####
getLinkerFile = function(){ #get files linking experiments to srr numbers
  if(!file.exists(linkerSmallPath)){
    
    
    linkerFile = fread(linkerPath,sep = "\t")
    
    linkerFileSmall = linkerFile[grep("SRR",linkerFile$Accession)]
    linkerFileSmall = linkerFileSmall[grep("GSM",linkerFileSmall$Alias)]
    save(linkerFileSmall, file = linkerSmallPath)
    
  }else{
    print("loading premade filtered linker file")
    load(linkerSmallPath)
  }
  assign("linkerFileSmall",linkerFileSmall,envir = .GlobalEnv)
}

getSpecificStudyAndSpecies = function(){
  #Choose study and species of study
  dat = read.csv(file = experimentsIDPath,sep = "\t")
  SpeciesGroup = dat[dat$Species == speciesName,]
#   for(studyName in unique(SpeciesGroup$Study)){studyGroup = SpeciesGroup[SpeciesGroup$Study == studyName,]
#   print(studyGroup)
#   }
  
  SpeciesStudyGroup = SpeciesGroup[SpeciesGroup$Study == studyName,]
  
  
  
  #linkerFileSmall[grep(SpeciesStudyGroup$Sample_ID[1],linkerFileSmall$Alias)]$Accession
  
  sorted = SpeciesStudyGroup[order(SpeciesStudyGroup$Sample_description),]
  
  if((nrow(sorted) %% 2 != 0) || (length(unique(sorted$Sample_description)) != nrow(sorted)/2) || ( sum(sorted$Sample_Type == "RPF") != (sum(sorted$Sample_Type == "RNA") ))){
    stop("Not equal number of matching rna and rfp files!")
  }
  #Make srr
  getSRRs = function(x){
    linkerFileSmall[grep(x,linkerFileSmall$Alias)]$Accession
  }
  SRR = lapply(sorted$Sample_ID, function(x) getSRRs(x))
  
  
  
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
  sorted
}


getNameVariantsForCageTissue = function(sorted){
  cellLine = sorted$Tissue_or_CellLine[1]
  if(cellLine == "HeLa"){
    cellLine = "Ovary"
  }else if(cellLine == "HEK293" ){
    cellLine = "kidney"
  }else if(cellLine == "THP-1"){
    cellLine = "macrophage"
  }
  
  currentCageFiles = grep(cellLine,list.files(cageFiles),ignore.case = T,value = T)
  
  currentCageFiles = grep(".bed",currentCageFiles,ignore.case = T,value = T)
  currentCageFiles = paste0(cageFiles,currentCageFiles)
}
getCageFiles = function(sorted,speciesName){
  ###Make cage files specific for rna seq and rfp
  if(speciesName == "Human"){
    cageFiles = cageFolder
  }else if(speciesName == "Mouse"){
    cageFiles = "/export/valenfs/projects/uORFome/DATA/CAGE/mouse/"
  }
  assign("cageFiles",cageFiles,envir = .GlobalEnv)
  currentCageFiles = getNameVariantsForCageTissue(sorted)
  if(length(currentCageFiles) != length(unique(currentCageFiles))){
    print("Warning duplicate cage files!")
  }
  currentCageFiles
}

runExperiments = function(sorted,currentCageFiles){
  #Do all combinations of cage, rna-seq and rfp and call bash script to run parallel
  a = 0 #number of cores
  tooManyRunning = F # bool for number of cores exeeds max allowed
  
  usedExperiments = sorted
  
  for(i in 1:(nrow(usedExperiments)/2)){
    if(usedExperiments[i*2-1,"Sample_Type"] == "RPF"){ #choose which is rfp and rna
      rfpIndex = i*2-1
      rnaIndex = i*2
    }else{
      rnaIndex = i*2-1
      rfpIndex = i*2
    }
    rnaName = as.character(usedExperiments[rnaIndex,]$RnaRfpFolders)
    rfpName = as.character(usedExperiments[rfpIndex,]$RnaRfpFolders)
    for(n in currentCageFiles){
      a = a+1
      if(a > maxCores)
        tooManyRunning = T
      system(paste(bashScriptLocation,T,n,rnaName,rfpName),wait = tooManyRunning) 
    }
  }
  cat("Running ",a," cores, dont do too many!!")
}

