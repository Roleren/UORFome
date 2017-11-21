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

#studyName = "Hsieh AC,2012"  #"Gonzalez C,2014"
speciesName = "Human"
rnaRFPMiddlePathName = "/final_results/aligned_GRCh38/" #Change this if mouse is needed!!!
bashScriptLocation = "./../prepareDataWithoutPredefinedLeaders.sh"

maxCores = as.integer(detectCores()-(detectCores()/20)) #dont use too many, 60 on kjempetuja

filterBadStudies = function(){
  dat = read.csv(file = experimentsIDPath,sep = "\t")
  SpeciesGroup = dat[dat$Species == speciesName,]
  
  #filter micro rna studies
  SpeciesGroup = SpeciesGroup[-as.numeric(grep("mir",SpeciesGroup$Sample_description,ignore.case = T)),]
  SpeciesGroup =  SpeciesGroup[!(SpeciesGroup$Study == "Guo H,2010"),]
  for(study in unique(SpeciesGroup$Study)){ #for each study, verify it
    thisStudy = SpeciesGroup[SpeciesGroup$Study == study,]
    if(( sum(thisStudy$Sample_Type == "RPF") != (sum(thisStudy$Sample_Type == "RNA") ))){
      cat("Not equal number of matching rna and rfp files, in !",study,"\n")
      
      drop = SpeciesGroup$Study == study
      SpeciesGroup = SpeciesGroup[!drop,]
    }
  }
  
  
  assign("SpeciesGroup",SpeciesGroup,envir = .GlobalEnv)
}
####FUNCTIONS####
getLinkerFile = function(){ #get files linking experiments to srr numbers
  if(!file.exists(linkerSmallPath)){ # this does not include srs studies
    
    linkerFile = fread(linkerPath,sep = "\t")
    
    linkerFileSmall = linkerFile[grep("SRR",linkerFile$Accession)]
    d = linkerFileSmall[grep("SRP029589",linkerFileSmall$study)] # cheat to add stumpf, needs to be automated in package
    linkerFileSmall = linkerFileSmall[grep("GSM",linkerFileSmall$Alias)]
    linkerFileSmall = rbind(linkerFileSmall, d[d$Study == "SRP029589",])
    save(linkerFileSmall, file = linkerSmallPath)
    
  }else{
    print("loading premade filtered linker file")
    load(linkerSmallPath)
  }
  assign("linkerFileSmall",linkerFileSmall,envir = .GlobalEnv)
}

makeCorrectExperimentNames = function(sorted,rnaRows,rpfRows){
  print("matching") 
  tempSorted = sorted
  a = 1
  for(i in 1:nrow(sorted)){
    if(i%%2==1)
      sorted[i,] = tempSorted[tempSorted$Sample_description == rnaRows[a],]
    else{
      sorted[i,] = tempSorted[tempSorted$Sample_description == rpfRows[a],]
      a = a+1
    }
  }
  return(sorted)
}

checkExperimentNamesAreCorrect = function(sorted){
  
  #start with mRNA seqs
  #rnaSeqs = sorted[sorted[,"Sample_Type"] =="RNA",]
  if(length(unique(sorted$Sample_description)) != nrow(sorted)/2){#fix naming order
    naming = sorted[sorted[,"Sample_Type"] =="RNA",]$Sample_description
    rnaRows = naming[grep("mRNA",naming,ignore.case = T)]
    
    if((length(rnaRows) == nrow(sorted)/2) & identical(naming,rnaRows) ){
      print("hit rna")
      rnaRemoved = gsub(pattern = "mRNA",replacement = "",x = rnaRows,ignore.case = T)
    }else{
      rnaRows = naming[grep("RNA",naming,ignore.case = T)]
      if((length(rnaRows) == nrow(sorted)/2) & identical(naming,rnaRows) ){
        print("hit rna")
        rnaRemoved = gsub(pattern = "RNA",replacement = "",x = rnaRows,ignore.case = T)
      }
      else{
        stop(cat("could not find rna seqs for study",as.character(sorted$Study[1])))
      }
    }
      
      naming = sorted[sorted[,"Sample_Type"] =="RPF",]$Sample_description
      
      rpfRows = naming[grep("Ribosome profile",naming,ignore.case = T)]
      if((length(rpfRows) == nrow(sorted)/2) & identical(naming,rpfRows) ){ #if matching
        print("hit rfp")
        rpfRemoved = gsub(pattern = "ribo",replacement = "",x = rpfRows,ignore.case = T)
      }else{
        rpfRows = naming[grep("ribo",naming,ignore.case = T)]
        if((length(rpfRows) == nrow(sorted)/2) & identical(naming,rpfRows) ){ #if matching
          print("hit rfp")
          rpfRemoved = gsub(pattern = "ribo",replacement = "",x = rpfRows,ignore.case = T)
        }else{
          rpfRows = naming[grep("footprints",naming,ignore.case = T)]
          if((length(rpfRows) == nrow(sorted)/2) & identical(naming,rpfRows) ){
            print("hit rfp")
            rpfRemoved = gsub(pattern = "footprints",replacement = "",x = rpfRows,ignore.case = T)
          }else{
            rpfRows = naming[grep("footprint",naming,ignore.case = T)]
            if((length(rpfRows) == nrow(sorted)/2) & identical(naming,rpfRows) ){
              print("hit rfp")
              rpfRemoved = gsub(pattern = "footprint",replacement = "",x = rpfRows,ignore.case = T)
            }else{
              rpfRows = naming[grep("RPF",naming,ignore.case = T)]
              if((length(rpfRows) == nrow(sorted)/2) & identical(naming,rpfRows) ){
                print("hit rfp")
                rpfRemoved = gsub(pattern = "RPF",replacement = "",x = rpfRows,ignore.case = T)
              }else{
                stop(cat("could not find ribo seqs for study",as.character(sorted$Study[1])))
              }
            }
          }
        }
      }
      
      if(identical(rnaRemoved,rpfRemoved)){ #ready to match, else remove more
        return( makeCorrectExperimentNames(sorted,rnaRows,rpfRows) )
      }else{
        #try to fix rna
        rnaRemoved = gsub(pattern = "-seq",replacement = "",x = rnaRemoved)
        if(identical(rnaRemoved,rpfRemoved)){
          return( makeCorrectExperimentNames(sorted,rnaRows,rpfRows) )
        }else{
          rnaRemoved = gsub(pattern = "seq",replacement = "",x = rnaRemoved)
          if(identical(rnaRemoved,rpfRemoved)){
            return( makeCorrectExperimentNames(sorted,rnaRows,rpfRows) )
          }
        }
        stop(cat("could not sort study",as.character(sorted$Study[1])))
      }
      
    
    #resort
  }
  return(sorted)
}

getSpecificStudyAndSpecies = function(studyName){
  #Choose study and species of study
  
  SpeciesStudyGroup = SpeciesGroup[SpeciesGroup$Study == studyName,]
  
  sorted = SpeciesStudyGroup[order(SpeciesStudyGroup$Sample_description),]
  
  if((nrow(sorted) %% 2 != 0) || ( sum(sorted$Sample_Type == "RPF") != (sum(sorted$Sample_Type == "RNA") ))){
    stop("Not equal number of matching rna and rfp files!")
  }
  
  sorted = checkExperimentNamesAreCorrect(sorted)
  
  #Make srr
  getSRRs = function(x){ #It now chooses first
    linkerFileSmall[grep(x,linkerFileSmall$Alias)]$Accession[1]
  }
  SRR = lapply(sorted$Sample_ID, function(x) getSRRs(x))
  if(is.na(SRR[1])){
    SRR = lapply(sorted$Sample_ID, function(x) {linkerFileSmall[grep(x,linkerFileSmall$Sample)]$Accession[1]})
  }
  if(is.na(SRR[1])){
    stop(cat("Could not find SRR ids for study: ",studyName,"\n"))
  }
  
  study = strsplit(as.character(sorted$Study[1])," ")[[1]][1]
  sorted$SRR = SRR
  
  sorted = checkExperimentNamesAreCorrect(sorted)
  
  
  ####Make rna seq and rfp paths
  currentRnaPath = paste0(rnaPath,grep(study,list.files(rnaPath),ignore.case = T,value = T),rnaRFPMiddlePathName)
  currentRfpPath = paste0(rfpPath,grep(study,list.files(rfpPath),ignore.case = T,value = T),rnaRFPMiddlePathName)
  
  folders = data.frame(matrix("A",  nrow = nrow(sorted),ncol = 1))
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

#'Supports only specific tissue conversions:
#'Hela -> Ovary, HEK293 ->kidney, THP-1 ->macrophage, PC3 -> prostate
getNameVariantsForCageTissue = function(sorted){
  if(length(unique(sorted$Tissue_or_CellLine)) != 1){
    stop("warning! different tissue in experiment, not allowed!")
  }
  cellLine = sorted$Tissue_or_CellLine[1]
  if(cellLine == "HeLa"){
    cellLine = "Ovary"
  }else if(cellLine == "HEK293" | cellLine == "HEK293T" ){
    cellLine = "kidney"
  }else if(cellLine == "THP-1"){
    cellLine = "macrophage"
  }else if(cellLine == "PC3"){
    cellLine = "prostate"
  }
  
  currentCageFiles = grep(cellLine,list.files(cageFiles),ignore.case = T,value = T)
  if(length(currentCageFiles) == 0)
    stop(cat("found no matches for cage type: ",as.character(cellLine),"\n"))
  
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
#'The acctual multithreaded run of the experiments
#'Need to optimize this!!!, it will stop using cores
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
      system(paste(bashScriptLocation,T,n,rnaName,rfpName,as.character(sorted$Tissue_or_CellLine[1])),wait = tooManyRunning) 
    }
  }
  cat("Running ",a," cores, dont do too many!!")
}

