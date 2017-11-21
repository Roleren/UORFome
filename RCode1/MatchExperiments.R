#Current plan:
#test the new features with gonzales data

arcs = commandArgs(trailingOnly = T)


if(length(arcs) == 1){
  setwd(arcs[1])
}else{
  setwd("/export/valenfs/projects/uORFome/RCode1")
}
source("./MatchExperimentsHeader.R")




getLinkerFile()
filterBadStudies()

for(study in unique(SpeciesGroup$Study)){ #for each study, run it
  sorted = getSpecificStudyAndSpecies(study)
  
  currentCageFiles =  getCageFiles(sorted,speciesName)
  
  runExperiments(sorted,currentCageFiles)
  
}
