###SCRIPT FOR FINDING UORFS in dataset and doing statistics on the result.
##Plot are made and output files saved to same folder as input data
##All input data must be in the same folder!

###########OVERVIEW / INFO################
#1. folder location with all the files(gtf, bed and bam files)(you can also include rangesOfUORFs)
#2. Using 2 bam files = T, using 1 bam and 1 bed = F
#3. Using cage data for improved 5prime utr recognition, T = using. 


#test data:
#Sidrauski_C_2015.Human.HEK293.RPF.GRCh38.SRR1795435.bam
#Andreev_DE_2015.Human.HEK293.RNA.GRCh38.SRR1173915.bam


#################INPUT READ FROM SHELL#############
arcs <- commandArgs(trailingOnly = T)
lArcs <- length(arcs)
if(lArcs == 5){
  # setwd(arcs[1])
  setwd("/export/valenfs/projects/uORFome/RCode1/")
}else if(lArcs == 1){ #Load from saved session
  load(file = paste0(RdataFolder,arcs[1]))
  setwd("/export/valenfs/projects/uORFome/RCode1/")
}else
  setwd("/export/valenfs/projects/uORFome/RCode1/")

if(lArcs != 1) #if using preload, dont do this!
  source("./uorfomeGeneratorHelperFunctions.R")

################END INPUT READ#####################


###MAIN FUNCTION###
getMatrix <- function(leaderBed = NULL,  doubleBAM = F, usingNewCage = F, cageName = standardCage ,rnaSeq = NULL,rfpSeq = NULL,doPreLoadings = T, tissueUsed = NULL){
  assign("tissueUsed",tissueUsed,envir = .GlobalEnv)
  infoPrints(doubleBAM,usingNewCage,cageName,leaderBed,rnaSeq,rfpSeq) #print info about run
  ##################PRE LOADINGS##############
  if( (exists("RFP") == F) & doPreLoadings ){
    
    getGTF()
    getCDS()
    
    getLeaders(leaderBed,usingNewCage,cageName) #get five prime utrs with or without cage
    getThreeUTRs()
    
    getRNAseq(rnaSeq) #get rna seq file
    getRFP(rfpSeq)# get ribozomal foot prints file
    
    cat("finished loading objects\n")
  }
  ##############FINISHED PRE LOADINGS#########
  
  ###GET te values and save them in a matrix with all values
  ##get 5', uorf, cds and 3' TE's and lengths of them all
  getGeneralTEValues(usingNewCage,leaderBed)
  
  decideHowToGetUORFRanges(assignUorf = T, thisCage) #Load ranges of Uorfs
  
  cat("finding Te-UORFs\n")
  if(exists("teUORF") == F){
    teUORF <- getTE(rangesOfuORFs,rna,RFP,allLengths$tx_len,sapply(rangesOfuORFs, function(x) sum(width(x))),specificTE = "teUORF")
    
    assign("teUORF",teUORF,envir = .GlobalEnv)
  }
  if(lArcs != 1)#if presaved, dont save again, fix this!!!!
    saveRData() #Save RData of project
  
  plotting(matrix,paste0(plottingFolder,gsub("%","_",detailedFullName),".pdf")) #plot results
}

###Script starting point
##Either run with reference to input folder, or use standard folder UORFome/test_data on kjempetuja

startUORFomeGenerator(arcs)
