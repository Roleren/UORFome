###SCRIPT FOR FINDING UORFS in dataset and doing statistics on the result.
##Plot are made and output files saved to same folder as input data
##All input data must be in the same folder!

###########INPUTS################
#1. folder location with all the files(gtf, bed and bam files)(you can also include rangesOfUORFs)
#2. Using 2 bam files = T, using 1 bam and 1 bed = F
#3. Using cage data for improved 5prime utr recognition, T = using. 


arcs = commandArgs(trailingOnly = T)
source("/export/valenfs/projects/uORFome/RCode1/uorfomeGeneratorHelperFunctions.R")


###MAIN FUNCTION###
getMatrix = function(data = dataFolder, leaderBed = NULL,  doubleBAM = F, usingNewCage = F, cageName = standardCage ,rnaSeq = NULL,rfpSeq = NULL){
  infoPrints(data,doubleBAM,usingNewCage,cageName,leaderBed,rnaSeq,rfpSeq) #print info about run
  ###PRE LOADINGS####
  if(exists("RFP") == F){
    
    getGTF()
    getCDS()
    
    getLeaders(leaderBed,usingNewCage,cageName) #get five prime utrs with or without cage
    threeUTRs = threeUTRsByTranscript(Gtf,use.names = T)
    
    getRNAseq(rnaSeq) #get rna seq file
    getRFP(rfpSeq)# get ribozomal foot prints file
      
    print("saving objects to global")
    
    assign("cds",cds,envir = .GlobalEnv)
    assign("threeUTRs",threeUTRs,envir = .GlobalEnv)
    
    cat("finished loading objects\n")
  }
  
  ###GET te values and save them in a matrix with all values
  
  ##get leader, cds, 3' TE's and lengths of them all
  getGeneralTEValues(usingNewCage,leaderBed)
  
  decideHowToGetUORFRanges() #Load ranges of Uorfs
  
  cat("finding Te-UORFs\n")
  if(exists("teUORF") == F){
    teUORF = getTE(rangesOfuORFs,rna,RFP,allLengths$tx_len,sapply(rangesOfuORFs, function(x) sum(width(x))),specificTE = "teUORF")
    
    assign("teUORF",teUORF,envir = .GlobalEnv)
  }
  saveRData() #Save RData of project
  plotting(matrix,paste0(plottingFolder,detailedFullName,".pdf")) #plot results
}

###Script starting point
##Either run with reference to input folder, or use standard folder UORFome/test_data on kjempetuja
startUORFomeGenerator(arcs)
