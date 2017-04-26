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
getMatrix = function(data = "/export/valenfs/projects/uORFome/test_results/Old_Tests/test_data",  doubleBAM = F, usingCage = F, cageName = "/export/valenfs/projects/uORFome/DATA/CAGE/human/brain%2c%20adult%2c%20donor1.CNhs11796.10084-102B3.hg38.nobarcode.ctss.bed.gz" ){
  infoPrints(data,doubleBAM,usingCage) #print info about run
  ####PRE LOADINGS####
  if(exists("RFP") == F){
    cat("starting loading objects\n")

    Gtf = makeTxDbFromGFF(findFF("gtf"))
    cds = cdsBy(Gtf,"tx",use.names = T)
    fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
    threeUTRs = threeUTRsByTranscript(Gtf,use.names = T)
    
    if(usingCage){
      cat("Using cage..")
      fiveUTRs = getNewfivePrimeUTRs(fiveUTRs,cageName)
      cat("finished new 5' UTRs")
    }

    
    rna = getRNAseq() #get rna seq file
    RFP = getRFP()# get ribozomal foot prints file
      
    
    if(saveToGlobal){
      print("saving objects to global")
      assign("rna",rna,envir = .GlobalEnv)
      assign("Gtf",Gtf,envir = .GlobalEnv)
      assign("RFP",RFP,envir = .GlobalEnv)
      assign("cds",cds,envir = .GlobalEnv)
      assign("fiveUTRs",fiveUTRs,envir = .GlobalEnv)
      assign("threeUTRs",threeUTRs,envir = .GlobalEnv)
      
      
      
    }
    cat("finished loading objects\n")
  }
  
  ###GET te values and save them in a matrix with all values
  
  ##get leader, cds, 3' TE's and lengths of them all
  getGeneralTEValues(usingCage)
  
  cat("started finding UORFS\n")
  if(UorfRangesNotExists())
    rangesOfuORFs = scanUORFs(fiveUTRs,saveToFile = T)
  
  cat("finding Te-UORFs\n")
  if(exists("teUORF") == F){
    teUORF = getTE(rangesOfuORFs,rna,RFP,allLengths$tx_len,sapply(rangesOfuORFs, function(x) sum(width(x))),specificTE = "teUORF")
    if(saveToGlobal)
      assign("teUORF",teUORF,envir = .GlobalEnv)
  }
  
  save.image("results.Rdata")
  plotting(matrix)
}

###Script starting point
##Either run with reference to input folder, or use standard folder UORFome/test_data on kjempetuja
saveToGlobal = T
startUORFomeGenerator(arcs)
