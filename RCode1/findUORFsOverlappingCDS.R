#arcsRFU = commandArgs(trailingOnly = T)

source("./createfasta.R")
library(GenomicFeatures)
library(GenomicAlignments)
require(data.table)


###Find all uorf into the cds,filter bad ones

# If called without path, need rangesofUORFs in global scope!!

removeFalseUORFs = function(rangesOfuORFs, loadPath = NULL,saveToFile = F,outputFastaAndBed = F, nameSave = NULL){
  
  ###########PRE LOADINGS#############
  print("starting to filter out bad ourfs...")
  if(!is.null(loadPath))
    load(loadPath,envir = .GlobalEnv)
  
  ################FILTERING############
  #rangesOfuORFs = rangesOfuORFs[width(rangesOfuORFs) > 5] #filter out 0
  
  #check start upstream, ending downstream
  rangesOfuORFs = removeUORFsThatAreAcctualyCDSs(rangesOfuORFs)
  
  
  print("finished filtering bad ourfs")
  
  ################SAVING###############
  ####Save overlaps of cds and uorfs for plotting later
  
  if(saveToFile){ #Not working!
    cat("Saving uorfs rdata to: ",loadPath)
    save(rangesOfuORFs, file = loadPath)
  }
  if(outputFastaAndBed){
    nameU = chooseUORFName(loadPath,nameSave)
    createFastaAndBedFile(rangesOfuORFs,loadPath = NULL,nameUsed = nameU)
  }
  return(rangesOfuORFs)
}
### Use findOverlaps to find equal start sites, this work for - strand ?
removeUORFsThatAreAcctualyCDSs = function(rangesOfuORFs){
 
  getGTF()
  getCDS()
  cdsUsed = unlist(cds[names(rangesOfuORFs)])
  gr = unlist(rangesOfuORFs,use.names = F)
  overlaps =findOverlaps(query = IRanges(start(gr),end = end(gr)) ,
                              subject = IRanges(start(cdsUsed),end=end(cdsUsed)), type = "equal")
  #remove and recombine to grangeslist
  gr = gr[-overlaps@from]
  
  l = Rle(names(gr))
  t = unlist(lapply(1:length(l@lengths),function(x){ rep(x,l@lengths[x])}))
  grl = split(gr,t)
  names(grl) = unique(names(gr))
  print("Removed uorfs that were acctualy cdss")
  return(grl)
}

saveOverlaps = function(rangesOfuORFs,cds){
  numberOfOverlaps = getUOrfOverlaps()
  assign("numberOfOverlaps",numberOfOverlaps,envir = .GlobalEnv)
}
