
getCageInfoTable <- function(){
  if (tableNotExists("cageInformation")){
    require(xlsx)
    cageTable <- read.xlsx("../HumanSamples2.0.sdrf.xlsx", sheetName = "Sheet1")
    #filter bed tissues
    cageTable <- as.data.table(cageTable)
    cageTable[is.na(Characteristics.Tissue.)] <- "unclassifiable"
    cageTable[Characteristics.Tissue. == "Urethra"]$Characteristics.Tissue. <- "urethra" 
    cageTable[Characteristics.Tissue. == "adipose tissue"]$Characteristics.Tissue. <-"adipose"
    
    insertTable(Matrix = cageTable,tableName = "cageInformation")
  } else{
    cageTable <- readTable("cageInformation")
  }
  return(cageTable)
}

#' get ribo-seq info table
#' 
#' If not already made, make it and save it
#' @param name what should the table be called
#' @param rfpList a character vector of relative paths
getRiboInfoTable <- function(tableName = "RiboSeqInfo", rfpList){
  if (tableNotExists(tableName)){
    library(splitstackshape)
    library(stringr)
    dtRFPList <- data.table(rfpList)
    
    colNames <- c("Study", "Species", "Tissue/Cell_line", "Seq-Tech",
                  "Annotation", "SRR", "file_type", "link")
    
    splitList <- splitstackshape::cSplit(indt = dtRFPList,
      stripWhite = F, sep = ".", splitCols = "rfpList")
    
    #cheat: fix ruthkowski:
    ruths <- splitList[splitList$rfpList_01 == "Rutkowski_AJ_2015",]
    ruths[,5:8] <- ruths[,c(8,10:12)]
    splitList[splitList$rfpList_01 == "Rutkowski_AJ_2015",] <- ruths
    splitList <- splitList[, c(1:6, 8)]
    splitList[, "full_name"] <- rfpList
    names(splitList) <- colNames
    
    insertTable(Matrix = splitList,tableName = tableName)
  } else{
    splitList <- readTable(tableName)
  }
  return(splitList)
}

#' get rna-seq info table
#' 
#' If not already made, make it and save it
#' @param name what should the table be called
#' @param rfpList a character vector of relative paths
getRNASeqInfoTable <- function(tableName = "RNASeqInfo", rnaList){
  if (tableNotExists(tableName)){
    library(splitstackshape)
    library(stringr)
    rnaList <- getRelativePathName(rnaList)
    dtRNAList <- data.table(rnaList)
    
    colNames <- c("Study", "Species", "Tissue/Cell_line", "Seq-Tech",
                  "Annotation", "SRR", "file_type", "link")
    
    splitList <- splitstackshape::cSplit(indt = dtRNAList,
                                         stripWhite = F, sep = ".", splitCols = "rnaList")
    
    #cheat: fix ruthkowski:
    ruths <- splitList[splitList$rnaList_01 == "Rutkowski_AJ_2015",]
    ruths[, 5:7] <- ruths[, 11:13]
    splitList[splitList$rnaList_01 == "Rutkowski_AJ_2015",] <- ruths
    splitList <- splitList[, 1:7]
    splitList[, "full_name"] <- rnaList
    names(splitList) <- colNames
    
    insertTable(Matrix = splitList,tableName = tableName)
  } else{
    splitList <- readTable(tableName)
  }
  return(splitList)
}

#' Match ribo-seq and rna-seq by SRR numbers
#' 
#' If you have downloaded a Sequence Read Archiv
#'  from ncbi, this can help you avoid a lot
#'  of manual matching of experiments
#'  The experiments must have SRR numbers
#'  
#'  Uorfome used:
#'  --2017-05-10 16:47:12--  
#'  ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
#'  => “SRA_Accessions.tab”
matchRNA_RFPInfo <- function(tableName = "linkRnaRfp"){
  
  if(!tableNotExists(tableName)){
    return(readTable(tableName))
  }
  
  getMatchingTable()
  filtered <- SpeciesGroup[SpeciesGroup$Species == "Human",]
  # pre filter, after we know how
  
  rfpInfo <- readTable("RiboSeqInfo")
  rnaInfo <- readTable("RNASeqInfo")
  rfpInfoTemp <- rfpInfo
  rnaInfoTemp <- rnaInfo
  # per study
  rnaStudies <- unique(rnaInfo$Study)
  rfpStudies <- unique(rfpInfo$Study)
  
  if(length(rnaStudies) <= length(rfpStudies)){
    rnaStudies <- rnaStudies[rnaStudies %in% rfpStudies]
    rfpStudies <- rfpStudies[rfpStudies %in% rnaStudies]
  } else {
    rfpStudies <- rfpStudies[rfpStudies %in% rnaStudies]
    rnaStudies <- rnaStudies[rnaStudies %in% rfpStudies]
  }
  rfpInfo <- rfpInfo[rfpInfo$Study %in% rnaStudies,]
  rnaInfo <- rnaInfo[rnaInfo$Study %in% rnaStudies,]
  
  if ((length(rnaStudies) == 0) || (length(rfpStudies) == 0)) {
    stop("got 0 matching rna-seq and ribo-seq experiments, check table")
  }
  filtered$Study <- gsub(pattern = ",", replacement = "_", filtered$Study)
  filtered$Study <- gsub(pattern = " ", replacement = "_", filtered$Study)
  realTable <- filtered[filtered$SRR %in% c(rfpInfo$SRR, rnaInfo$SRR),]
  
  # Now here write your integer vector, i.g. c(1,1,2,2,3,3,NA,4,4) etc.
  matching <- c(1,1,2,2,3,3,4,4,NA,  5,6,5,6,  7,NA,8,9,NA,NA,10,11, NA, NA, 7, NA,8,9,10,11,  12,13,12,13,14,15,14,15,  16,17,18,19,20,21,22,23,16,17,18,19,20,21,22,23,  24,25,26,27,28,29,24,25,26,27,28,29, 30,31,32,33,34,35,30,31,32,33,34,35)
  
  realTable$matching <- matching
  View(cbind(as.character(realTable$Sample_description), realTable$matching))
  
  # remove NAs
  realTable <- realTable[!is.na(realTable$matching),]
  rfpInfo <- rfpInfo[rfpInfo$SRR %in% realTable$SRR]
  rnaInfo <- rnaInfo[rnaInfo$SRR %in% realTable$SRR]
  if (nrow(rfpInfo) != nrow(rnaInfo)) {
    stop("nrow rna and nrow rfp is not equal, check again")
  }
  
  originalIndex <- sapply(unlist(1:length(realTable$SRR), use.names = F), function(x){
    if(realTable$Sample_Type[x] == "RNA"){
      return(which(rnaInfoTemp$SRR == realTable$SRR[x]))
    } else {
      return(which(rfpInfoTemp$SRR == realTable$SRR[x]))
    }
  })
  realTable$originalIndex <- originalIndex
  realTable <- as.data.table(realTable)
  
  
  # if this doesnt work, some of the rows contain empty fields, cleanup
  realTable[,RnaRfpFolders:=NULL]
  realTable[Study == "Gonzalez_C_2014"]$Study_title <- "Ribosome profiling reveals a cell-type-specific translational landscape in brain tumors."
  realTable$SRR <- unlist(realTable$SRR, use.names = F)
  
  insertTable(Matrix = realTable, tableName = tableName)
  return(realTable)
}