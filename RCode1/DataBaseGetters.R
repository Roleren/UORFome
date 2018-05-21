
#' Get orf names from the orf data-base
#' 
#' This is the primary key for most tables in the data-base
#' @param with.transcript a logical(F), should the transcript be included, this makes the
#'  list have duplicated orfs
#' @param only.transcripts a logical(F), should only the transcript and not orfId be included
#' @param asCharacter a logical(T), should it return as character or data.table(F)
getORFNamesDB <- function(with.transcript = F, only.transcripts = F, asCharacter = T){
  if (with.transcript) {
    dt <- readTable("linkORFsToTx")
    if(only.transcripts){
      if (asCharacter) {
        return(as.character(unlist(dt[, 2], use.names = F)))
      }
      return(dt[, 2])
    }
    return(dt)    
  }
  
  dt <-readTable("uniqueIDs")
  if (asCharacter) {
    dt <- as.character(unlist(dt, use.names = F))
  }
 return(dt)
}

#' get GrangesList from uniqueIDs as strings
#' 
uniqueIdsAsGR <- function(makeBed = T){
  if(tableNotExists("SplittedByExonsuniqueUORFs")){
    uniqueIDs <- readTable("uniqueIDs")
    if (sum(duplicated(uniqueIDs)) > 0 ) stop("duplicated uorf names in uniqueIDs")
    grl <- toGRFromUniqueID(uniqueIDs)
    insertTable(Matrix = grl,tableName = "SplittedByExonsuniqueUORFs", rmOld = T)
    #grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
    
    #now make uscs bed 12 format
    if (makeBed)
      bed12(grl, "bedUniqueUorfs.bed", T)
    
  } else {
    grl <- readTable("SplittedByExonsuniqueUORFs", asGR = T)
  }
  return(grl)
}

# get table matching ribo-seq and rna-seq
getMatchingTable <- function(){
  fimport("/export/valenfs/projects/uORFome/test_results/
          Old_Tests/test_data/unfilteredSpeciesGroup.rdata")
}

#' Takes two tables from the database and extracts the rows of toBeMatched
#' that matches the txNames in referenced.
#' Both must have a column called txNames
#' @return the toBeMatched object matched by txNames
matchByTranscript <- function(toBeMatched, referenced){
  
  Indices <- data.table(txNames = toBeMatched$txNames, ind = 1:length(toBeMatched$txNames))
  merged <- merge(Indices, data.table(txNames = referenced$txNames),
                  by = "txNames", all.y = T, sort = F) 
  return(toBeMatched[merged$ind, ])
}

getIDColumns <- function(dt, allowNull = F){
  nIDs <- 0
  if (!is.numeric(dt[1,1][[1]])) {
    nIDs = nIDs + 1
    if (!is.numeric(dt[1,2][[1]])) {
      nIDs = nIDs + 1
    }
  }
  if(!nIDs){
    if (allowNull) {
      return(NULL)
    } else {
      stop("No id columns found for dt")
    }
  }
  return(dt[, nIDs, with = FALSE])
}
#' fix this to work on string tables
removeIDColumns <- function(dt){
  if (!is.numeric(dt[1,1][[1]])) {
    dt <- dt[, -1]
    if (!is.numeric(dt[1,1][[1]])) {
      dt <- dt[, -1]
    }
  }
  return(dt)
}

#' get the uorfs in the database
#' @param withExons should the uorfs be splitted by exons
#' @param withTranscripts should the uorfs have transcript information, 
#' warning, this will duplicate some uorfs.
#' @return a GRangesList or data.table, if(F, F)
getUorfsInDb <- function(withExons = T, withTranscripts = T){
  if (withExons && withTranscripts) {
    if(file.exists(p(dataBaseFolder, "/uORFsAsGR.rdata"))) {
      load(p(dataBaseFolder, "/uORFsAsGR.rdata"))
      return(grl)
    } else if(!tableNotExists("uorfsAsGRWithTx")) {
      grl <- readTable("uorfsAsGRWithTx", asGR = T)
      gr <- unlist(grl, use.names = F)
      names(gr) <- gsub("_[0-9]*", "", names(gr))
    } else {
      stop("uORFs could not be found, check that they exist")
    }
    return(groupGRangesBy(gr, gr$names))
  } else if (!withExons) {
    return(readTable("uniqueIDs"))
  } else if (withExons && !withTranscripts) {
    return(readTable("SplittedByExonsuniqueUORFs", asGR = T))
  } else {
    stop("not supported way of getting uorfs")
  }
}