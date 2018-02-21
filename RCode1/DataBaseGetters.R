
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