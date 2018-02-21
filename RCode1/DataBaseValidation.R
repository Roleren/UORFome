
#' Validate uorfs of data-base
#' 
#' This is a check to see that pipeline have done everything correctly
#' if redoing the findOverlaps does not find all orfs within fiveUTRs
#' it means that some orfs are outside the mapping area
#' this should not happen!
validateExperiments <- function(grl){
  
  fiveUTRs <- leaderAllSpanning()
  a <- findOverlaps(query = unlist(grl, use.names = F), fiveUTRs)
  a <- a[!duplicated(from(a))]
  if(length(a) != length(unlist(grl))){ 
    stop("Not all orfs where within the FiveUTRs used
         to make them, something is wrong!")
  } else { print("experiments look valid")}
}

