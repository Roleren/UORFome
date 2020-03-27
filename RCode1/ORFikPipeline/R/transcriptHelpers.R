# Transcript helper functions

#' Get filtered transcript names by groups
#' @param splitList the list per element is character vector of transcripts
#' @param extension what is inpute name extension on parts, 1 per txdb
#' names = c("validMir", "validNames")
#' extensions = c("R", "B")
#' @export
filterTranscriptsSplit <- function(splitList, names = "txNames", extensions = "",
                                   l = 100, c = 100, t = 100, longestPerGene = FALSE) {
  x = 1
  for (e in extensions) {
    accepted <- filterTranscripts(get(paste0("txdb", e)), l, c, t, longestPerGene)
    assign(paste0(names[1], e), accepted[accepted %in% splitList[[x]]], envir = .GlobalEnv)
    if (length(splitList) * 2 == length(names)) {
      assign(paste0(names[2], e), accepted[!(accepted %in% splitList[[x]])], envir = .GlobalEnv)
    }
    x = x + 1
  }
}

#' Split regions into groups
#' Assigns to .GlobalEnv all variables needed
#' @param splitList the list per element is character vector of transcripts to keep
#' @param extension what is inpute name extension on parts
#' @param splitExt what is output name extension on parts
splitRegions <- function(parts = c("mrna", "leaders", "cds", "trailers"),
                         splitList, extensions = "", splitExt = "") {
  if (length(splitList) != length(splitExt)) stop("stop!")

  x = 1
  for (e in extensions) {
    for (s in splitExt) {
      for (i in parts) {
        part <- paste0(i, e)
        grl <- get(part)
        assign(x = paste0(part, s), value = grl[names(grl) %in% splitList[[x]]],
               envir = .GlobalEnv)
      }
      x = x + 1
    }
  }
  return(NULL)
}
