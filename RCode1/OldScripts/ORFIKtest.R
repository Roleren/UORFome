#library(ORFik)
require(testthat)
require(GenomicFeatures)
require(GenomicAlignments)
require("stringi")
require("IRanges")

#find first startcodon, find end in frame, go to next base, check for new startcodon etc..

find_in_frame_ORFs <- function(fastaSeq,
                               startCodon = "ATG",
                               stopCodon = "TAA|TAG|TGA",
                               longestORF = F,
                               minimumLength = 0) {
  longestORFRegex <- paste0("(?<!(", startCodon, "))")
  mainRegex <- paste0("(?=(", startCodon, "))(?=(?:([ATGCN]{3}))*?(", stopCodon, "))")
  codpos <- paste0(if (longestORF) longestORFRegex, mainRegex)
  a <- gregexpr(codpos, fastaSeq, perl = TRUE)[[1]]
  gr <- IRanges(start = as.vector(a), end = attr(a,"capture.start")[, if (longestORF) 4 else 3] + 2)
  
  gr <- gr[width(gr) >=  6 + minimumLength*3]
  return(gr[order(start(gr), end(gr))])
}


ge