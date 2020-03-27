#' Get difference between wild type and target variant
#' Per library type like RNA-seq or Ribo-seq.
#' @param dt a data.table of counts or fpkm etc.
#' @return a data.table with differences
SEdif <- function(dt, df) {
  dif <- data.table()
  bamVars <- colnames(dt)
  
  bamVarsT <- bamVarName(df, skip.experiment = T, skip.condition = T)
  dists <- which(!is.na(chmatch(bamVarsT, bamVarsT[1])))
  if (length(dists) != 2) stop("Wrong naming!")
  seperator <- dists[2] - dists[1]
  
  for(i in 1:(ncol(dt)/2)) {
    a <- dt[, log2((get(bamVars[(i+seperator)]) + 0.00001) / (get(bamVars[(i)]) + 0.00001))] 
    dif <- cbind(dif, a)
  }
  colnames(dif) <- bamVarsT[1:(ncol(dt)/2)]
  return(dif)
}

SESplit <- function(dif, splitCol, score = 15, dtS, df) {
  whichMatch <- rowSums(dtS > score) == ncol(dtS)
  
  vdif <- dif[whichMatch, ]
  ismiRNA <- splitCol[whichMatch]
  
  vdifm <- melt(vdif)
  vdifm$ismiRNA <- rep(ismiRNA, ncol(dif))
  vdifm$stage <- gsub(".*_", x = vdifm$variable, replacement = "")
  
  if(!is.null(df$experiment)) {
    vdifm$type <- sub("_.*", sub("..._", x = vdifm$variable,
                                 replacement = "", perl = T), replacement = "")
  } else {
    vdifm$type <- ORFik:::libraryTypes(vdifm$variable)
  }
  
  d <- data.table()
  libTypes <- ORFik:::libraryTypes(df)
  d <- vdifm[type == libTypes[1],]
  
  for(i in libTypes[-1]) {
    d <- cbind(d, vdifm$value[vdifm$type == i])
  }
  
  
  colnames(d) <- c(c("variable", libTypes[1], "ismiRNA", "stage", "type"), libTypes[-1])
  return(d)
}
