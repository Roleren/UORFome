# Load data needed

dfC <- getMcManus17()
outputLibs(dfC)
bamVars <- bamVarName(dfC)
CAGE <- GRanges()

for (v in bamVars) {
  CAGE <- c(CAGE, convertToOneBasedRanges(get(v), addScoreColumn = T))
}
export.bed(CAGE, paste0(dirname(dfC$filepath[1]),"/merged_CAGE_5prime.bed"))

dfC <- getWery15()
outputLibs(dfC)
bamVars <- bamVarName(dfC)
CAGE <- GRanges()

for (v in bamVars) {
  CAGE <- c(CAGE, convertToOneBasedRanges(get(v), addScoreColumn = T))
}
export.bed(CAGE, paste0(dirname(dfC$filepath[1]),"/merged_CAGE_5prime.bed"))
