# all features of Riboseq and RNAseq
allFeaturesAtlas <- function(){
  #validate ribo and rna
  validateRiboRNAPairs()

  # first sequence features
  getSequenceFeatures()

  # then computeFeatures
  getAllFeaturesFromUorfs()
  # then rfpTables
  riboAtlasFPKMTissue() # fix this is wrong

  # then rna tables
  getRNAFpkms()
  rnaAtlasFPKMTissue()
  # then do te tables
  getTeFeatures()
  teAtlasTissue()

  #cds and 3'
  getCDSFeatures()
  getFeaturesThreeUTRs()
}

#' Main function to fill uORF database
createCatalogueDB <- function(){
  dataBaseFolder <- p(mainFolder,"/dataBase")
  setwd(dataBaseFolder)
  uorfDB <- createDataBase("uorfCatalogue.sqlite")
  createUniqueIDs()
  createGRObjects()
  createUORFAtlas()
  getTissueTable()
  allFeaturesAtlas()
}
