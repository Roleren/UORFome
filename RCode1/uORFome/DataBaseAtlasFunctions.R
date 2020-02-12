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

createCatalogueDB <- function(){
  #uorfDB <- createDataBase(databaseName)
  createUniqueIDs()
  createGRObjects()
  createUORFAtlas()
  getTissueTable()
  allFeaturesAtlas()
}
