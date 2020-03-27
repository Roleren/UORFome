### This is the main script to run to get the features

### For many data-sets use a server, to speed up.


setwd("/export/valenfs/projects/uORFome/RCode1/uORFome/")
source("./DataBaseSetup.R")

pipelineCluster(10) # number of cores

rfpList <- grep(pattern = "merged",x = list.files(rfpFolder), value = T)
nrfpList <- length(rfpList)
# all rfp features
foreach(i=1:nrfpList, .inorder = F) %dopar% {
  dataBaseFolder <- "/export/valenfs/projects/uORFome/dataBase"
  if (file.exists(paste0( dataBaseFolder, "/featureTablesTemp/dt_",i,".rdata")))
    return(print(i))
  setwd("/export/valenfs/projects/uORFome/RCode1/uORFome/")
  source("./DataBaseSetup.R")
  
  rfpList <- grep(pattern = "merged",
                  x = list.files(rfpFolder), value = T)
  RFPPath <- p(rfpFolder, rfpList[i])
  
  getAllFeatures(grl = NULL, RFPPath, i = i)
  return(print(i))
}

# list <- foreach(i=1:nrfpList, .inorder = F) %do% {
#   dataBaseFolder <- "/export/valenfs/projects/uORFome/dataBase"
#   if (!file.exists(paste0( dataBaseFolder, "/featureTablesTemp/dt_",i,".rdata")))
#     return(i)
# 
#   return(0)
# }
# list <- unlist(list)
# setwd("/export/valenfs/projects/uORFome/RCode1/")
# source("./DataBaseCreator.R")
# setwd("/export/valenfs/projects/uORFome/dataBase/")
# 
# uorfID <- getORFNamesDB(T, F, F)
# floss <- uorfID
# entropyRFP <- uorfID        
# disengagementScores <- uorfID
# RRS <- uorfID                
# RSS <- uorfID           
# fpkmRFP <- uorfID           
# ORFScores <- uorfID
# ioScore <- uorfID
# #rm(list = c("floss","entropyRFP","disengagementScores", "RRS", "RSS", "fpkmRFP", "ORFScores", "ioScore" ))
# 
# # all dt features, split them, save them seperatly
# foreach(i=1:nrfpList) %do% {
#   
#   dtFolder <- "featureTablesTemp/"
#   dtPath <- paste0(dtFolder, "dt_",i,".rdata")
#   
#   load(dtPath)
#   if(ncol(dt) != 8) stop(paste("not correct ncol of", i))
#   floss[, p("floss_",i)] <- dt$floss
#   entropyRFP[, p("entropyRFP_",i)] <- dt$entropyRFP
#   disengagementScores[, p("disengagementScores_",i)] <- dt$disengagementScores 
#   RRS[, p("RRS_",i)] <- dt$RRS
#   RSS[, p("RSS_",i)] <- dt$RSS
#   fpkmRFP[, p("fpkmRFP_",i)] <- dt$fpkmRFP
#   ORFScores[, p("ORFScores_",i)] <- dt$ORFScores
#   ioScore[, p("ioScore_",i)] <- dt$ioScore
#   print(i)
# }
# # insert all the ribo features tables
# insertTable(floss, "floss")
# insertTable(entropyRFP, "entropyRFP")
# insertTable(disengagementScores, "disengagementScores")
# insertTable(RRS, "RRS")
# insertTable(RSS, "RSS")
# insertTable(fpkmRFP, "Ribofpkm")
# insertTable(ORFScores, "ORFScores")
# insertTable(ioScore, "ioScore")
# 
# # insert info
# getRiboInfoTable(rfpList = rfpList)

stopCluster(cl)
rm(cl)