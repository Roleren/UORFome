setwd("/export/valenfs/projects/uORFome/RCode1/")

pipelineCluster()

# load RNA list, this varies how you do, but get a characterVector of
# General paths
getRNAFpkms()

stopCluster(cl)
rm(cl)
