# ORFik QC report
# Welcome, update the steps you need to update, marked with 1., 2., 3. and 4
# Then run the ORFikQC function
#devtools::install_github("Roleren/ORFik@RELEASE_3_10")
rm(list=ls())
source("/export/valenfs/projects/uORFome/RCode1/loadUorfome.R")

# First you must create your ORFik experiment, 5 steps:
# 1. Update path to experiment data
exp_dir = "/export/valenfs/data/processed_data/Ribo-seq/Beaudoin_2018_zebrafish/aligned_GRCz10/merged"
# 2. Set a 5 character name for experiment
exper_name = "Beu18"
# 3. Create a template experiment, and check that it is correct
temp <- create.experimentl(dir = exp_dir, exper_name) # Make sure each row(sample) is unique
temp[5, 4] <- "PatA" # Did not find a way to split 24h samples
save.experimentl(temp)
# 5. Load experiment, and validate that it is correct
df <- read.experimentl(exper_name)

ORFikQC(df, out.dir = "~")
