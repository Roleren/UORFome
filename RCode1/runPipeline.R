# Parts:
#1. set up parameters
#2. Find new cage leaders
#3. Find uORFs
#4. Create feature database



##### First set up parameters
# 1. If output dir is not created yet, run:

# orfikDirs(mainPath = , makeDatabase = )

#set working dir correctly to ./RCode1/ location
setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./uorfomeGeneratorHelperFunctions.R")


# set up multithreading options
pipelineCluster()


##### Second Find new cage leaders


cageList = list.files(cageFolder) # make sure this folder is correct
nCageList = length(cageList)

getLeadersFromCage(nCageList)



##### Third Find uORFs

leadersList = list.files(leadersFolder)
nLeadersList = length(leadersList)
#clusterCall(cl, function(x) .libPaths(x), .libPaths())

getUorfsFromLeaders(nLeadersList)


### Fourth Create feature database
#4. 

setwd("/export/valenfs/projects/uORFome/RCode1/")
source("./DataBaseCreator.R")

createCatalogueDB() # fix this to work

# stop cluster
stopCluster(cl)
rm(cl)





# library(Rcpp)
# sourceCpp(code = '
#   #include <fstream>
#   #include <iostream>
#   #include <vector>
#   #include <sstream>
#   #include <algorithm>
#   #include <cassert>
#   
#   #include <Rcpp.h>
#   
#   
#   using namespace Rcpp;
#   
#   Function GRangesC("GRanges", Environment::namespace_env("GenomicRanges"));
#   Function IRangesC("IRanges", Environment::namespace_env("IRanges"));
#   // [[Rcpp::export]]
#   S4 findORFs_fasta(){
#     std::vector<int> start;
#     std::vector<int> stop;
#     std::vector<int> Seqnames;
#     std::vector<int> strands;
# 
#     for(int i = 0; i < 10; i++){
#       start.push_back(i);
#       stop.push_back(i);
#       Seqnames.push_back(1);
#       strands.push_back(-1);
#     }
# 
#     return (GRangesC(Seqnames,
#                    IRangesC(wrap(start),
#                             wrap(stop)),
#                             strands));
# 
#   }')

