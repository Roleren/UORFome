setwd("/export/valenfs/projects/uORFome/test_results/rangesOfUORFs/")
sizes = list()
index = 1
for(file in list.files()){
  load(file)
  nameUsed = gsub("%.*","",file)
  
  sizes[index] = length(unique(unlist(rangesOfuORFs)[,"names"]))
  names(sizes)[index] = nameUsed 
  index = index+1
}
averageValue = mean(unlist(sizes))


# create unique identifiers
# filter uorfs, cds
# Integration, how should we place and name files ?
#type 1 2 3
#brain Cage Ribo rna