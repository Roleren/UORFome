#########Variables used by uorfomeGenerator############
#Change them if needed on other folders!
#####SLASH IS ALWAYS ADDED IN START; NEVER END!!!!

p = function(nameA,nameB){
  paste0(nameA,nameB)
}

mainFolder = "./.."

resultsFolder = p(mainFolder,"/test_results") #output folder
helperMainFolder = p(resultsFolder,"/Old_Tests/test_data") #main folder for gtf,fasta and .fai
dataMainFolder = p(mainFolder,"/DATA") #data folder for cage, rna seq and rfp
codeFolder = p(mainFolder,"/RCode1")

dataFolder = helperMainFolder
fastaName = p(dataFolder,"/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
faiName = p(dataFolder,"/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
gtfName = p(dataFolder,"/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf")

cageFolder = p(dataMainFolder,"/CAGE/human")
standardCage = p(cageFolder,"/brain%2c%20adult%2c%20donor1.CNhs11796.10084-102B3.hg38.nobarcode.ctss.bed.gz")

matrixFolder = p(resultsFolder,"/Matrices/")
RdataFolder = p(resultsFolder,"/Rdata/")
plottingFolder = p(resultsFolder,"/Plotting/Single result Plots")
leadersFolder = p(resultsFolder,"/New_Cage_Leaders")
fastaFolder = p(resultsFolder,"/fasta")
uorfBedFolder = p(resultsFolder,"/bedUORFS")
bamFolder = p(resultsFolder,"/SortedAndIndexedBams")


