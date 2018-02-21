#########Variables used by uorfome data-base############
#Change them if needed on other folders!
#####SLASH IS ALWAYS ADDED IN START; maybe end, think about this!!!

p = function(nameA,nameB){
  paste0(nameA,nameB)
}

mainFolder = "./.." # main folder is one back from RCode1/ folder
#primary folders
resultsFolder = p(mainFolder,"/test_results") #output folder
helperMainFolder = p(resultsFolder,"/Old_Tests/test_data") #main folder for gtf,fasta and .fai
dataMainFolder = p(mainFolder,"/DATA") #data folder for cage, rna seq and rfp
codeFolder = p(mainFolder,"/RCode1")
#data folders
dataFolder = helperMainFolder
fastaName = p(dataFolder,"/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
faiName = p(dataFolder,"/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
gtfName = p(dataFolder,"/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf")
gtfdb = p(dataFolder,"/Gtf.db")  ### This is wrong!!!

#output folders
matrixFolder = p(resultsFolder,"/Matrices/")
RdataFolder = p(resultsFolder,"/Rdata/")
plottingFolder = p(resultsFolder,"/Plotting/Single_result_Plots/")
leadersbedFolder = p(resultsFolder,"/New_Cage_bedLeaders/")
leadersFolder = p(resultsFolder,"/New_Cage_Leaders/")
fastaFolder = p(resultsFolder,"/fasta/")
uorfBedFolder = p(resultsFolder,"/bedUORFS/")
uorfFolder = p(resultsFolder,"/rangesOfUORFs/")
bamFolder = p(resultsFolder,"/SortedAndIndexedBams/")
idFolder = p(resultsFolder,"/uorfIDs/")

# input folders (cage, ribo and rna)
cageFolder = p(dataMainFolder,"/CAGE/human/")
rfpFolder <- "/export/valenfs/data/processed_data/Ribo-seq/fantom_human_bed/per_length/merged/"
rnaFolder <- "/export/valenfs/data/processed_data/RNA-seq/"
#specific run names
if(exists("detailedFullName") == F){
  detailedFullName = ""
  thisCage = NULL
}


#debug/test variables
standardCage <- p(cageFolder,"brain%2c%20adult%2c%20donor1.CNhs11796.10084-102B3.hg38.nobarcode.ctss.bed.gz")
standardRFP <- p(rfpFolder,"Andreev_DE_2015.Human.HEK293.RPF.GRCh38.SRR1173905.reads_merged.bed")
standardRNA <-p(rnaFolder, "gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/Gonzalez_C_2014.Human.brain.RNA.GRCh38.SRR1562544.bam")