#########Variables used by uorfome data-base############
#Change them if needed on other folders!
#####SLASH IS ALWAYS ADDED IN START; maybe end, think about this!!!

source("./Init_Variables.R")

mainFolder = "./.." # main folder is one back from RCode1/ folder
#primary folders
resultsFolder = p(mainFolder,"/test_results") #output folder

dataMainFolder = p(mainFolder,"/DATA") #data folder for cage, rna seq and rfp
codeFolder = p(mainFolder,"/RCode1")
#data folders

fastaName = p(dataFolder,"/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")

# speed up data
gtfdb = p(dataFolder,"/Gtf.db")  ### a speed up for Gtf
# make the ones for fiveUTRs, cds, threeUTRs and tx

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


#specific run names
if(exists("detailedFullName") == F){
  detailedFullName = ""
  thisCage = NULL
}

#debug/test variables
standardCage <- p(cageFolder,"brain%2c%20adult%2c%20donor1.CNhs11796.10084-102B3.hg38.nobarcode.ctss.bed.gz")
standardRFP <- "/export/valenfs/data/processed_data/Ribo-seq/fantom_human_bed/per_length/merged/Andreev_DE_2015.Human.HEK293.RPF.GRCh38.SRR1173914.reads_merged.bed"
standardRNA <-p(rnaFolder, "gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/Gonzalez_C_2014.Human.brain.RNA.GRCh38.SRR1562544.bam")