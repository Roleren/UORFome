#########Variables used by uorfome data-base############
#Change them if needed on other folders!
#####SLASH IS ALWAYS ADDED IN START;

# input folders
source("./uORFome/Init_Variables.R")

# output folders
plottingFolder = p(resultsFolder,"/Plotting/Single_result_Plots/")
regionUORFsFolder = p(resultsFolder,"/regionUORFs/")
leadersFolder = p(resultsFolder,"/New_Cage_Leaders/")
fastaFolder = p(resultsFolder,"/fasta/")
uorfFolder = p(resultsFolder,"/rangesOfUORFs/")
idFolder = p(resultsFolder,"/uorfIDs/")

if (!dir.exists(leadersFolder)) stop("could not find results folders, run orfikDirs()")

#debug/test variables
# standardCage <- p(cageFolder,"brain%2c%20adult%2c%20donor1.CNhs11796.10084-102B3.hg38.nobarcode.ctss.bed.gz")
# standardRFP <- "/export/valenfs/data/processed_data/Ribo-seq/fantom_human_bed/per_length/merged/Andreev_DE_2015.Human.HEK293.RPF.GRCh38.SRR1173914.reads_merged.bed"
# standardRNA <-p(rnaFolder, "gonzalez_C_2014_human_mouse/final_results/aligned_GRCh38/Gonzalez_C_2014.Human.brain.RNA.GRCh38.SRR1562544.bam")