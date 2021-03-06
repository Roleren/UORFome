#########Variables used by uorfome data-base############
#Change them if needed on other folders!
#####SLASH IS ALWAYS ADDED IN START#####

# set these 4 directories
mainFolder = "/export/valenfs/projects/uORFome" # main folder is one back from RCode1/ folder

codeFolder = p(mainFolder,"/RCode1") # the rcode location

resultsFolder = p(mainFolder,"/results") #output folder

dataFolder = p(mainFolder,"/Annotations") #location of gtf, fasta and .fai
faiName = p(dataFolder,"/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
gtfName = p(dataFolder,"/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf")
gtfdb = p(dataFolder,"/Gtf.db")  ### a speed up for Gtf, remove if not used

# now validate all that directories exist
if(!all(dir.exists(c(codeFolder, resultsFolder, dataFolder)))){
  stop(p("Could not find directory: ", c(codeFolder, resultsFolder, dataFolder)[!file.exists(c(codeFolder, resultsFolder, dataFolder))]))
}

# now validate all files exist
if(!all(file.exists(c(faiName, gtfName, gtfdb)))){
  stop(p("Could not find file: ", c(faiName, gtfName, gtfdb)[!file.exists(c(faiName, gtfName, gtfdb))]))
}

# input folders (cage, ribo and rna) Optional (set to NULL if not needed)
cageFolder = "/export/valenfs/projects/uORFome/DATA/CAGE/human/"
rfpFolder <- "/export/valenfs/data/processed_data/Ribo-seq/fantom_human_bed/per_length/merged/"
rnaFolder <- "/export/valenfs/data/processed_data/RNA-seq/"

if(!is.null(cageFolder)){
  if(!dir.exists(cageFolder)){
    stop("cage folder not found")
  }
}

if(!is.null(rfpFolder)){
  if(!dir.exists(rfpFolder)){
    stop("ribo-seq folder not found")
  }
}

if(!is.null(rfpFolder)){
  if(!dir.exists(rfpFolder)){
    stop("rna-seq folder not found")
  }
}

# results

#########Variables used by uorfome data-base############
#Change them if needed on other folders!
#####SLASH IS ALWAYS ADDED IN START;

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
