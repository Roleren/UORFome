####FOR FUNCTION, might need to remove nan in NORMRNA######!!

library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(Biostrings)

library(ggplot2)
#setwd("/export/valenfs/projects/uORFome/RProject/")
setwd("/export/valenfs/projects/uORFome/")
# # if data is already loaded
if(exists("fasta") == F){
  
  #1. create vaiables for rna gtf and bed
  
  #mRna alligned to the transcriptome
  rna = readGAlignmentPairs("test_data/total_RNA.mock.Rep1.GRCh38.bam")
  #Gene description
  Gtf = makeTxDbFromGFF("test_data/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf","gtf")
  #Ribosome protected footprint
  RFP = import.bed("test_data/RP_mock_Rep1.all_runs.bam_shifted_footprints.bed")
  #Fasta file with sequences readDNAStringSet("D:/Hakon/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
  fasta = readDNAStringSet("test_data/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
  
  #2. get reads
  
  #Coding sequences
  cds = cdsBy(Gtf,"tx",use.names = T)
  #
  allLengths = transcriptLengths(Gtf,with.cds_len = T,with.utr5_len = T,with.utr3_len = T)
  
  #Find length of the coding sequences
  cdsLengths = allLengths$cds_len
  
  
  
  #Count overlaps between coding sequences and the mRna
  cdsRnaOverlaps = countOverlaps(cds,rna)
  
  # Count overlaps for cds and ribozomal footprint
  cdsRPFOverlaps = countOverlaps(cds,RFP)
  
  #get the five utrs from the coding sequences
  fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
  
  #Find lenth of the five utrs
  fiveUtrLengths = allLengths$utr5_len
  
  #Count overlaps(footprints for ribozomal utr)
  fiveUTRRiboOverlaps = countOverlaps(fiveUTRs,RFP)
  
  #Count overlaps(rna utr)
  fiveUTRRnaOverlaps = countOverlaps(fiveUTRs,rna)
  
  libraryRna = length(rna)
  libraryRPF = length(RFP)
  
  
  correctPositionCdsRna = cdsRnaOverlaps[match(as.character( allLengths$tx_name),names(cdsRnaOverlaps))]
  correctPositionCdsRFP = cdsRPFOverlaps[match(as.character( allLengths$tx_name),names(cdsRPFOverlaps))]
  correctPositionUtrRFP = fiveUTRRiboOverlaps[match(as.character( allLengths$tx_name),names(fiveUTRRiboOverlaps))]
  correctPositionUtrRna = fiveUTRRnaOverlaps[match(as.character( allLengths$tx_name),names(fiveUTRRnaOverlaps))]
}


#Normalization needed for rna
normCDSRNA = FPKMNorlization(unlist(correctPositionCdsRna),allLengths$tx_len,libraryRna)
#Normalization needed for ribozome footprint
normCDSRPF = FPKMNorlization(unlist(correctPositionCdsRFP),unlist(cdsLengths),libraryRPF)


#Nomalization needed for fiveUTRs
normUTRRNA = FPKMNorlization(unlist(correctPositionUtrRna),allLengths$tx_len,libraryRna)
#Normalization needed for ribozome footprint
normUTRRFP = FPKMNorlization(unlist(correctPositionUtrRFP),unlist(fiveUtrLengths),libraryRPF)


#Translational efficiancy cds
teCDS = normCDSRPF/normCDSRNA
#Translational efficiancy utrs
teUTR = normUTRRFP/normUTRRNA

####STEPS AFTER SCAN uoRFS#########

rangesOfuORFs= scanUORFs(fiveUTRs,saveToFile = T)
uorfRNAOverlaps = countOverlaps(unlist(rangesOfuORFs),rna)
uorfRiboOverlaps = countOverlaps(unlist(rangesOfuORFs),RFP)



grl = as.data.frame(rangesOfuORFs)

agUORFRNA = aggregate(uorfRNAOverlaps,by = list(grl$names), sum)
agUORFRibo = aggregate(uorfRiboOverlaps,by = list(grl$names), sum)
#SHIFT

#ORFLengths = sapply(rangesOfuORFs, function(x) sum(width(x)))

ORFLengths = aggregate(grl[,"width"],by=list(grl[,"names"]),sum) 

  
transcriptNames = gsub("_[0-9]*","", agUORFRNA$Group.1)

vector = allLengths$tx_len
names(vector) = allLengths$tx_name
index = match(transcriptNames,names(vector))
vector = vector[index]


normUORFRNA = FPKMNorlization(agUORFRNA$x,vector,libraryRna)
normUORFRFP = FPKMNorlization(agUORFRibo$x,ORFLengths$x,libraryRPF)

teUORF = normUORFRFP/normUORFRNA

#create matrix 9 x genes + normalizations and te's
matrixA = cbind(allLengths,normCDSRNA,normCDSRPF,normUTRRNA,normUTRRFP,teCDS,teUTR)


matrixB = cbind(transcriptNames,agUORFRibo$Group.1, teUORF, normUORFRNA, normUORFRFP)
colnames(matrixB)[2] <- "namesOFuorfs"
matrix = merge(matrixA, matrixB, by.x = "tx_name", by.y = "transcriptNames", all = T)

###Add transcript names


###Index 1 remove all bad values from teCDS
index = matrix[,13] == 0 | is.na(matrix[,13] ) | is.infinite(matrix[,13] ) | is.nan(matrix[,13] )
matrix = matrix[!index,]
###Index 2 remove all bad values from normUORFRNA
index2 = matrix[,"normUORFRNA"] == 0 | is.na(matrix[,"normUORFRNA"] ) | is.infinite(matrix[,"normUORFRNA"] ) | is.nan(matrix[,"normUORFRNA"] )
matrix = matrix[!index2,]

write.csv(matrix, file = "matrix.csv")
plot(log(matrix[,"teCDS"]),log(as.numeric(matrix[,"teUORF"])))

###Find max te uorf per tx


maxUORFTEList = aggregate(as.numeric(as.matrix(matrix$teUORF)),by = list(matrix$tx_name),max)
maxCDSTEList = aggregate(as.numeric(as.matrix(matrix$teCDS)),by = list(matrix$tx_name),max)
###plot of teuorfMax vs teCDS with 0 values for uorf removed

maxTEUORFPLOT = log(maxUORFTEList$x)
maxTEUORFPLOT[which(maxTEUORFPLOT == -Inf)] = NA
maxTECDSPLOT = log(maxCDSTEList$x)
plot(maxTEUORFPLOT,maxTECDSPLOT, main = "teCDS vs maxTEUORF, 0 values for uorf removed")
abline(lm(maxTEUORFPLOT ~ maxTECDSPLOT))

###plot of teuorfmax vs teCDS with 0 values increased to min value
minValueTEUORF = min(maxTEUORFPLOT, na.rm = T)
maxTEUORFPLOT[which(is.na(maxTEUORFPLOT))] = minValueTEUORF
plot(maxTEUORFPLOT,maxTECDSPLOT, main = "teCDS vs maxTEUORF, 0 values increased to min value")
abline(lm(maxTEUORFPLOT ~ maxTECDSPLOT))

###plot normalizations against eachother
qplot(log(matrix$normCDSRNA), y = log(matrix$normCDSRPF))
qplot(log(matrix$normUTRRNA), y = log(matrix$normUTRRFP))
qplot(log(matrix$normUORFRNA), y = log(matrix$normUTRRFP))


#3. creates normalzation of the count
FPKMNorlization = function(counts, lengthSize, librarySize){
  
  result = (as.numeric(counts)*(10^9)) / (as.numeric(lengthSize)*as.numeric(librarySize))
  return(result)
}