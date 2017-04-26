library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(Biostrings)

#1. create vaiables for rna gtf and bed
if(exists("rna") == F){
  #mRna alligned to the transcriptome
  rna = readGAlignmentPairs("test_data/total_RNA.mock.Rep1.GRCh38.bam")
  #Gene description
  Gtf = makeTxDbFromGFF("test_data/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf","gtf")
  #Ribosome protected footprint
  BED = import.bed("test_data/RP_mock_Rep1.all_runs.bam_shifted_footprints.bed")
  #Fasta file with sequences readDNAStringSet("D:/Hakon/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
  fasta = readDNAStringSet("test_data/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
  
}


#2. get reads

  
  
  #Coding sequences
  cds = cdsBy(Gtf,"tx",use.names = T)
  #
  allLengths = transcriptLengths(Gtf,with.cds_len = T,with.utr5_len = T,with.utr3_len = T)
  cdsLengths = allLengths$cds_len
  
  #Find length of the coding sequences
  
  #Count overlaps between coding sequences and the mRna
  cdsRnaOverlaps = countOverlaps(cds,rna)
  #get the five utrs from the coding sequences
  fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
  #Find lenth of the five utrs
  fiveUtrLengths = lapply(fiveUTRs, function(x) sum(width(x)))
  # Count overlaps for cds and ribozomal footprint
  cdsBEDOverlaps = countOverlaps(cds,BED)
  
  #Count overlaps(footprints for ribozomal utr)
  fiveUTRRiboOverlaps = countOverlaps(fiveUTRs,BED)
  
  #Count overlaps(rna utr)
  fiveUTRRnaverlaps = countOverlaps(fiveUTRs,rna)
  
  libraryRna = length(rna)
  libraryBED = length(BED)
  #Normalization needed for rna
  normRNA = FPKMNorlization(unlist(cdsRnaOverlaps),unlist(cdsLengths),libraryRna)
  #Normalization needed for ribozome footprint
  normBED = FPKMNorlization(unlist(cdsBEDOverlaps),unlist(cdsLengths),libraryBED)
  #Nomalization needed for fiveUTRs
  normUTR = FPKMNorlization(unlist(fiveUTRs),unlist(fiveUtrLengths),libraryRna)
  #Normalization needed for ribozome footprint
  normUTRBED = FPKMNorlization(unlist(fiveUTRs),unlist(fiveUtrLengths),libraryBED)
  
  #Translational efficiancy cds
  teRNA = NormBED/normRNA
  #Translational efficiancy utrs
  teUTR = normUTRBED/normUTR
  
  
  #Correlation 5'utr and cds and plot
  correlation = cor(normUTR,normUTRBED)
  print(correlation)
  plot(normUTR,normUTRBED)
  
  #Correlation te's
  correlation = cor(teRNA,teUTR)
  print(correlation)
  plot(teRNA,teUTR)
  #cor.test(teRNA,teUTR)
  
  ### TODO LIST ###
  #1. statistics:
  # Test if the utrs have a high translational efficiancy (t test, wilcoxen)
  hist(fiveUTRs) #get index
  index = 1
  
  #t-test = t.test(cds[0:index],cds[index+1:end])
  #a. See if they have homogeneous variance
  # if(var.test(fiveUTRs,cds)$p.value < 0.05){
  #   T = t.test(x = fiveUTRs,y = cds,alternative = "two.sided",var.equal = T)
  #   W = wilcox.test(x = fiveUTRs,y = cds,alternative = "two.sided", var.equal = T)
  # }
  #Print all results or save to file
  # print("Translational efficiancy is:")
  # print(te)
  # print("t-test result is:")
  # print(T)
  # print("W-test result is:")
  # print(W)
  # 
  # scanUORFs(fiveUTRs,fasta)
  
  #create matrix 9xgenes
  correctPosition = cdsRnaOverlaps[match(as.character( allLengths$tx_name),names(cdsRnaOverlaps))]
  
  matrix = cbind(allLengths,cdsRnaOverlaps,cds)


#3. creates normalzation of the count
FPKMNorlization = function(counts, lengthSize, librarySize){
  result = (counts*(10^9)) / (lengthSize*librarySize)
  return(result)
}
#4. Find uorfs from sequences in fastafile combined with positions from utrs
scanUORFs = function(fiveUTRs,Fastafile){
  #Read fiveutrs and find all atg that also has a stop codon in frame %3 = 0
  
}

