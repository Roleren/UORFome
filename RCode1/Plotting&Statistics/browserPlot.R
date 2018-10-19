setwd('/export/valenfs/projects/uORFome/')
# Gviz example
library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
library(gridExtra)

# bioconductor
library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
library(Biostrings)
library(Rsamtools)

transluORFs<-  GRanges(c('chr22','chr22','chr22','chr22'),
                       IRanges(c(39520586,39520648, 39520747, 39521354),width = c(6,12,5,175)))
names(transluORFs)<- c('uORF_1','uORF_2','uORF_3','uORF_3')
strand(transluORFs)<- '+'
mcols(transluORFs)$feature<- c('not-predicted','predicted','predicted','predicted')
transluORFTrack<- GeneRegionTrack(transluORFs,
                                  showId=TRUE,
                                  id=names(transluORFs),
                                  name='uORF predictions',
                                  transcript=names(transluORFs),
                                  symbol=names(transluORFs),
                                  feature = transluORFs$feature)


# define region of interest:
mychrom='chr22'
mystart=min(start(transluORFs))
myend=max(end(transluORFs))
myregion_extend=GRanges(mychrom,IRanges(mystart-1000,myend+1000))

options(ucscChromosomeNames=FALSE)
genomeTrack<- GenomeAxisTrack()

orfFile='/export/valenfs/projects/uORFome/dataBase/bedUniqueUorfs.bed'
orfs<- import.bed(orfFile)
orfs<- unique(orfs)
orfTrack<- AnnotationTrack(orfs[orfs %over% myregion_extend], name="uORFs based on sequence", col='darkgrey', fill='darkgrey')



gtfFile='/export/valenfs/data/references/GRCh38_human/Homo_sapiens.GRCh38.79.chr.gtf'
txdb<- makeTxDbFromGFF(gtfFile) 
geneTrack<- GeneRegionTrack(txdb,name="gene models",showId=TRUE,geneSymbol=TRUE) # can get a txdb object

# get the sequence in region
fastaFile='/export/valenfs/data/references/GRCh38_human/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa'
fa=FaFile(fastaFile)
refseq<- readDNAStringSet(fastaFile)
#refseq<- getSeq(fa, myregion_extend)
names(refseq)<- gsub(" .*","",names(refseq))
seqTrack<- SequenceTrack(refseq, chromosome = mychrom)


# read in ribo-seq wig files here, either shifted all, or just a few
riboWigFileF='/export/valenfs/data/processed_data/Ribo-seq/fantom_human_tracks/tracks/Stern-Ginossar_N_2012.Human.fibroblasts.RPF.GRCh38.SRR592954-forward.wig'
riboWigFileR='/export/valenfs/data/processed_data/Ribo-seq/fantom_human_tracks/tracks/Stern-Ginossar_N_2012.Human.fibroblasts.RPF.GRCh38.SRR592954-reverse.wig'
wigF<- import.wig(riboWigFileF)
riboTrackF<- DataTrack(chromosome = seqnames(wigF),
                       start=start(wigF),
                       end=end(wigF),
                       data=score(wigF),
                       genome='',
                       name='Ribo-seq Stern-Ginossar 2012 \nSRR592954 (fibroblast) P-site aligned')
wigR<- import.wig(riboWigFileR)
riboTrackR<- DataTrack(chromosome = seqnames(wigR),
                       start=start(wigR),
                       end=end(wigR),
                       data=score(wigR),
                       genome='',
                       name='Ribo-seq Reverse')


# make a cage track # find the right fibroblast cell line
#cageFile=''
displayPars(riboTrackF)<- list(showAxis=TRUE,
                               type='h',
                                   col='darkblue',
                                   fontcolor.title='black',
                                   col.axis='black',
                                   cex.axis=0.5,
                                   background.title='white')

displayPars(geneTrack)<- list(
                               fontcolor.title='black',
                               col.axis='black',
                               cex.axis=0.5,
                               background.title='white')

displayPars(transluORFTrack)<- list(
  fontcolor.title='black',
  col.axis='black',
  fill='darkred',
  cex.group=0.7,
  cex.axis=0.5,
  background.title='white')

displayPars(orfTrack)<- list(
  fontcolor.title='black',
  col.axis='black',
  fill='darkgrey',
  cex.axis=0.5,
  background.title='white')

displayPars(seqTrack)<- list(fontcolor=biovizBase::getBioColor(),noLetters=FALSE,fontsize=6)


# use this to highlight
mytracks<- list(genomeTrack,
                #cageTrack,
                riboTrackF,
                #riboTrackR,
                transluORFTrack,
                #orfTrack,
                #seqTrack,
                geneTrack
)

ht <- HighlightTrack(trackList = mytracks,  start = start(transluORFs['uORF_2'])-80,
                     end = end(transluORFs['uORF_2'])+20,
                     chromosome = seqnames(transluORFs['uORF_2']))
pdf('browser_snapshot.pdf', width=10, height=6)
plotTracks(ht,
           from = mystart - 1300,
           to = myend + 500,
           chromosome = mychrom)

dev.off()


# use this to highlight
mytracks<- list(genomeTrack,
                #cageTrack,
                riboTrackF,
                #riboTrackR,
                transluORFTrack,
                #orfTrack,
                seqTrack,
                orfTrack
                #geneTrack
)

pdf('browser_close.pdf',width=10,height=6)
plotTracks(mytracks,
           from = start(transluORFs['uORF_2'])-85,
           to = end(transluORFs['uORF_2'])+10,
           chromosome = seqnames(transluORFs['uORF_2']))

dev.off()


