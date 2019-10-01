#!/bin/bash

#HT 12/02/19
# script to process and align genomic non paired fastq datasets 
#1 STAR INDICES (If not existing, must be done in seperate script for now)
#2 make directories
#3 Trim adaptor
#4 Remove PhiX
#5 Remove rRNA (from silva)
#6 Remove organism specific ncRNA (from ensembl, zebrafish does not contain tRNA)
#7 Remove organism specific tRNA (from tRNAscan-SE)
#8 Map to organism specific reference

usage(){
cat << EOF
usage: $0 options

script to process and align genomic single end or paired end fastq datasets 

OPTIONS:
	Important options:
	-f	path to input folder (only fasta/q files in folder allowed for now)
		Must be file types of: fasta, fa, fastq, fq or gz (single or paired end reads)
	-o	path to output dir
	-p	paired end? (yes, defualt: no)
	-l	minimum length of reads
	-g	genome dir for all indices (Standard is zebrafish: danrerio10, change to human index if needed etc)
	-s	steps of depletion and alignment wanted:
		(a string: which steps to do? (default: "tr-ge", write "all" to get all: "tr-ph-rR-nc-tR-ge")
			 tr: trim, ph: phix, rR: rrna, nc: ncrna, tR: trna, ge: genome) 
		Write your wanted steps, seperated by "-". Order does not matter.
		To just do trim and alignment to genome write -s "tr-ge"
	-a	adapter sequence for trim (found automaticly if not given), also you can write -a "disable", 
		to disable it
	-t	trim front (default 3) How many bases to pre trim reads 5' end ?
	-A	Alignment type: (default Local, EndToEnd (Local is Local, EndToEnd is force Global))

	Less important options:
	-r	resume?: a character (defualt n) (n for new start fresh with file f from point s, 
			             (if you want a continue from crash specify the step you want to start 
				      from, as defined in -s, start on genome, do -r "ge")
	-m	max cpus allowed (defualt 90)
	-S	include subfolders (defualt n, for no), if you want subfolder do y, for yes.
	-h	this help message

fastp location must be: ~/bin/fastp
STAR location must be: ~/bin/STAR-2.7.0c/source/STAR

NOTE: if STAR is stuck on load, run this line:
STAR --genomeDir /export/valenfs/projects/uORFome/STAR_INDEX_ZEBRAFISH/genomeDir/ --genomeLoad Remove

example usage: RNA_Align_pipeline_folder.sh -f <in.fastq.gz> -o <out_dir>

EOF
}

# Default arguments:
echo $'\nArguments for folder run are the following:'
min_length=15
gen_dir=/export/valenfs/projects/uORFome/STAR_INDEX_ZEBRAFISH
allSteps="tr-ge"
steps=$allSteps
resume="n"
alignment="Local"
adapter="auto"
maxCPU=90
subfolders="n"
trim_front=3
paired="no"
while getopts ":f:o:p:l:g:s:a:t:A:r:m:S:h" opt; do
    case $opt in 
    f)
        in_dir=$OPTARG
        echo "-f input folder $OPTARG"
	;;
    o)
        out_dir=$OPTARG
        echo "-o output folder $OPTARG"
        ;;
    p)
	paired=$OPTARG
        echo "-p output folder $OPTARG"
        ;;
    l)
        min_length=$OPTARG
        echo "-l minimum length of reads $OPTARG"
        ;;
    g)
        gen_dir=$OPTARG
        echo "-g genome dir for all indices $OPTARG"
        ;;
    s)
        steps=$OPTARG
        echo "-s steps to do: $OPTARG"
        ;;
    a)
        adapter=$OPTARG
        echo "-s adapter sequence $OPTARG"
        ;;
    t)
        trim_front=$OPTARG
        echo "-s trim front number $OPTARG"
        ;;
    A)
        alignment=$OPTARG
        echo "-s alignment type $OPTARG"
        ;;
    r)
	resume=$OPTARG
	echo "-r resume (r or new n): $OPTARG"
        ;;
    m)
	maxCPU=$OPTARG
	echo "-m maxCPU: $OPTARG"
        ;;
    S)
	subfolders=$OPTARG
	echo "-S subfolders: $OPTARG"
        ;;
    h)
        usage
        exit
        ;;
    ?) 
        echo "Invalid option: -$OPTARG"
        usage
        exit 1
        ;;
    esac
done

echo $'\n'

if [ $steps == "all" ]; then
	steps="tr-ph-rR-nc-tR-ge"
fi
if [ -z $out_dir ]; then
	echo "Error, out directory (-o) must be speficied!"
	exit 1
fi

#TODO: Fix checks, incase files are not in order
# 1 check if there exists a grouping, and the are equal number in each grouping
# 2 make a list for each end
# send each end into the paired end function	
function findPairs()
{
	f=($1) && shift
	myArray=($@)
	echo "${#myArray[@]}"
	echo $(( ${#myArray[@]} % 2 ))

	if [ $(( ${#myArray[@]} % 2 )) == 1 ]; then
		echo "folder must have even number of fasta/q files, for paired end run!"
		exit 1
	fi
	
	suffix="_001.fastq.gz"
        # remove suffix to get matched pairs
	#for i in "${!myArray[@]}"; do

		#bn=${myArray[i]} | sed  -e "s/$suffix$//"
		#echo ${myArray[i]} | sed  -e "s/$suffix$//"
        	#echo $(basename ${bn}) 

    	#done

	for ((x=0; x<${#myArray[@]}; x = x + 2));
	do
    		echo "running paired end for files: $f/${myArray[x]} and $f/${myArray[x+1]}"
		a="$f/${myArray[x]}"
		b="$f/${myArray[x+1]}"
		/export/valenfs/projects/uORFome/RNA_Align_pipeline.sh -o "$out_dir" -f "$a" -F "$b"  -a "$adapter" -s "$steps" -r "$resume" -l "$min_length" -g "$gen_dir" -m "$maxCPU" -A "$alignment" -t "$trim_front" 
        		
	done			
} 
#TODO: find a way to remove this
function findPairsSub()
{
	myArray=($@)
	echo "${#myArray[@]}"
	echo $(( ${#myArray[@]} % 2 ))

	if [ $(( ${#myArray[@]} % 2 )) == 1 ]; then
		echo "folder must have even number of fasta/q files, for paired end run!"
		exit 1
	fi
	
	for ((x=0; x<${#myArray[@]}; x = x + 2));
	do
    		echo "running paired end with subfolders for files: ${myArray[x]} and ${myArray[x+1]}"
		a="${myArray[x]}"
		b="${myArray[x+1]}"
		/export/valenfs/projects/uORFome/RNA_Align_pipeline.sh -o "$out_dir" -f "$a" -F "$b"  -a "$adapter" -s "$steps" -r "$resume" -l "$min_length" -g "$gen_dir" -m "$maxCPU" -A "$alignment" -t "$trim_front" 
        		
	done			
} 

if [ $paired == "yes" ]; then
	if [ $subfolders == "n" ]; then
		findPairs $in_dir $(ls ${in_dir} | grep '\.fasta\|\.fa\|\.fastq\|\.fq')
	else
		findPairsSub $(find ${in_dir} | grep '\.fasta\|\.fa\|\.fastq\|\.fq\|\.gz' | sort)
	fi
	
else
	for f in $(ls ${in_dir} | grep '\.fasta\|\.fa\|\.fastq\|\.fq')
	do
		echo "running single end for file: $f"
		/export/valenfs/projects/uORFome/RNA_Align_pipeline.sh -o "$out_dir" -f "$f"  -a "$adapter" -s "$steps" -r "$resume" -l "$min_length" -g "$gen_dir" -m "$maxCPU" -A "$alignment" -t "$trim_front" 
	done
fi
