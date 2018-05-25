This document describes the uORF data-base and how to use it.




# Introduction
The uORF data-base is a sql-lite data-base with tables from result of human uorf predictions with tissue variance.
It consists of uORF ids, features and many more tables. 

It is a pipeline wrapper around our bioconductor package ORFik.
https://bioconductor.org/packages/release/bioc/html/ORFik.html

# how to use:

1. First download UORFome to your R library:

devtools::install_github("UORFome")
Note: If you want our sqllite data-base with uORFs (25GB) you need to send an email.

To use the database put it relative to the RCODE1 folder as following.
Create a folder called database in same folder as RCODE1,
put the sqlite database there.

Path relative to RCODE1 is then:
/RCODE1/../database/database.sqlite

2. To start in R source the files with a uORF database included:

setwd("/export/valenfs/projects/uORFome/RCode1/") # your own location of the DataBaseCreator.R script
source("./DataBaseSetup.R") # here you will get error if you do not have all packages needed, then install them

# to see available tables in the data-base:

listTables()

You will now see the tables available, something like this:

[1] "ORFScores"                  "RRS"                       
 [3] "RSS"                        "Ribofpkm"                  
 [5] "SplittedByExonsuniqueUORFs" "cageInformation"           
 [7] "disengagementScores"        "distORFCDS"                
 [9] "entropyRFP"                 "floss"                     
[11] "fractionLengths"            "inFrameCDS"                
[13] "ioScore"                    "isOverlappingCds"          
[15] "kozak"                      "linkORFsToTx"              
[17] "rankInTx"                   "riboAll"                   
[19] "tissueAtlasByCage"          "uniqueIDs"                 
[21] "uorfsAsGRWithTx"     

# Extract ORFs as GRanges

If we want the uorfs, as GRangesList, each group is a uORF:

grl <- getUORFsInDb()  

Now lets pick out the orfs and extract them as GRangesList, each group is uORFs per transcript:
grl <- getUORFsInDb()  
g <- unlistGrl(grl)
grl <- groupGRangesBy(g)

Now you have the uORFs by transcript

You see the grl object have a column called names
This column is the ORF identifier in the transcript, so if two exons in a transcript group 
have the same name, it means they are the exons of one uORF
i.g. 
_1 
_1 # <- exon two of uORF one
_2
_3 # <- only exon of uORF three in that transcript


The primary key of the data-base is the uorf ID. 
It is defined as a string of this syntax:

chromosome, strand, start width

i.g:
"chr1, +, 18124 154"

to extract unique uorfIDs do:

uorfIDs <- readTable("uniqueIDs")  

In the database there is a table called linkORFsToTx

linkUoRFTx <- readTable("linkORFsToTx")

This table gives 2 columns, 1st is uORF id, 2nd is transcript that uORF was found on. 

# Extracting predictied uORFs

The logical table uorfPredictions contains the predicted uORFs.
To get them as GRangesList, do:

grl <- getUORFsInDb()  
pred <- readTable("uorfPredictions")

predictedUORFs <- grl[pred]

3. Making a new database:

To make a new data-base, copy the RCode1 folder to your prefered location.
In the script HelperVariables.R you need to set all the General paths, to these data:
Cage-data (bed files)
RiboSeq-data (as bed files shifted to p-site)
RNASeq-data (bam files)

Then run the script
