#!/bin/bash

#This gives the absolute paths for each read file, these are needed for the manifest file

wrkdir=//home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Cmb_Trimmed_QualFilt_Sequences
outdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/

##This gets the name of the file. R is recursive and the name of the file is the 9th object in the full string deliminated by /
#for File in $(ls $wrkdir/*R2* | awk -F"/" '{print $9}') ; do
#echo ${File} 
#done

#this puts all the r1 and r2s in their own csv
readlink -f $wrkdir/*R1* > $outdir/R1.csv
readlink -f $wrkdir/*R2* > $outdir/R2.csv

#This gives the sample names
readlink -f $wrkdir/*R1* | awk -F"/" '{print $9}'| awk -F"_" '{print $1}' > $outdir/All_Sample_Names.csv

#This give the fastq file names
readlink -f $wrkdir/*R1* | awk -F"/" '{print $9}' > $outdir/R1_Fastq_Names.csv
readlink -f $wrkdir/*R2* | awk -F"/" '{print $9}' > $outdir/R2_Fastq_Names.csv
