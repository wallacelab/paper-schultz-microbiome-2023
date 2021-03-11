#!/bin/bash

mkdir /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/GH_CS/FastQC_Trimmed_All

fastqc=/home/coreyschultz/1.Programs.Installs/fastqc_v0.11.8/FastQC/fastqc

wrkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Cmb_Trimmed_QualFilt_Sequences

for File in $(ls $wrkdir/*R*);do

$fastqc --extract $File --outdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/GH_CS/FastQC_Trimmed_All

done


#cd into this directory and run multiqc . to generate a report for EVERY sample
