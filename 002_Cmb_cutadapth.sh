#!/bin/bash

wrkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/All_Sequences
outputdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Cmb_Trimmed_QualFilt_Sequences


#The awk command says take this and break it into objects with the -F "." delimiter given after the -F
#then the '{print $1,$2,$3}' part says take the first, second and third objects of those seperated by the #previous delimter and it will print with a space between them, then sed replaces that space with _
#so this line is saying: for vvariable file in list all the files in the work dir,
#deliminate file name by _ and  gives us the first three parts, which are everything before _R1/2. But they are separated and need to be concatinated.
#sed: | pipes everything into sed. s stands for substitute, so substitute a space, specifically the spaces inbetween objects 1 2 and 3, and merges them with an underscore _ globally (g) 


for File in $(ls $wrkdir/ | awk -F"_" '{print $1,$2,$3}' | sed "s/ /_/g") ; do
#echo ${File} 


#echo and then put " around the command will show you what the computer sees

#-m 1 says it needs at least 1 base to look at, otherwise it includes blanks that qiime doesnt like
#-q 26 says that the phred quality score needs to be at least 26
#-u -20 removes the last 20 bases Matt had this in his comments but not his actual script? Might have to make it 100 based on jasons slack

cutadapt -m 1 -q 26 -g NNNNNNNNGAGTGYCAGCMGCCGCGGTAA -G NNNNNNNNGAGGACTACNVGGGTWTCTAAT -o $outputdir/${File}_R1_001.trimmed.fastq.gz -p $outputdir/${File}_R2_001.trimmed.fastq.gz $wrkdir/${File}_R1_001.fastq.gz $wrkdir/${File}_R2_001.fastq.gz
done 


