source activate qiime2-2019.10
wrkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS
qzaDir=Combined_qza_files
qiime feature-table summarize \
  --i-table $wrkdir/$qzaDir/Combined_deblur_table.qza \
  --o-visualization $wrkdir/$qzaDir/Combined_deblur_table.qzv \
  --m-sample-metadata-file $wrkdir/Combined_Key.csv
qiime feature-table tabulate-seqs \
  --i-data $wrkdir/$qzaDir/Combined_rep_seqs.qza \
  --o-visualization $wrkdir/$qzaDir/Combined_rep_seqs.qzv

#This creates the sequences we will use to determine what is in our samples and how the samples segregate. 
