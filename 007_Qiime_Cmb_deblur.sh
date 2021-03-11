source activate qiime2-2019.10
wrkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files
qiime deblur denoise-16S \
  --i-demultiplexed-seqs $wrkdir/Combined_CS_joined.qza \
  --p-trim-length 200 \
  --o-representative-sequences $wrkdir/Combined_rep_seqs.qza \
  --o-table $wrkdir/Combined_deblur_table.qza \
  --p-sample-stats \
  --o-stats $wrkdir/Combined_deblur_stats.qza

