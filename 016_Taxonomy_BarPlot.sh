source activate qiime2-2019.10

taxdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/TaxonomyStuff

qzadir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files

qiime taxa barplot \
  --i-table $qzadir/Combined_deblur_table.qza \
  --i-taxonomy $taxdir/Combined.taxonomy.qza  \
  --m-metadata-file /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Key.csv \
  --o-visualization $taxdir/Combined-taxa-bar-plots.qzv
