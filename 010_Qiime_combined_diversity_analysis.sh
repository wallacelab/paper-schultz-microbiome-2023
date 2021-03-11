source activate qiime2-2019.10

wkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $wkdir/Combined_rooted_tree.qza \
  --i-table $wkdir/Combined_deblur_table.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Key.csv\
  --output-dir /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/Combined_core_metrics
