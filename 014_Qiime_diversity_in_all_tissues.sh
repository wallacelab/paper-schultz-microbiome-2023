source activate qiime2-2019.10

wkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files
key_file=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Key.csv
analysis_results=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $wkdir/Combined_rooted_tree.qza \
  --i-table $wkdir/Combined_Rhizosphere_OTU_table.qza \
  --p-sampling-depth 1133 \
  --m-metadata-file $key_file \
  --output-dir $analysis_results/Combined_Rhizosphere_Metrics_PreFilt2

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $wkdir/Combined_rooted_tree.qza \
  --i-table $wkdir/Combined_Root_OTU_table.qza \
  --p-sampling-depth 2449 \
  --m-metadata-file $key_file \
  --output-dir $analysis_results/Combined_Root_Metrics_PreFilt2

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $wkdir/Combined_rooted_tree.qza \
  --i-table $wkdir/Combined_Stalks_OTU_table.qza \
  --p-sampling-depth 50 \
  --m-metadata-file $key_file \
  --output-dir $analysis_results/Combined_Stalk_Metrics_PreFilt2
