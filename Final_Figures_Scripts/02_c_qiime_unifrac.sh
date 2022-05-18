source activate qiime2-2019.10

wkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/FinalFigs_qza

qiime diversity core-metrics-phylogenetic \
  --i-table $wkdir/phy_data__rare_features-table.qza \
  --i-phylogeny $wkdir/phy_data_rare_tree-rooted.qza \
  --p-sampling-depth 334 \
  --m-metadata-file $wkdir/Combined_Metadata.csv \
  --output-dir $wkdir/phy_data_rare-core-metrics-results
