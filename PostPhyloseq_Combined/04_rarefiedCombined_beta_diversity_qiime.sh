source activate qiime2-2019.10

wrkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_filtered_OTUs"
resdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $wrkdir/rhizosF_tree.qza \
  --i-table $wrkdir/rhizosF-rarefied-table.qza \
  --p-sampling-depth 1133 \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --output-dir $resdir/Cmb_rhizosphere-core-metrics-rarefied-results \
  
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $wrkdir/rootsF_tree.qza \
  --i-table $wrkdir/rootF-rarefied-table.qza \
  --p-sampling-depth 2449 \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --output-dir $resdir/Cmb_roots-core-metrics-rarefied-results \

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $wrkdir/stalksF_tree.qza \
  --i-table $wrkdir/stalksF-rarefied-table.qza \
  --p-sampling-depth 50 \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --output-dir $resdir/Cmb_stalks-core-metrics-rarefied-results \
