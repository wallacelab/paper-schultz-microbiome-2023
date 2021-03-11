source activate qiime2-2019.10

wrkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied"

qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_rhizosphere-core-metrics-rarefied-results/faith_pd_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_rhizosphere-core-metrics-rarefied-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_rhizosphere-core-metrics-rarefied-results/evenness_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_rhizosphere-core-metrics-rarefied-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_rhizosphere-core-metrics-rarefied-results/shannon_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_rhizosphere-core-metrics-rarefied-results/shannon-group-significance.qzv


qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_roots-core-metrics-rarefied-results/faith_pd_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_roots-core-metrics-rarefied-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_roots-core-metrics-rarefied-results/evenness_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_roots-core-metrics-rarefied-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_roots-core-metrics-rarefied-results/shannon_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_roots-core-metrics-rarefied-results/shannon-group-significance.qzv


qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_stalks-core-metrics-rarefied-results/faith_pd_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_stalks-core-metrics-rarefied-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_stalks-core-metrics-rarefied-results/evenness_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_stalks-core-metrics-rarefied-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $wrkdir/Cmb_stalks-core-metrics-rarefied-results/shannon_vector.qza \
  --m-metadata-file $wrkdir/Combined_Key.csv \
  --o-visualization $wrkdir/Cmb_stalks-core-metrics-rarefied-results/shannon-group-significance.qzv
