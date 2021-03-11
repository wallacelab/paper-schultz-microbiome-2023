source activate qiime2-2019.10

taxdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/TaxonomyStuff

qzadir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files

qiime feature-classifier classify-sklearn \
  --i-classifier $taxdir/silva-132-99-nb-classifier.qza \
  --i-reads $qzadir/Combined_rep_seqs.qza \
  --o-classification $taxdir/Combined.taxonomy.qza

qiime metadata tabulate \
  --m-input-file $taxdir/Combined.taxonomy.qza \
  --o-visualization $taxdir/Combined.taxonomy.qzv
