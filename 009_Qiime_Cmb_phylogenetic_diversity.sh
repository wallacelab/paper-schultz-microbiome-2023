source activate qiime2-2019.10
wkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $wkdir/Combined_rep_seqs.qza \
  --o-alignment $wkdir/Combined_aligned_rep_seqs.qza \
  --o-masked-alignment $wkdir/Combined_masked_aligned_rep_seqs.qza \
  --o-tree $wkdir/Combined_unrooted_tree.qza \
  --o-rooted-tree $wkdir/Combined_rooted_tree.qza

