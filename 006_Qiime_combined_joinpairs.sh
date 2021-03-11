source activate qiime2-2019.10
qiime vsearch join-pairs \
  --i-demultiplexed-seqs /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_CS.qza \
  --o-joined-sequences /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_CS_joined.qza \
