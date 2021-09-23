source activate qiime2-2019.10

clusterDir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/FinalFigs_qza
wkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/FinalFigs_qza

qiime tools import \
--input-path $clusterDir/phy2_features-table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path $wkdir/phy2_features-table.qza


qiime tools import \
--input-path $clusterDir/phy2_tree-rooted.newick \
--type 'Phylogeny[Rooted]' \
--output-path $wkdir/phy2_tree-rooted.qza



#website says you need an argument "--source-format BIOMV210Format \", it does not work if
#you include that argument
