source activate qiime2-2019.10

clusterDir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_OTU_tables
wkdir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files

qiime tools import \
--input-path $clusterDir/Combined_Rhizosphere_OTU_table.biom \
--type 'FeatureTable[Frequency]' \
--output-path $wkdir/Combined_Rhizosphere_OTU_table.qza


qiime tools import \
--input-path $clusterDir/Combined_Stalks_OTU_table.biom \
--type 'FeatureTable[Frequency]' \
--output-path $wkdir/Combined_Stalks_OTU_table.qza

qiime tools import \
--input-path $clusterDir/Combined_Root_OTU_table.biom \
--type 'FeatureTable[Frequency]' \
--output-path $wkdir/Combined_Root_OTU_table.qza

#website says you need an argument "--source-format BIOMV210Format \", it does not work if
#you include that argument
