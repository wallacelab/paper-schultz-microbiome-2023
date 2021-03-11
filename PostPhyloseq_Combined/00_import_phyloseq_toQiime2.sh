source activate qiime2-2019.10

wrkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_filtered_OTUs"

qiime tools import \
--input-path $wrkdir/rhizosF_features-table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path $wrkdir/rhizosF-otu-table.qza

qiime tools import \
--input-path $wrkdir/rootsF_features-table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path $wrkdir/rootsF-otu-table.qza

qiime tools import \
--input-path $wrkdir/stalksF_features-table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path $wrkdir/stalksF-otu-table.qza

