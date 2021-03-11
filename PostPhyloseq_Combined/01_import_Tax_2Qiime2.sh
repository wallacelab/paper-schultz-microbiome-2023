source activate qiime2-2019.10

wrkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_filtered_OTUs"

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $wrkdir/rhizosF_tax.txt \
--output-path $wrkdir/rhizosF_taxonomy.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $wrkdir/rootsF_tax.txt \
--output-path $wrkdir/rootsF_taxonomy.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $wrkdir/stalksF_tax.txt \
--output-path $wrkdir/stalksF_taxonomy.qza
