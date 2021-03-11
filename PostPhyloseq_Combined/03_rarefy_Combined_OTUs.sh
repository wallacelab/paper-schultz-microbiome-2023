source activate qiime2-2019.10

wrkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_filtered_OTUs"

qiime feature-table rarefy \
--i-table $wrkdir/rhizosF-otu-table.qza \
--p-sampling-depth 1133 \
--o-rarefied-table $wrkdir/rhizosF-rarefied-table.qza

qiime feature-table rarefy \
--i-table $wrkdir/rootsF-otu-table.qza \
--p-sampling-depth 2449 \
--o-rarefied-table $wrkdir/rootF-rarefied-table.qza

qiime feature-table rarefy \
--i-table $wrkdir/stalksF-otu-table.qza \
--p-sampling-depth 50 \
--o-rarefied-table $wrkdir/stalksF-rarefied-table.qza
