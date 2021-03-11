source activate qiime2-2019.10

wrkdir="/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_filtered_OTUs"

qiime tools import \
--type 'Phylogeny[Rooted]' \
--input-path $wrkdir/rhizosF_tree-rooted.newick \
--output-path $wrkdir/rhizosF_tree.qza

qiime tools import \
--type 'Phylogeny[Rooted]' \
--input-path $wrkdir/rootsF_tree-rooted.newick \
--output-path $wrkdir/rootsF_tree.qza

qiime tools import \
--type 'Phylogeny[Rooted]' \
--input-path $wrkdir/stalksF_tree-rooted.newick \
--output-path $wrkdir/stalksF_tree.qza
