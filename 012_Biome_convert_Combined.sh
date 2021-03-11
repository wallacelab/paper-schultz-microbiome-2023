source activate qiime2-2019.10

#This takes tsv files and converts them to biome
clusterDir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_OTU_tables


biom convert -i $clusterDir/Combined_Rhizosphere_OTU_table.tsv -o $clusterDir/Combined_Rhizosphere_OTU_table.biom --to-hdf5 --table-type="OTU table" 


biom convert -i $clusterDir/Combined_Stalks_OTU_table.tsv -o $clusterDir/Combined_Stalks_OTU_table.biom --to-hdf5 --table-type="OTU table" 

biom convert -i $clusterDir/Combined_Root_OTU_table.tsv -o $clusterDir/Combined_Root_OTU_table.biom --to-hdf5 --table-type="OTU table" 
