source activate qiime2-2019.10

clusterDir=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files


qiime tools export --input-path $clusterDir/Combined_deblur_table.qza --output-path $clusterDir/
mv $clusterDir/feature-table.biom $clusterDir/Combined_table.biom

biom convert -i $clusterDir/Combined_table.biom -o $clusterDir/Combined_table.tsv --to-tsv
sed "s/#//g" $clusterDir/Combined_table.tsv > $clusterDir/Combined_table.forR.tsv


#this creates a readable OTU table, where we can see how many otus and what they are in each sample. 
