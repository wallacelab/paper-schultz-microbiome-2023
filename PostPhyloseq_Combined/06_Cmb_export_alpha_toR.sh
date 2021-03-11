source activate qiime2-2019.10

inputDir1=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied/Cmb_rhizosphere-core-metrics-rarefied-results

qiime tools export --input-path $inputDir1/faith_pd_vector.qza --output-path $inputDir1/Cmb_rhizosphere_faithPD


inputDir2=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied/Cmb_roots-core-metrics-rarefied-results

qiime tools export --input-path $inputDir2/faith_pd_vector.qza --output-path $inputDir2/Cmb_roots_faithPD



inputDir3=/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied/Cmb_stalks-core-metrics-rarefied-results

qiime tools export --input-path $inputDir3/faith_pd_vector.qza --output-path $inputDir3/Cmb_stalks_faithPD
