source activate qiime2-2019.10
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Manifest.tsv \
  --output-path /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_CS.qza\
  --input-format PairedEndFastqManifestPhred33V2 
