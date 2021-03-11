# This script is ran after you make the tsv file from the deblur table. This subsets OTU table by sample type

OTU_table = read.table("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_table.forR.tsv", sep = "\t", skip = 1)
key_file = read.csv("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Key.csv", sep = "\t")

#This gets the otu table and key file in the right format to merge with 
OTU_Transposed=t(OTU_table)
colnames(OTU_Transposed)=OTU_Transposed[1,]
OTU_proper_header=OTU_Transposed[-1,]

#Take the columns we need from key file
Important_cols=key_file[c("SampleID", "Sample_Type")]

names(Important_cols)=c("OTU_ID", "Sample_Type")

Merged_file=OTU_Key_merge=merge(Important_cols,OTU_proper_header, by.x = "OTU_ID", by.y = "OTU ID")

#THis subsets the otu information by sample type
OTU_Stalk=subset(Merged_file,Sample_Type=="Stalk")
OTU_Rhizosphere=subset(Merged_file,Sample_Type=="Rhizosphere")
OTU_Root=subset(Merged_file,Sample_Type=="Root")

#Drop the sample type column
OTU_Stalk_no_sample_type=OTU_Stalk[,-2]
OTU_Rhizospher_no_sample_type=OTU_Rhizosphere[,-2]
OTU_Root_no_sample_type=OTU_Root[,-2]

#Formating for writing to file
Transposed_stalks=t(OTU_Stalk_no_sample_type)
Transposed_stalks_DF=as.data.frame(Transposed_stalks)
write.table(Transposed_stalks_DF, file = '/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_OTU_tables/Combined_Stalks_OTU_table.tsv', quote = F, col.names = F, sep = "\t")

Transposed_RZ=t(OTU_Rhizospher_no_sample_type)
Transposed_RZ_DF=as.data.frame(Transposed_RZ)
write.table(Transposed_RZ_DF,file = "//home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_OTU_tables/Combined_Rhizosphere_OTU_table.tsv", quote = F, col.names = F, sep = "\t")

Transposed_Root=t(OTU_Root_no_sample_type)
Transposed_Root_DF=as.data.frame(Transposed_Root)
write.table(Transposed_Root_DF,file = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_OTU_tables/Combined_Root_OTU_table.tsv", quote = F, col.names = F, sep = "\t")



