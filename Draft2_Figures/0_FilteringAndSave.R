# Generate Phyloseq Objects 
library(ggplot2)
library(phyloseq)
library(plotrix)
library(vegan)
library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(tidyr)
library(microbiome)
library(ape)

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS")

# Read in the data

metadata <- read.csv("Combined_Key.csv", header = TRUE, sep = "\t")
metadata$Group <- paste(metadata$Sample_Type_Blanks_differentiated, metadata$Experiment, sep = "_")

SVs <- read_qza("Combined_qza_files/Combined_deblur_table.qza")

taxonomy <-read_qza("Combined_Results/TaxonomyStuff/Combined.taxonomy.qza")
taxtable <-taxonomy$data %>% as.data.frame() 
taxtable <- taxtable %>% separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") #convert the table into a tabular split version

tree <- read_qza("Combined_qza_files/Combined_rooted_tree.qza")

phy <- qza_to_phyloseq("Combined_qza_files/Combined_deblur_table.qza", "Combined_qza_files/Combined_rooted_tree.qza", "Combined_Results/TaxonomyStuff/Combined.taxonomy.qza", "Combined_Metadata.csv")
phy

##### FILTERING

##    Taxonomic Filtering: Supervised
# creat taxonomy table: number of features for each phyla
table(tax_table(phy)[, "Phylum"], exclude = NULL)

# Remove phyla for which only one feature was observed. And ambiguous calls


# compute prevalence of each feature, store as data frame
prevdf = apply(X = otu_table(phy),
               MARGIN = ifelse(taxa_are_rows(phy), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phy),
                    tax_table(phy))
#Computes the total and average prevalences of the features in each phylum
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Filter out Chlorophlasts and Mitochondria
phy
psm = subset_taxa(phy,Family!="Mitochondria")
psm
psc = subset_taxa(psm,Order!="Chloroplast")
psc

#     Prevalence Filtering: Unsupervised

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(psc, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / phyloseq::nsamples(phy),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

prevalenceThreshold = 0.05 * phyloseq::nsamples(psc)
prevalenceThreshold

keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phy2 = prune_taxa(keepTaxa, psc)
phy2

#Filter out all blank OTU's           #When should I do this??
blanks = subset_samples(phy2,Sample_Type_Blanks_differentiated=="Blank-Stalk" |
                          Sample_Type_Blanks_differentiated=="Blank-Rhizosphere" |
                          Sample_Type_Blanks_differentiated=="Blank-Root" |
                          Sample_Type_Blanks_differentiated=="Blank-Soil") 

######This works - our main data
phy_noblank <- prune_samples(!(sample_names(phy2) %in% sample_names(blanks)),phy2)
phy_data <- prune_taxa(taxa_sums(blanks) < 1, phy_noblank)


### Shared OTUs
# Split by experiment, so you can drop OTUs that are not found in all three experiments

GH <- subset_samples(phy_data, Experiment=="GH")
END <- subset_samples(phy_data, Experiment=="END")
MMH <- subset_samples(phy_data, Experiment=="MMH")


#If its not in all three experiments, drop it. We have to do it step by step to not get errors. 
phy3 <- prune_taxa(taxa_sums(GH) > 0, phy_data)
MMH2 <- prune_taxa(taxa_sums(GH) > 0, MMH)
phy4 <- prune_taxa(taxa_sums(MMH2) > 0, phy3)
END2 <- prune_taxa(taxa_sums(GH) > 0, END)
END3 <- prune_taxa(taxa_sums(MMH2) > 0, END2)

phy_shared <- prune_taxa(taxa_sums(END3) > 0, phy4)


# # Class, Genus, or ASV level
# phyClass <- tax_glom(phyCmbFilt, taxrank = "Class")
# phyGenus <- tax_glom(phyCmbFilt, taxrank = "Genus")
# phyHigh <- phyCmbFilt
# 
# # Use ASV level for Volcano Plot
# phyHigh
# 
# #Before filtering out by experiments how many taxa do I have?
# phyClassPreFilt <- tax_glom(phy_data, taxrank = "Class")

#Tissue
###############
stalks <- subset_samples(phy_data, Sample_Type_Blanks_differentiated=="Stalk")
rhizos <- subset_samples(phy_data, Sample_Type_Blanks_differentiated=="Rhizosphere")
roots <- subset_samples(phy_data, Sample_Type_Blanks_differentiated=="Root")

###### Datasets
phy2 #Raw
phy_data #Main Data
phy_shared #Only in all three experiments 
######

# Filter and rareify in other scripts 
# phy_data_rare <- rarefy_even_depth(phy_data, sample.size = min(sample_sums(phy_data)),
#                                    replace = TRUE)


setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/Phyloseq_Objects")
# export than test import - otu table, taxonomy, tree, metadata

#Export phy2
write.csv(otu_table(phy2),"phy2_otu.csv")
write.csv(tax_table(phy2),"phy2_taxonomy.csv")
write.csv(sample_data(phy2),"phy2_sampdat.csv")
ape::write.tree(phy_tree(phy2), "phy2_tree.tree")

# import test
otu = otu_table(as.matrix(read.table("phy2_otu.csv", row.names = 1, header = TRUE, sep = ",")),taxa_are_rows = TRUE)
taxa = tax_table(as.matrix(read.table("phy2_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")))
meta = sample_data(data.frame(read.table("phy2_sampdat.csv", row.names = 1, header = TRUE, sep = ",")))
tree = read.tree("phy2_tree.tree")

phy2_test = phyloseq(otu,taxa)
phy2_test = merge_phyloseq(phy2_test, meta, tree)
phy2_test


#Export phy_data
write.csv(otu_table(phy_data),"phy_data_otu.csv")
write.csv(tax_table(phy_data),"phy_data_taxonomy.csv")
write.csv(sample_data(phy_data),"phy_data_sampdat.csv")
ape::write.tree(phy_tree(phy_data), "phy_data_tree.tree")

#Export phy_shared
write.csv(otu_table(phy_shared),"phy_shared_otu.csv")
write.csv(tax_table(phy_shared),"phy_shared_taxonomy.csv")
write.csv(sample_data(phy_shared),"phy_shared_sampdat.csv")
ape::write.tree(phy_tree(phy_shared), "phy_shared_tree.tree")

