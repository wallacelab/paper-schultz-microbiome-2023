# Script to identify % Similarity of microbiomes

library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(microbiome)
library(microbiomeutilities)
library(ggVennDiagram)

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS")

# Read in the data

metadata <- read.csv("Group_Metadata.csv", header = TRUE, sep = ",")

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

# Split by experiment, so you can drop OTUs that are not found in all three experiments

GH <- subset_samples(phy2, Experiment=="GH")
END <- subset_samples(phy2, Experiment=="END")
MMH <- subset_samples(phy2, Experiment=="MMH")

get_sample(GH, "b8a52db1ff47d14bd94a7c36dbe4f92d")

#Combined Data
phy2

Roots_phy <- subset_samples(phy2, Sample_Type=="Root")
Rhizos_phy <- subset_samples(phy2, Sample_Type=="Rhizosphere")
Stalks_phy <- subset_samples(phy2, Sample_Type=="Stalk")

### By Inbred/Hybrid

table(meta(phy2)$Inbred_or_Hybrid)

# All tissues
BackgroundGroups <- c("Inbred", "Open_Pollinated", "Hybrid")      # don't want soil or blanks
print(BackgroundGroups)
list_core <- c() # empty object
for (n in BackgroundGroups){
  print(as.name(n))
  core.sub <- subset_samples(phy2, BackgroundGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .00001,            # .001 in atleast 90% of samples
                         prevalence = 0.0001)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing Background OTUs: All Tissues")


# Roots
BackgroundGroups <- c("Inbred", "Open_Pollinated", "Hybrid")      # don't want soil or blanks
print(BackgroundGroups)
list_core <- c() # empty object
for (n in BackgroundGroups){
  print(as.name(n))
  core.sub <- subset_samples(Roots_phy, BackgroundGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .0001,            # .001 in atleast 90% of samples
                         prevalence = 0.0001)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing Background OTUs: Roots")

# Rhizos
BackgroundGroups <- c("Inbred", "Open_Pollinated", "Hybrid")      # don't want soil or blanks
print(BackgroundGroups)
list_core <- c() # empty object
for (n in BackgroundGroups){
  print(as.name(n))
  core.sub <- subset_samples(Rhizos_phy, BackgroundGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .0001,            # .001 in atleast 90% of samples
                         prevalence = 0.0001)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing Background OTUs: Rhizos")

# Stalks
BackgroundGroups <- c("Inbred", "Open_Pollinated", "Hybrid")      # don't want soil or blanks
print(BackgroundGroups)
list_core <- c() # empty object
for (n in BackgroundGroups){
  print(as.name(n))
  core.sub <- subset_samples(Stalks_phy, BackgroundGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .0001,            # .001 in atleast 90% of samples
                         prevalence = 0.0001)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing Background OTUs: Stalks")

