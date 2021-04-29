# Take pre proccessed qiime objects, filter with phyloseq, then insert back into qiime for diversity metrics!
# Why did I create a second version of this?

library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS")

# Read in the data

metadata <- read.csv("Combined_Key.csv", header = TRUE, sep = "\t")

SVs <- read_qza("Combined_qza_files/Combined_deblur_table.qza")

taxonomy <-read_qza("Combined_Results/TaxonomyStuff/Combined.taxonomy.qza")
taxtable <-taxonomy$data %>% as.data.frame() 
taxtable <- taxtable %>% separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") #convert the table into a tabular split version

tree <- read_qza("Combined_qza_files/Combined_rooted_tree.qza")

phy <- qza_to_phyloseq("Combined_qza_files/Combined_deblur_table.qza", "Combined_qza_files/Combined_rooted_tree.qza", "Combined_Results/TaxonomyStuff/Combined.taxonomy.qza", "Combined_Key.csv")
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
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(phy),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

prevalenceThreshold = 0.05 * nsamples(psc)
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
#If its not in all three experiments, drop it. We have to do it step by step to not get errors. 
phy3 <- prune_taxa(taxa_sums(GH) > 0, phy2)
MMH2 <- prune_taxa(taxa_sums(GH) > 0, MMH)
phy4 <- prune_taxa(taxa_sums(MMH2) > 0, phy3)
END2 <- prune_taxa(taxa_sums(GH) > 0, END)
END3 <- prune_taxa(taxa_sums(MMH2) > 0, END2)
phyCmbFilt <- prune_taxa(taxa_sums(END3) > 0, phy4)

#Split by samples then remove blank OTUs and do phyloseq to qiime

stalks <- subset_samples(phyCmbFilt, Sample_Type_Blanks_differentiated=="Stalk")
rhizos <- subset_samples(phyCmbFilt, Sample_Type_Blanks_differentiated=="Rhizosphere")
roots <- subset_samples(phyCmbFilt, Sample_Type_Blanks_differentiated=="Root")

#Filter out all blank OTU's           #When should I do this??
blank_stalks = subset_samples(phy2,Sample_Type_Blanks_differentiated=="Blank-Stalk")
blank_rhizos = subset_samples(phy2,Sample_Type_Blanks_differentiated=="Blank-Rhizosphere")
blank_roots = subset_samples(phy2,Sample_Type_Blanks_differentiated=="Blank-Root")


#This is all taxa that have a non zero value in the blanks
otu_table(prune_taxa(taxa_sums(blank_stalks) > 1, blank_stalks))

sbgood <-  prune_taxa(taxa_sums(blank_stalks) < 1, blank_stalks)
sbad <- prune_taxa(taxa_sums(blank_stalks) > 1, blank_stalks)
goodlist <- taxa_names(sbgood)
badlist <- taxa_names(sbad)

#This works 
stalksF <- prune_taxa(goodlist, stalks)

Rhizegood <-  prune_taxa(taxa_sums(blank_rhizos) < 1, blank_rhizos)
Rgoodlist <- taxa_names(Rhizegood)
rhizosF <- prune_taxa(Rgoodlist, rhizos)

Rootgood <-  prune_taxa(taxa_sums(blank_roots) < 1, blank_roots)
Rootgoodlist <- taxa_names(Rootgood)
rootsF <- prune_taxa(Rootgoodlist, roots)

#data sets - make compositional?
stalksF
rhizosF
rootsF

stalksR <- microbiome::transform(stalksF, "compositional")
rhizosR <- microbiome::transform(rhizosF, "compositional")
phyCmbComp <- microbiome::transform(phyCmbFilt, "compositional")
phy2com <- microbiome::transform(phy2, "compositional")
##############################################
#devtools::install_github("microbiome/microbiome")
library(microbiome)
library(RColorBrewer)
library(reshape)
library(eulerr)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)


# test with stalks first
stalksR # 660 taxa 79 samples - these are below 50% so ignore it
head(prevalence(stalksR, detection = 1/100, sort = TRUE)) # Most Abundant OTUs and what percentage samples they are in
head(prevalence(stalksR, detection = 1/100, sort = TRUE, count = TRUE)) # How many samples are they actually in 

#rhizos
rhizosR
head(prevalence(rhizosR, detection = 1/100, sort = TRUE)) # Most Abundant OTUs and what percentage samples they are in
head(prevalence(rhizosR, detection = 1/100, sort = TRUE, count = TRUE)) # How many samples are they actually in 
#rhizos.core <- core(rhizosR, detection = 1/100, prevalence = .5)

prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

prevalences <- seq(.5, 1, .1)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

plot_core(rhizosR, plot.type = "heatmap", colours = gray,
          prevalences = prevalences, detections = detections, min.prevalence = .5) +
  labs(x = "Detection Threshold (Relative Abundance (%))")

#roots
rootsF
head(prevalence(rootsF, detection = 1/100, sort = TRUE)) # Most Abundant OTUs and what percentage samples they are in
head(prevalence(rootsF, detection = 1/100, sort = TRUE, count = TRUE)) # How many samples are they actually in 


#### Venn Diagram
library("ggVennDiagram")
#https://microbiome.github.io/tutorials/core_venn.html
# how many samples for metadata groups

### Experiment
# For ultra filtered PhyCmbFilt - only taxa that they all share
table(meta(phyCmbComp)$Experiment)
ExperimentGroups <- unique(as.character(meta(phyCmbComp)$Experiment))
print(ExperimentGroups)
list_core <- c() # empty object
for (n in ExperimentGroups){
  print(as.name(n))
  core.sub <- subset_samples(phyCmbComp, Experiment == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

GH_core_tax = as.data.frame(tax_table(phyCmbComp)[list_core$GH, ], row.names = NULL)
END_core_taxa = as.data.frame(tax_table(phyCmbComp)[list_core$END, ], row.names = NULL)
MMH_core_taxa = as.data.frame(tax_table(phyCmbComp)[list_core$MMH, ], row.names = NULL)

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Ultra Filtered Core Microbiome")

# How to identify the actual taxa!!!!!
itercept1 = merge(GH_core_tax, MMH_core_taxa, by="row.names")
merge(GH_core_tax, END_core_taxa, by="row.names")


# For less filtered phy2com - all taxa
table(meta(phy2com)$Experiment)
ExperimentGroups <- unique(as.character(meta(phy2com)$Experiment))
print(ExperimentGroups)
list_core <- c() # empty object
for (n in ExperimentGroups){
  print(as.name(n))
  core.sub <- subset_samples(phy2com, Experiment == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .0001,            # had to decrease stringency because there were so many taxa
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
# Get taxa info
GH_core_tax = as.data.frame(tax_table(phy2com)[list_core$GH, ], row.names = NULL)

END_core_taxa = as.data.frame(tax_table(phy2com)[list_core$END, ], row.names = NULL)

MMH_core_taxa = as.data.frame(tax_table(phy2com)[list_core$MMH, ], row.names = NULL)

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "forestgreen") + ggtitle("All Taxa Core Microbiome")

itercept1 = merge(GH_core_tax, MMH_core_taxa, by="row.names")


### By Inbred/Hybrid
# For ultra filtered PhyCmbFilt - only taxa that they all share
table(meta(phyCmbComp)$Inbred_or_Hybrid)
#BackgroundGroups <- unique(as.character(meta(phyCmbComp)$Inbred_or_Hybrid))
BackgroundGroups <- c("Inbred", "Open_Pollinated", "Hybrid")      # don't want soil or blanks
print(BackgroundGroups)
list_core <- c() # empty object
for (n in BackgroundGroups){
  print(as.name(n))
  core.sub <- subset_samples(phyCmbComp, BackgroundGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
#Inbred_core_tax = data.frame()

# Not sure what this issue is but had to change to this format and it takes a while
Inbred_core_tax = as.data.frame(phyCmbComp@tax_table)[list_core$Inbred, ]
Open_Pollinated_core_tax = as.data.frame(phyCmbComp@tax_table)[list_core$Open_Pollinated, ]
Hybrid_core_tax = as.data.frame(phyCmbComp@tax_table)[list_core$Hybrid, ]

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Ultra Filtered Core Microbiome")

# How to identify the actual taxa!!!!!
itercept1 = merge(GH_core_tax, MMH_core_taxa, by="row.names")
merge(GH_core_tax, END_core_taxa, by="row.names")


# # For less filtered phy2com - all taxa
table(meta(phy2com)$Inbred_or_Hybrid)
#BackgroundGroups <- unique(as.character(meta(phyCmbComp)$Inbred_or_Hybrid))
BackgroundGroups <- c("Inbred", "Open_Pollinated", "Hybrid")      # don't want soil or blanks
print(BackgroundGroups)
list_core <- c() # empty object
for (n in BackgroundGroups){
  print(as.name(n))
  core.sub <- subset_samples(phy2com, BackgroundGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .0001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
#Inbred_core_tax = data.frame()

# Not sure what this issue is
Inbred_core_tax = as.data.frame(phy2com@tax_table)[list_core$Inbred, ]
Open_Pollinated_core_tax = as.data.frame(phy2com@tax_table)[list_core$Open_Pollinated, ]
Hybrid_core_tax = as.data.frame(phy2com@tax_table)[list_core$Hybrid, ]

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "forestgreen") + ggtitle("All Taxa Core Microbiome")

# How to identify the actual taxa!!!!!
itercept1 = merge(GH_core_tax, MMH_core_taxa, by="row.names")
merge(GH_core_tax, END_core_taxa, by="row.names")


### By Sample Type
# For ultra filtered PhyCmbFilt - only taxa that they all share
table(meta(phyCmbComp)$Sample_Type)
#BackgroundGroups <- unique(as.character(meta(phyCmbComp)$Inbred_or_Hybrid))
SampleGroups <- c("Rhizosphere", "Root", "Stalk")      # don't want soil or blanks
print(SampleGroups)
list_core <- c() # empty object
for (n in SampleGroups){
  print(as.name(n))
  core.sub <- subset_samples(phyCmbComp, SampleGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
#Inbred_core_tax = data.frame()

# Not sure what this issue is
Rhizo_core_tax = as.data.frame(phyCmbComp@tax_table)[list_core$Rhizosphere, ]
Roots_Pollinated_core_tax = as.data.frame(phyCmbComp@tax_table)[list_core$Root, ]
Stalks_core_tax = as.data.frame(phyCmbComp@tax_table)[list_core$Stalk, ]

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Ultra Filtered Core Microbiome")

# How to identify the actual taxa!!!!!
itercept1 = merge(GH_core_tax, MMH_core_taxa, by="row.names")
merge(Stalks_core_tax, Rhizo_core_tax, by="row.names")


# # For less filtered phy2com - all taxa
table(meta(phy2com)$Sample_Type)
#BackgroundGroups <- unique(as.character(meta(phyCmbComp)$Inbred_or_Hybrid))
SampleGroups <- c("Rhizosphere", "Root", "Stalk")      # don't want soil or blanks
print(SampleGroups)
list_core <- c() # empty object
for (n in SampleGroups){
  print(as.name(n))
  core.sub <- subset_samples(phy2com, SampleGroups == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .0001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
#Inbred_core_tax = data.frame()

# Not sure what this issue is
Rhizo_core_tax = as.data.frame(phy2com@tax_table)[list_core$Rhizosphere, ]
Roots_Pollinated_core_tax = as.data.frame(phy2com@tax_table)[list_core$Root, ]
Stalks_core_tax = as.data.frame(phy2com@tax_table)[list_core$Stalk, ]

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "forestgreen") + ggtitle("All Taxa Core Microbiome")






