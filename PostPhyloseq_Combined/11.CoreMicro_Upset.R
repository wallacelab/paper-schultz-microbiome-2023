# Make an Upset Plot of the Core Microbiome in each compartment for each experiment

# Take pre proccessed qiime objects, filter with phyloseq, then insert back into qiime for diversity metrics!
# Why did I create a second version of this?

library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)
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

phyCmbFilt <- tax_glom(phyCmbFilt, taxrank = "Genus")
phy2 <- tax_glom(phy2, taxrank = "Genus")

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
library(HMP2Data)
library(phyloseq)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(ggplot2)
library(UpSetR)
library(tidyr)
library(data.table)
library(ggVennDiagram)

phy2com

mergedGroups <- merge_samples(phy2com, "Group", fun = sum)

upset_object <- as.data.frame(t(otu_table(mergedGroups)))

binary_full <- upset_object

# Make it all zeros

zero_fun <- function(x){
  return(x*0)
}

binary_full[] <- lapply(binary_full,zero_fun)

table(meta(phy2com)$Group)
ExperimentGroups <- unique(as.character(meta(phy2com)$Group))

list_core <-c() #empty object
for (n in ExperimentGroups){
  print(as.name(n))
  core.sub <- subset_samples(phy2com, Group == n)
 
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.9)  
  list_core[[n]] <- core_m
  
  for(i in core_m){
    #print(i)
    binary_full[i,toString(n)] <- 1
  }
}

binary_full 

upset_order <- colnames(binary_full)

noblanks = subset(binary_full[c(9:13,15:17)])
undergroundmatrix = subset(binary_full[c(9:13)])
above_matrix = subset(binary_full[c(15:17)])

# table(meta(phy2com)$Group)
# ExperimentGroups <- unique(as.character(meta(phy2com)$Group))

all_upset <- upset(noblanks, nsets = 8, nintersects = NA, order.by = c("degree", "freq"))
all_upset

under_upset <- upset(undergroundmatrix, nsets = 5, nintersects = NA, order.by = c("degree", "freq"))
under_upset

above_upset <- upset(above_matrix, nsets = 3, nintersects = NA, order.by = c("degree", "freq"), empty.intersections = TRUE)
above_upset

## Get taxa from all 5 underground groups!
undergroundmatrix

# tidyverse
library(dplyr)
core_otus <- undergroundmatrix %>% filter_all(all_vars(.>0))
nrow(core_otus)

core_tax <- rownames(core_otus)

core_under_phy <- subset(otu_table(phy2com), rownames(otu_table(phy2com)) %in% 
                           c(core_taxa))
core_under_res <- merge_phyloseq(core_under_phy, tax_table(phy2com), sample_data(phy2com))

core2write <- tax_table(core_under_res)

write.csv(core2write,file = "core_underground_taxa.csv", row.names = TRUE)
getwd()

### Explore whether inbreds and their hybrid cross have similar core microbiomes

# Compare B73xMO17 in the greenhouse and field
B73xMO17_Family = c("B73","Mo17","B73xMo17","Mo17xB73")
GH_cross1 = subset_samples(GH, Genotype %in% B73xMO17_Family)
GH_cross1_com <- microbiome::transform(GH_cross1, "compositional")
Genos <- unique(as.character(meta(GH_cross1_com)$Genotype))
print(Genos)

list_core <-c() #empty object
for (n in Genos){
  print(as.name(n))
  core.sub <- subset_samples(GH_cross1_com, Genotype == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing shared Inbred and Hybrid Taxa: Greenhouse")


END_cross1 = subset_samples(END, Genotype %in% B73xMO17_Family)
END_cross1_com <- microbiome::transform(END_cross1, "compositional")
Genos <- unique(as.character(meta(END_cross1_com)$Genotype))
print(Genos)

list_core <-c() #empty object
for (n in Genos){
  print(as.name(n))
  core.sub <- subset_samples(GH_cross1_com, Genotype == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing shared Inbred and Hybrid Taxa: END")

# Compare B73 and PH207 in the greenhouse and field
B73xPH207_Family = c("B73","Ph207","B73xPh207","Ph207xB73")
GH_cross2 = subset_samples(GH, Genotype %in% B73xPH207_Family)
GH_cross2_com <- microbiome::transform(GH_cross2, "compositional")
Genos <- unique(as.character(meta(GH_cross2_com)$Genotype))
print(Genos)

list_core <-c() #empty object
for (n in Genos){
  print(as.name(n))
  core.sub <- subset_samples(GH_cross2_com, Genotype == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing shared Inbred and Hybrid Taxa: Greenhouse")

# MMH
B73xPH207_Family = c("B73","Ph207","B73xPh207","Ph207xB73")
MMH_cross2 = subset_samples(MMH, Genotype %in% B73xPH207_Family)
MMH_cross2_com <- microbiome::transform(MMH_cross2, "compositional")
Genos <- unique(as.character(meta(MMH_cross2_com)$Genotype))
print(Genos)

list_core <-c() #empty object
for (n in Genos){
  print(as.name(n))
  core.sub <- subset_samples(MMH_cross2_com, Genotype == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing shared Inbred and Hybrid Taxa: MMH")

# Other GH Crosses
# Mo17xPh207
B73xPH207_Family = c("Mo17","Ph207","Mo17xPh207","Ph207xMo17")
GH_cross2 = subset_samples(GH, Genotype %in% B73xPH207_Family)
GH_cross2_com <- microbiome::transform(GH_cross2, "compositional")
Genos <- unique(as.character(meta(GH_cross2_com)$Genotype))
print(Genos)

list_core <-c() #empty object
for (n in Genos){
  print(as.name(n))
  core.sub <- subset_samples(GH_cross2_com, Genotype == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing shared Inbred and Hybrid Taxa: Greenhouse")

# B73 x CML247
B73xPH207_Family = c("B73","CML247","B73xCML247","CML247xB73")
GH_cross2 = subset_samples(GH, Genotype %in% B73xPH207_Family)
GH_cross2_com <- microbiome::transform(GH_cross2, "compositional")
Genos <- unique(as.character(meta(GH_cross2_com)$Genotype))
print(Genos)

list_core <-c() #empty object
for (n in Genos){
  print(as.name(n))
  core.sub <- subset_samples(GH_cross2_com, Genotype == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing shared Inbred and Hybrid Taxa: Greenhouse")

# B73 x Oh43
B73xPH207_Family = c("B73","Oh43","B73xOh43","Oh43xB73")
GH_cross2 = subset_samples(GH, Genotype %in% B73xPH207_Family)
GH_cross2_com <- microbiome::transform(GH_cross2, "compositional")
Genos <- unique(as.character(meta(GH_cross2_com)$Genotype))
print(Genos)

list_core <-c() #empty object
for (n in Genos){
  print(as.name(n))
  core.sub <- subset_samples(GH_cross2_com, Genotype == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)          # 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}

ggVennDiagram(list_core) + scale_fill_gradient(low = "white", high = "firebrick") + ggtitle("Comparing shared Inbred and Hybrid Taxa: Greenhouse")






