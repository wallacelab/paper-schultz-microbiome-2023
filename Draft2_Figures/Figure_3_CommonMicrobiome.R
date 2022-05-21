# Common Microbiome Figure
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
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(ggpubr)
library(rbiom)
library(microbiome)
library(RColorBrewer)
library(reshape)
library(eulerr)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)
library(HMP2Data)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(UpSetR)
library(data.table)
library(ggVennDiagram)


######## Load Data
getwd()
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/Phyloseq_Objects")

# Use main data - phy_data

# import 
metadata <- data.frame(read.table("Group_Metadata.csv", row.names = 1, header = TRUE, sep = ","))

otu = otu_table(as.matrix(read.table("phy2_otu.csv", row.names = 1, header = TRUE, sep = ",")),taxa_are_rows = TRUE)
taxa = tax_table(as.matrix(read.table("phy2_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")))
meta = sample_data(data.frame(read.table("phy2_sampdat.csv", row.names = 1, header = TRUE, sep = ",")))
tree = read.tree("phy2_tree.tree")

# replace "." with "-" otu dfs
colnames(otu)
colnames(otu) <- gsub("\\.", "-", colnames(otu))
colnames(otu)


phy_data = phyloseq(otu,taxa, meta, tree)
phy_data = merge_phyloseq(phy_data, meta, tree)

sample_sums(subset_samples(phy_data, Sample_Type=="Stalk"))
sample_sums(subset_samples(phy_data, Sample_Type=="Root"))
sample_sums(subset_samples(phy_data, Sample_Type=="Rhizosphere"))

phy_data = prune_samples(sample_sums(phy_data)>=10,phy_data)

phy_data

# Start core microbiome


##################### Separate tissue by inbred or hybrid and do that
phy_data

phy_tissue <- subset_samples(phy_data,Sample_Type_Blanks_differentiated=="Stalk" |
                                        Sample_Type_Blanks_differentiated=="Rhizosphere" |
                                        Sample_Type_Blanks_differentiated=="Root") 
phy_tissue
phy_tissue <- subset_samples(phy_tissue, Inbred_or_Hybrid != "Open_Pollinated")
sample_data(phy_tissue)$T_and_B <- paste(sample_data(phy_tissue)$Inbred_or_Hybrid,
                                         sample_data(phy_tissue)$Sample_Type_Blanks_differentiated)

phy_tissue <- tax_glom(phy_tissue, taxrank = "Genus") 

phy_tissue_comp <- microbiome::transform(phy_tissue, "compositional")

mergedGroups <- merge_samples(phy_tissue_comp, "T_and_B", fun = sum)

upset_object <- as.data.frame(t(otu_table(mergedGroups)))

binary_full <- upset_object

# Make it all zeros

zero_fun <- function(x){
  return(x*0)
}

binary_full[] <- lapply(binary_full,zero_fun)

table(meta(phy_tissue_comp)$T_and_B)
ExperimentGroups <- unique(as.character(meta(phy_tissue_comp)$T_and_B))

list_core <-c() #empty object
for (n in ExperimentGroups){
  print(as.name(n))
  core.sub <- subset_samples(phy_tissue_comp, T_and_B == n)
  
  core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
                         detection = .001,            # .001 in atleast 90% of samples
                         prevalence = 0.5)  
  list_core[[n]] <- core_m
  
  for(i in core_m){
    #print(i)
    binary_full[i,toString(n)] <- 1
  }
}

binary_full 

#Spiffy it up with ComplexUpset instead
library(ComplexUpset)
Tissue_and_Background <- colnames(binary_full)

## Shared table
tax_table(phy_tissue_comp)

taxa_matrix <- as.data.frame(as(tax_table(phy_tissue_comp), "matrix"))


shared_df <- merge(binary_full,taxa_matrix, by = 'row.names', all = FALSE)

# Test without colors and queiries
#complex <- ComplexUpset::upset(binary_full, Tissue_and_Background, name = "Tissue and Background Genus")
the_intersections = list(c("Inbred Rhizosphere","Hybrid Rhizosphere", "Inbred Root", "Hybrid Root"),
                         c('Hybrid Rhizosphere','Inbred Rhizosphere'),
                         c('Hybrid Root','Inbred Root'),
                         c('Hybrid Stalk','Inbred Stalk'),
                         c('Hybrid Root','Hybrid Rhizosphere'),
                         c('Inbred Root','Inbred Rhizosphere'))

complex <- ComplexUpset::upset(shared_df, Tissue_and_Background, name = "Tissue and Background Common Genus",
                               max_size = 100, min_degree=2,width_ratio=0.1, intersections = the_intersections,base_annotations=list(
                                 'Intersection size'=intersection_size(
                                   text_colors=c(
                                     on_background='black', on_bar='white'
                                   )
                                 )
                                 + annotate(geom='text', x=1, y=25,label='All Underground', color = 'white', angle = 90, size = 6)
                                 + annotate(geom='text', x=2, y=25,label='Inbred Underground', color = 'royalblue3', angle = 90, size = 6)
                                 + annotate(geom='text', x=3, y=25,label='Hybrid Underground', color = 'firebrick', angle = 90, size = 6)
                                 + annotate(geom='text', x=4, y=25,label='Rhizosphere', color = 'purple', angle = 90, size = 6)
                                 + annotate(geom='text', x=5, y=25,label='Root', color = 'tan3', angle = 90, size = 6)
                                 + annotate(geom='text', x=6, y=25,label='Stalk', color = 'olivedrab', angle = 90, size = 6)),
                               queries = list(
                                 upset_query(
                                   intersect = c('Hybrid Rhizosphere',
                                                 'Inbred Rhizosphere'), color = 'purple', fill = 'purple'),
                                 upset_query(
                                   intersect = c('Hybrid Root',
                                                 'Inbred Root'), color = 'tan4', fill = 'tan4'),
                                 upset_query(
                                   intersect = c('Hybrid Stalk',
                                                 'Inbred Stalk'), color = 'olivedrab', fill = 'olivedrab'),
                                 upset_query(
                                   intersect = c('Hybrid Root',
                                                 'Hybrid Rhizosphere'), color = 'firebrick', fill = 'firebrick'),
                                 upset_query(
                                   intersect = c('Inbred Root',
                                                 'Inbred Rhizosphere'), color = 'royalblue3', fill = 'royalblue3'),
                                 upset_query(set='Inbred Rhizosphere', fill='purple'), upset_query(set='Hybrid Rhizosphere', fill='purple'),
                                 upset_query(set='Inbred Root', fill='tan4'), upset_query(set='Hybrid Root', fill='tan4'),
                                 upset_query(set='Inbred Stalk', fill='olivedrab'), upset_query(set='Hybrid Stalk', fill='olivedrab')
                               )
                               # , annotations = list(
                               #   'Phylum Breakdown'=(
                               #     ggplot(mapping=aes(fill=Phylum))
                               #     + geom_bar(stat='count', position='fill')
                               #     + scale_y_continuous(labels=scales::percent_format())
                               #   )
                               #   + ylab('Phylum Percentage') 
                               # # This is for the phylum breakdown
                               ) 
complex

ggsave("Fig3_Common_D2.png",
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
       complex, device = "png", width = 10, height = 10, dpi = 600)


# ###########################################
# 
# upset_order <- colnames(binary_full)
# 
# tissue_and_background_upset <- UpSetR::upset(binary_full, nsets = 9, nintersects = NA, order.by = "freq")
# tissue_and_background_upset
# 
# phy_data <- tax_glom(phy_data, taxrank = "Genus") 
# phy_comp <- microbiome::transform(phy_data, "compositional")
# 
# mergedGroups <- merge_samples(phy_comp, "Group", fun = sum)
# 
# upset_object <- as.data.frame(t(otu_table(mergedGroups)))
# 
# binary_full <- upset_object
# 
# # Make it all zeros
# 
# zero_fun <- function(x){
#   return(x*0)
# }
# 
# binary_full[] <- lapply(binary_full,zero_fun)
# 
# table(meta(phy_comp)$Group)
# ExperimentGroups <- unique(as.character(meta(phy_comp)$Group))
# 
# list_core <-c() #empty object
# for (n in ExperimentGroups){
#   print(as.name(n))
#   core.sub <- subset_samples(phy_comp, Group == n)
#   
#   core_m <- core_members(core.sub,                   # core.sub is only samples in experiment
#                          detection = .001,            # .001 in atleast 90% of samples
#                          prevalence = 0.9)  
#   list_core[[n]] <- core_m
#   
#   for(i in core_m){
#     #print(i)
#     binary_full[i,toString(n)] <- 1
#   }
# }
# 
# binary_full 
# 
# upset_order <- colnames(binary_full)
# 
# noblanks = subset(binary_full[c(9:13,15:17)])
# undergroundmatrix = subset(binary_full[c(9:13)])
# above_matrix = subset(binary_full[c(15:17)])
# 
# # table(meta(phy2com)$Group)
# # ExperimentGroups <- unique(as.character(meta(phy2com)$Group))
# 
# all_upset <- UpSetR::upset(noblanks, nsets = 8, nintersects = NA, order.by = "freq")
# all_upset