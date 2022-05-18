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

phy_data <- tax_glom(phy_data, taxrank = "Genus") 
phy_comp <- microbiome::transform(phy_data, "compositional")

mergedGroups <- merge_samples(phy_comp, "Group", fun = sum)

upset_object <- as.data.frame(t(otu_table(mergedGroups)))

binary_full <- upset_object

# Make it all zeros

zero_fun <- function(x){
  return(x*0)
}

binary_full[] <- lapply(binary_full,zero_fun)

table(meta(phy_comp)$Group)
ExperimentGroups <- unique(as.character(meta(phy_comp)$Group))

list_core <-c() #empty object
for (n in ExperimentGroups){
  print(as.name(n))
  core.sub <- subset_samples(phy_comp, Group == n)
  
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

all_upset <- upset(noblanks, nsets = 8, nintersects = NA, order.by = "freq")
all_upset









