# Beta Diversity Figure and PERMANOVA
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

######## Load Data
getwd()
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/Phyloseq_Objects")

# Use main data - phy_data

# import 
metadata <- data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ","))
# 
# otu = otu_table(as.matrix(read.table("phy_data_otu.csv", row.names = 1, header = TRUE, sep = ",")),taxa_are_rows = TRUE)
# taxa = tax_table(as.matrix(read.table("phy_data_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")))
# meta = sample_data(data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ",")), errorIfNULL = FALSE)
# tree = read.tree("phy_data_tree.tree")


otu = read.table("phy_data_otu.csv", row.names = 1, header = TRUE, sep = ",")
taxa = read.table("phy_data_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")
meta = sample_data(data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ",")), errorIfNULL = FALSE)
tree = read.tree("phy_data_tree.tree")

# replace "." with "-" otu dfs
colnames(otu)
colnames(otu) <- gsub("\\.", "-", colnames(otu))
colnames(otu)

df <- otu %>%
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) 

df2 <- df[,-1]
row.names(df2) <- df[,1]

remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)


# Transpose data frame as a dataframe, merge metadata, make long!


