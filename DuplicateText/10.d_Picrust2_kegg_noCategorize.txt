library(tidyverse)
library(devtools)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)

### using picrust table generated in 10.b - if I like it fill in the steps to get to picrust.py

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS")
metadata <- read.csv("Combined_Key.csv", header = TRUE, sep = "\t")

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Picrust2")
KO_table <- read.csv("KO_pred_metagenome_unstrat_descript.tsv", sep = "\t")
head(KO_table)[1:10]

# Create a phyloseq object out of the otu table and the metadata. 
colnames(KO_table) <- gsub(x = colnames(KO_table), pattern = "\\.", replacement = "-")
Meta_EC <- metadata
row.names(Meta_EC) <- Meta_EC$SampleID
row.names(KO_table) <- do.call(paste,c(KO_table[c("description","function-")],sep="-"))

# Drop function and description cols
KO_table <- KO_table[, -c(1:2)]

# Make phyloseq Object
EC_phy <- phyloseq(otu_table(KO_table, taxa_are_rows = TRUE), sample_data(Meta_EC))
EC_phy # this is our functional phyloseq object - call it ec even though its kegg for simplicity

# Deseq
library("DESeq2")
library("ggplot2")

stalks_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Stalk")
rhizos_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Rhizosphere")
roots_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Root")

###########################################################################
# Differential Abundance Functions:

Diff_abun_func_location <- function(phyloseq_obj, metadata_category, alpha_value){
  compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
  phydesq = phyloseq_to_deseq2(compartment_t, design = metadata_category)
  deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
  
  res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Location", "GH", "IH"))
  alpha = alpha_value
  sigtab = res[which(res$padj < alpha), ]
  Functional_Group <- (c(rownames(sigtab)))
  sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
  
  return(dim(sigtab))
}

plot_diff_funct_location <- function(phyloseq_obj, metadata_category, title){
  compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
  phydesq = phyloseq_to_deseq2(compartment_t, design = metadata_category)
  deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
  
  res = results(deseq_obj, cooksCutoff = FALSE,contrast = c("Location", "GH", "IH"))
  alpha = 0.001
  sigtab = res[which(res$padj < alpha), ]
  Functional_Group <- (c(rownames(sigtab)))
  sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
  
  st <- as.data.frame(sigtab)
  st['Functional_Group'] = row.names(st)
  st <- base::transform(st, Functional_Group = reorder(Functional_Group, log2FoldChange))
  plot <- ggplot(st, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
    geom_bar(stat = 'identity') + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, size = 12)) + 
    scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1")) + coord_flip()
  
  return(plot)
}

Diff_abun_func_IvH <- function(phyloseq_obj, metadata_category, alpha_value){
  compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
  phydesq = phyloseq_to_deseq2(compartment_t, design = metadata_category)
  deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
  
  res = results(deseq_obj, cooksCutoff = FALSE, contrast = c("Inbred_or_Hybrid", "Inbred", "Hybrid"))
  alpha = alpha_value
  sigtab = res[which(res$padj < alpha), ]
  Functional_Group <- (c(rownames(sigtab)))
  sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
  
  return(dim(sigtab))
}

plot_diff_funct_IvH <- function(phyloseq_obj, metadata_category, title){
  compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
  phydesq = phyloseq_to_deseq2(compartment_t, design = metadata_category)
  deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
  
  res = results(deseq_obj, cooksCutoff = FALSE,contrast = c("Inbred_or_Hybrid", "Inbred", "Hybrid"))
  alpha = 0.001
  sigtab = res[which(res$padj < alpha), ]
  Functional_Group <- (c(rownames(sigtab)))
  sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
  
  st <- as.data.frame(sigtab)
  st['Functional_Group'] = row.names(st)
  st <- base::transform(st, Functional_Group = reorder(Functional_Group, log2FoldChange))
  plot <- ggplot(st, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
    geom_bar(stat = 'identity') + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, size = 12)) + 
    scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1")) + coord_flip()
  
  return(plot)
}  



#### All Tissues
# Inbred vs hybrid
Diff_abun_func_IvH(EC_phy, ~Inbred_or_Hybrid, 0.001)
plot_diff_funct_IvH(EC_phy, ~Inbred_or_Hybrid, " All Tissues: Inbred vs Hybrid")

# Location
Diff_abun_func_location(EC_phy, ~Location, 0.001)
plot_diff_funct_location(EC_phy, ~Location, " All Tissues: Location")

### All Experiments
# Inbred vs Hybrid
Diff_abun_func_IvH(stalks_EC, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(rhizos_EC, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(roots_EC, ~Inbred_or_Hybrid, 0.001)

# plot all rhizos
plot_diff_funct_IvH(rhizos_EC, ~Inbred_or_Hybrid, " All Roots: IvsH")

# Location
Diff_abun_func_location(stalks_EC, ~Location, 0.001)
Diff_abun_func_location(rhizos_EC, ~Location, 0.001)
Diff_abun_func_location(roots_EC, ~Location, 0.001)
Diff_abun_func_location(EC_phy, ~Location, 0.001)

### By Experiment - Inbred vs Hybrid
# GH
All_T_GH <- subset_samples(EC_phy, Experiment=="GH")
stalks_GH <- subset_samples(stalks_EC, Experiment=="GH")
rhizos_GH <- subset_samples(rhizos_EC, Experiment=="GH")
roots_GH <- subset_samples(roots_EC, Experiment=="GH")

Diff_abun_func_IvH(All_T_GH, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(stalks_GH, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(rhizos_GH, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(roots_GH, ~Inbred_or_Hybrid, 0.001)

# END
All_T_END <- subset_samples(EC_phy, Experiment=="END")
stalks_END <- subset_samples(stalks_EC, Experiment=="END")
rhizos_END <- subset_samples(rhizos_EC, Experiment=="END")
roots_END <- subset_samples(roots_EC, Experiment=="END")

Diff_abun_func_IvH(All_T_END, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(stalks_END, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(rhizos_END, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(roots_END, ~Inbred_or_Hybrid, 0.001)

# MMH
All_T_MMH <- subset_samples(EC_phy, Experiment=="MMH")
stalks_MMH <- subset_samples(stalks_EC, Experiment=="MMH")
rhizos_MMH <- subset_samples(rhizos_EC, Experiment=="MMH")

Diff_abun_func_IvH(All_T_MMH, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(stalks_MMH, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(rhizos_MMH, ~Inbred_or_Hybrid, 0.001)















