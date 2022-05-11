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

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS")
# Read in the data
metadata <- read.csv("Combined_Key.csv", header = TRUE, sep = "\t")

### try to categorize by function similar to picrust1
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Picrust2")

KO_table <- read.csv("KO_pred_metagenome_unstrat_descript.tsv", header=TRUE, sep="\t", row.names=1)
phy2_KO_table <- read.csv("ph2_KOtable_descript.tsv", header=TRUE, sep="\t", row.names=1)
phy2NB_KO_table <- read.csv("ph2noblank_KOtable_descript.tsv", header=TRUE, sep="\t", row.names=1)

head(phy2_KO_table)[1:10]

kegg_brite_map <- read.table("picrust1_KO_BRITE_map.tsv",
                             header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)
#head(kegg_brite_map)
dim(kegg_brite_map)

# remove rows in brite map related to human diseases and a few other things
# Look at unique groups
unique(unlist(strsplit(as.character(kegg_brite_map$metadata_KEGG_Pathways), ";")))

brite_map_trim <- kegg_brite_map[!grepl("Human Diseases", kegg_brite_map$metadata_KEGG_Pathways),]
dim(brite_map_trim)



categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.
  
  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))
  
  colnames(out_pathway) <- c("pathway", colnames(in_ko))
  
  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}

### categorize by function similar to picrust1

### Run function to categorize all KOs by level 3 in BRITE hierarchy.
table_ko_L3 <- categorize_by_function_l3(phy2NB_KO_table, brite_map_trim)
# test_ko_L3_sorted <- test_ko_L3[rownames(orig_ko_L3), ]
head(table_ko_L3)[1:10]
ko_l3 <- table_ko_L3
# Now do the phyloseq object and diff abundance. 
colnames(ko_l3) <- gsub(x = colnames(ko_l3), pattern = "\\.", replacement = "-")
#ko_l3 <- cbind(rownames(ko_l3), data.frame(ko_l3, row.names = NULL))
#ko_l3 <- subset(ko_l3, select = -c(description))
descriptions = row.names(ko_l3)
ko_l3 = ko_l3[, -c(1)]

head(ko_l3)
colnames(ko_l3)
Meta_EC <- metadata
row.names(Meta_EC) <- Meta_EC$SampleID
EC_phy <- phyloseq(otu_table(ko_l3, taxa_are_rows = TRUE), sample_data(Meta_EC))
EC_phy # this is our functional phyloseq object - call it ec even though its kegg for simplicity

gplots::venn(list(metadata=rownames(Meta_EC), featuretable=colnames(ko_l3)))

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
    scale_fill_manual("More abundant in:", labels = c("Inbred Plants","Hybrid Plants"),values = c("turquoise", "indianred1")) + coord_flip()
  
  return(plot)
}  


# Figure for paper is All Experiments Inbred vs Hybrid Roots, everything else comment out

Diff_abun_func_IvH(roots_EC, ~Inbred_or_Hybrid, 0.001)
Fig_5 <- plot_diff_funct_IvH(roots_EC, ~Inbred_or_Hybrid, "Roots: Inbred vs Hybrid Functional Predictions")

Fig_5 

ggsave("Fig4_DiffAbund.png", 
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/PaperFigures",
       Fig_5, device = "png", width = 9, height = 6, dpi = 600)


###################################
# #### All Tissues
# Inbred vs hybrid
Diff_abun_func_IvH(EC_phy, ~Inbred_or_Hybrid, 0.001)
plot_diff_funct_IvH(EC_phy, ~Inbred_or_Hybrid, " All Tissues: Inbred vs Hybrid")

# Location
Diff_abun_func_location(EC_phy, ~Location, 0.001)
plot_diff_funct_location(EC_phy, ~Location, " All Tissues: Greenhouse vs Field")


### All Experiments

# Inbred vs Hybrid
Diff_abun_func_IvH(stalks_EC, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(rhizos_EC, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(roots_EC, ~Inbred_or_Hybrid, 0.001)
plot_diff_funct_IvH(roots_EC, ~Inbred_or_Hybrid, "Roots: Inbred vs Hybrid Functional Predictions")

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

plot_diff_funct(All_T_GH, ~Inbred_or_Hybrid, " All Tissues in GH: Inbred vs Hybrid/Open Pollinated")

# END
All_T_END <- subset_samples(EC_phy, Experiment=="END")
stalks_END <- subset_samples(stalks_EC, Experiment=="END")
rhizos_END <- subset_samples(rhizos_EC, Experiment=="END")
roots_END <- subset_samples(roots_EC, Experiment=="END")

Diff_abun_func_IvH(All_T_END, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(stalks_END, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(rhizos_END, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(roots_END, ~Inbred_or_Hybrid, 0.001)

plot_diff_funct(roots_END, ~Inbred_or_Hybrid, " Roots in END: Inbred vs Hybrid/Open Pollinated")


# MMH
All_T_MMH <- subset_samples(EC_phy, Experiment=="MMH")
stalks_MMH <- subset_samples(stalks_EC, Experiment=="MMH")
rhizos_MMH <- subset_samples(rhizos_EC, Experiment=="MMH")

Diff_abun_func_IvH(All_T_MMH, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(stalks_MMH, ~Inbred_or_Hybrid, 0.001)
Diff_abun_func_IvH(rhizos_MMH, ~Inbred_or_Hybrid, 0.001)


plot_diff_funct(All_T_MMH, ~Inbred_or_Hybrid, " All tissues in MMH: Inbred vs Hybrid/Open Pollinated")



#### Figure for the paper
library(ggpubr)
p1 <- plot_diff_funct(stalks_EC, ~Location, "All Stalk Samples: Greenhouse vs Field")
p2 <- plot_diff_funct(rhizos_EC, ~Location, "All Rhizo Samples: Greenhouse vs Field")
p3 <- plot_diff_funct(roots_EC, ~Location, "All Root Samples: Greenhouse vs Field")

pf <- ggarrange(p1,p2,p3, labels = c("A","B","C"), ncol = 3, nrow =1)

# 
# 
# # Cleaned up all this code with functions - should i delete this all?
# 
# 
# # #All tissues inbed vs. hybrid
# # alltis_ect = transform_sample_counts(EC_phy, function(OTU) OTU + 1)
# # deAll = phyloseq_to_deseq2(alltis_ect, ~ Inbred_or_Hybrid)
# # deStalk = DESeq(deAll, test="Wald", fitType = "parametric")
# # 
# # res = results(deStalk, cooksCutoff = FALSE)
# # alpha = 0.001
# # sigtab = res[which(res$padj < alpha), ]
# # Functional_Group <- (c(rownames(sigtab)))
# # sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# # head(sigtab)
# # 
# # dim(sigtab)
# # 
# # ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
# #   geom_bar(stat = 'identity') + ggtitle("All Tissue: Inbred vs Hybrid/Open Pollinated") +
# #   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
# #   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# # 
# # # All tissues Field vs GH
# alltis_ect = transform_sample_counts(EC_phy, function(OTU) OTU + 1)
# deAll = phyloseq_to_deseq2(alltis_ect, ~ Location)
# deStalk = DESeq(deAll, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE, contrast = c("Location", "GH", "IH"))
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- (c(rownames(sigtab)))
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)