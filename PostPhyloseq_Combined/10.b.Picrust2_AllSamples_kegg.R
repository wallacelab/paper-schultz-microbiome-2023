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

SVs <- read_qza("Combined_qza_files/Combined_deblur_table.qza")

taxonomy <-read_qza("Combined_Results/TaxonomyStuff/Combined.taxonomy.qza")
taxtable <-taxonomy$data %>% as.data.frame() 
taxtable <- taxtable %>% separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") #convert the table into a tabular split version

tree <- read_qza("Combined_qza_files/Combined_rooted_tree.qza")

phy <- qza_to_phyloseq("Combined_qza_files/Combined_deblur_table.qza", "Combined_qza_files/Combined_rooted_tree.qza", "Combined_Results/TaxonomyStuff/Combined.taxonomy.qza", "Combined_Key.csv")
phy

################################################################ no need to run filtering more than once
# Line 200 to start
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

#Split by samples then remove blank OTUs and do phyloseq to Picrust file formats:

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

phyCmbFilt1 <- prune_taxa(goodlist, phy2)

Rhizegood <-  prune_taxa(taxa_sums(blank_rhizos) < 1, blank_rhizos)
Rgoodlist <- taxa_names(Rhizegood)
phyCmbFilt2 <- prune_taxa(Rgoodlist, phyCmbFilt1)

Rootgood <-  prune_taxa(taxa_sums(blank_roots) < 1, blank_roots)
Rootgoodlist <- taxa_names(Rootgood)
phy2noblank <- prune_taxa(Rootgoodlist, phyCmbFilt2)

# Cleaned combined data - maybe even too cleaned ya know - export a biome file
phyCmbFiltClean

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Picrust2")

phyloseq2qiime2<-function(physeq){
  #take a phyloseq object,check for individual parts, write to files ready for qiime2 upload
  library(phyloseq)
  library(biomformat)
  library(ape)
  library(Biostrings)
  library(dada2)
  if(packageVersion("biomformat") < "1.7") {
    stop("This will only work with biomformat version > 1.7")
  }
  ps_name <-deparse(substitute(physeq))
  taxa_are_rows_logical<-taxa_are_rows(physeq)
  #write OTU table to biom file
  if(is.null(access(physeq,"otu_table"))==FALSE){
    if(taxa_are_rows_logical==TRUE) {
      otu<-as(otu_table(physeq),"matrix")
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_features-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    } else if (taxa_are_rows_logical==FALSE) {
      otu<-t(as(otu_table(physeq),"matrix"))
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_feature-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    }
  }
  #write sample data (metadata) to tsv
  if(is.null(access(physeq,"sam_data"))==FALSE){
    write.table(sample_data(physeq),file=paste0(ps_name,"_sample-metadata.txt"),
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    print(paste0("Writing sample metadata to ",ps_name,"_sample-metadata.txt"))
  }
  #write taxonomy table to qiime2 formatted taxonomy
  if(is.null(access(physeq,"tax_table"))==FALSE){
    tax<-as(tax_table(physeq),"matrix")
    tax_cols <- colnames(tax)
    tax<-as.data.frame(tax)
    tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
    for(co in tax_cols) tax[co]<-NULL
    write.table(tax, file=paste0(ps_name,"_tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
    print(paste0("Writing taxonomy table to ",ps_name,"_tax.txt"))
  }
  #write phylogenetic tree to newick formwat
  if(is.null(access(physeq,"phy_tree"))==FALSE){
    if(is.rooted(phy_tree(physeq))==TRUE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-rooted.newick"))
      print(paste0("Writing rooted tree to ",ps_name,"_tree-rooted.newick"))
    } else if (is.rooted(phy_tree(physeq))==FALSE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-unrooted.newick"))
      print(paste0("Writing unrooted tree to ",ps_name,"_tree-unrooted.newick"))
    }
  }
  #write representative sequences to fasta format
  if(is.null(access(physeq,"refseq"))==FALSE){
    writeXStringSet(refseq(physeq),filepath=paste0(ps_name,"_ref-seqs.fasta"))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  } else if (taxa_are_rows_logical==FALSE && unique(grepl("[^ATCG]",colnames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(t(otu), fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  } else if (taxa_are_rows_logical==TRUE && unique(grepl("[^ATCG]",rownames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(otu, fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  }
}

# Export to qiime qza and also a biome file
phyloseq2qiime2(phy2)
phyloseq2qiime2(phy2noblank)

# biome file = phyCmbFiltClean_features-table.biom

# fasta table = dna-sequences.fasta from /home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/Combined_rep_seqs.qza

# metadata is metadata



### Command line -> conda activate picrust2
### picrust2_pipeline.py -s dna-sequences.fasta -i phyCmbFiltClean_features-table.biom -o picrust2_out_pipeline -p 4
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Picrust2")

KO_table <- read.csv("KO_pred_metagenome_unstrat_descript.tsv", sep = "\t")
head(KO_table)[1:10]

# Create a phyloseq object out of the otu table and the metadata. 
#head(EC_table)[1:10]
# E_C table column names have . instead of - 
colnames(KO_table) <- gsub(x = colnames(KO_table), pattern = "\\.", replacement = "-")
Meta_EC <- metadata
row.names(Meta_EC) <- Meta_EC$SampleID
row.names(KO_table) <- KO_table$`function-`
# Drop the descriptors
descriptions = KO_table[,0:2]
EC_table <- KO_table[, -c(1:2)]
EC_phy <- phyloseq(otu_table(EC_table, taxa_are_rows = TRUE), sample_data(Meta_EC))
EC_phy # this is our functional phyloseq object


### try to categorize by function similar to picrust1
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Picrust2")

KO_table <- read.csv("KO_pred_metagenome_unstrat_descript.tsv", header=TRUE, sep="\t", row.names=1)
head(KO_table)[1:10]

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
table_ko_L3 <- categorize_by_function_l3(KO_table, brite_map_trim)
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
    scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1")) + coord_flip()
  
  return(plot)
}  
#### All Tissues
# Inbred vs hybrid
Diff_abun_func_IvH(EC_phy, ~Inbred_or_Hybrid, 0.001)
plot_diff_funct_IvH(EC_phy, ~Inbred_or_Hybrid, " All Tissues: Inbred vs Hybrid")

# Location
Diff_abun_func_location(EC_phy, ~Location, 0.001)
plot_diff_funct_location(EC_phy, ~Location, " All Tissues: Inbred vs Hybrid/Open Pollinated")


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



# Cleaned up all this code with functions - should i delete this all?


# #All tissues inbed vs. hybrid
# alltis_ect = transform_sample_counts(EC_phy, function(OTU) OTU + 1)
# deAll = phyloseq_to_deseq2(alltis_ect, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deAll, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- (c(rownames(sigtab)))
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
#   geom_bar(stat = 'identity') + ggtitle("All Tissue: Inbred vs Hybrid/Open Pollinated") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# 
# # All tissues Field vs GH
alltis_ect = transform_sample_counts(EC_phy, function(OTU) OTU + 1)
deAll = phyloseq_to_deseq2(alltis_ect, ~ Location)
deStalk = DESeq(deAll, test="Wald", fitType = "parametric")

res = results(deStalk, cooksCutoff = FALSE, contrast = c("Location", "GH", "IH"))
alpha = 0.001
sigtab = res[which(res$padj < alpha), ]
Functional_Group <- (c(rownames(sigtab)))
sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
head(sigtab)

dim(sigtab)
# 
# ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
#   geom_bar(stat = 'identity') + ggtitle("All Tissue: Inbred vs Hybrid/Open Pollinated") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# 
# # Stalks - inbred vs hybrid
# stalks_ECt = transform_sample_counts(stalks_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
# geom_bar(stat = 'identity') + ggtitle("Stalks: Inbred vs Hybrid/Open Pollinated") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# 
# # Rhizos - inbred vs hybrid
# stalks_ECt = transform_sample_counts(rhizos_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
#   geom_bar(stat = 'identity') + ggtitle("Rhizos: Inbred vs Hybrid/Open Pollinated") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# 
# 
# # Roots - inbred vs hybrid
# stalks_ECt = transform_sample_counts(roots_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
#   geom_bar(stat = 'identity') + ggtitle("Roots: Inbred vs Hybrid/Open Pollinated") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# 
# 
# ####### Notes: Interestingly when I subset by tissue type and then do deseq: 
# # even with alpha = .01, stalk is the only compartment to experience different functional groups?
# 
# stalks_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Stalk")
# rhizos_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Rhizosphere")
# roots_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Root")
# # Deseq
# 
# # Stalks - location
# stalks_ECt = transform_sample_counts(stalks_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Location)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
#   geom_bar(stat = 'identity') + ggtitle("Stalks: Field vs GH") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# 
# # Rhizos - Location
# stalks_ECt = transform_sample_counts(rhizos_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Location)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# # idk why we cant just plot sigtab now so here it goes
# st <- as.data.frame(sigtab)
# st['Functional_Group'] = row.names(st)
# st <- base::transform(st, Functional_Group = reorder(Functional_Group, log2FoldChange))
# ggplot(st, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
#   geom_bar(stat = 'identity') + ggtitle("Rhizos: Inbred vs Hybrid/Open Pollinated") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1")) + coord_flip()
# 
# 
# # Roots - Location
# stalks_ECt = transform_sample_counts(roots_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Location)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# ggplot(sigtab, aes(x = Functional_Group, y = log2FoldChange, fill = log2FoldChange < 0)) + 
#   geom_bar(stat = 'identity') + ggtitle("Roots: Inbred vs Hybrid/Open Pollinated") +
#   theme(axis.text.x = element_text(angle = 90, size = 12)) + 
#   scale_fill_manual("Down Regulated", values = c("turquoise", "indianred1"))
# 
# 
# ##################################################################################
# ##### Inbred and Hybrid for just one experiment at a time - GH
# stalks_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Stalk")
# rhizos_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Rhizosphere")
# roots_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Root")
# 
# stalks_EC <- subset_samples(stalks_EC, Experiment=="GH")
# rhizos_EC <- subset_samples(rhizos_EC, Experiment=="GH")
# roots_EC <- subset_samples(roots_EC, Experiment=="GH")
# 
# # Stalks - inbred vs hybrid
# stalks_ECt = transform_sample_counts(stalks_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# 
# # Rhizos - inbred vs hybrid
# stalks_ECt = transform_sample_counts(rhizos_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# # Roots - inbred vs hybrid
# stalks_ECt = transform_sample_counts(roots_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# 
# ##################################################################################
# ##### Inbred and Hybrid for just one experiment at a time - END
# stalks_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Stalk")
# rhizos_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Rhizosphere")
# roots_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Root")
# 
# stalks_EC <- subset_samples(stalks_EC, Experiment=="END")
# rhizos_EC <- subset_samples(rhizos_EC, Experiment=="END")
# roots_EC <- subset_samples(roots_EC, Experiment=="END")
# 
# # Stalks - inbred vs hybrid
# stalks_ECt = transform_sample_counts(stalks_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# 
# # Rhizos - inbred vs hybrid
# stalks_ECt = transform_sample_counts(rhizos_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# # Roots - inbred vs hybrid
# stalks_ECt = transform_sample_counts(roots_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# 
# ##################################################################################
# ##### Inbred and Hybrid for just one experiment at a time - MMH
# stalks_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Stalk")
# rhizos_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Rhizosphere")
# roots_EC <- subset_samples(EC_phy, Sample_Type_Blanks_differentiated=="Root")
# 
# stalks_EC <- subset_samples(stalks_EC, Experiment=="MMH")
# rhizos_EC <- subset_samples(rhizos_EC, Experiment=="MMH")
# #roots_EC <- subset_samples(roots_EC, Experiment=="MMH")
# 
# # Stalks - inbred vs hybrid
# stalks_ECt = transform_sample_counts(stalks_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# 
# # Rhizos - inbred vs hybrid
# stalks_ECt = transform_sample_counts(rhizos_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# # Roots - inbred vs hybrid
# stalks_ECt = transform_sample_counts(roots_EC, function(OTU) OTU + 1)
# deStalk = phyloseq_to_deseq2(stalks_ECt, ~ Inbred_or_Hybrid)
# deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")
# 
# res = results(deStalk, cooksCutoff = FALSE)
# alpha = 0.001
# sigtab = res[which(res$padj < alpha), ]
# Functional_Group <- as(descriptions[c(rownames(sigtab)),2], "matrix")
# sigtab = cbind(as(sigtab, "data.frame"), Functional_Group)
# head(sigtab)
# 
# dim(sigtab)
# 
# # So having these KOs get changed into functional groups is helpful, but I have a massive amount of false positives
# # Unfortunately, these are old KOs and one KO is in several groups, so I dont know how to create a custom 
# # mapping file in any real useful way. Obviously plants cant get parkinsons, but the same KO also counts for oxidative phosphorylation
# 
# # use grepl to remove things from kegg brite mapping file like human diseases above.
# 
# # Look at unique groups
# unique(unlist(strsplit(as.character(kegg_brite_map$metadata_KEGG_Pathways), ";")))







