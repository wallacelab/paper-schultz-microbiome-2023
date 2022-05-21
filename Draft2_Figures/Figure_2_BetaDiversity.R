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

otu = otu_table(as.matrix(read.table("phy_data_otu.csv", row.names = 1, header = TRUE, sep = ",")),taxa_are_rows = TRUE)
taxa = tax_table(as.matrix(read.table("phy_data_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")))
meta = sample_data(data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ",")), errorIfNULL = FALSE)
tree = read.tree("phy_data_tree.tree")

# replace "." with "-" otu dfs
colnames(otu)
colnames(otu) <- gsub("\\.", "-", colnames(otu))
colnames(otu)


phy_data = phyloseq(otu,taxa, meta, tree)
phy_data = merge_phyloseq(phy_data, meta, tree)
phy_data


sample_sums(subset_samples(phy_data, Sample_Type=="Stalk"))
sample_sums(subset_samples(phy_data, Sample_Type=="Root"))
sample_sums(subset_samples(phy_data, Sample_Type=="Rhizosphere"))

phy_data = prune_samples(sample_sums(phy_data)>=10,phy_data)

phy_data

set.seed(18)   # Rarefy and relative abundance

#phy_data_rare <- rarefy_even_depth(phy_data, sample.size = min(sample_sums(phy_data)), replace = TRUE)
phy_data_rare <- rarefy_even_depth(phy_data, sample.size = 50, replace = FALSE)

phy_data_rel <- transform_sample_counts(phy_data, function(x) ((x / sum(x)) * 10000))
phy_data_rel <- transform_sample_counts(phy_data, function(x) ((ceiling(x))))
##############################################

# Beta Diversity 

##############################################
# Generate Unifrac dfs from qiime
# Dont want to use phyloseq unifrac, so instead convert to biome file and use rbiome.

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/FinalFigs_qza")

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

######################
# Data sets 
phy_data
phy_data_rare
phy_data_rel

phyloseq2qiime2(phy_data)
phyloseq2qiime2(phy_data_rare)
phyloseq2qiime2(phy_data_rel)

# Run 02_b_qiime_import.sh the 02_c_qiime_unifrac.sh

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/FinalFigs_qza/phy_data-core-metrics-results")
metadata = read.csv("Combined_Metadata.csv", header = TRUE, sep = "\t")

wunifrac_cmb <- read_qza("weighted_unifrac_pcoa_results.qza")

unwunifrac_cmb <- read_qza("unweighted_unifrac_pcoa_results.qza")


# Change Experiment Names
levels(metadata$Experiment) <- c(levels(metadata$Experiment), "Field 1")
levels(metadata$Experiment) <- c(levels(metadata$Experiment), "Field 2")
levels(metadata$Experiment) <- c(levels(metadata$Experiment), "Greenhouse")

metadata$Experiment[metadata$Experiment == 'MMH'] <- ('Field 1')
metadata$Experiment[metadata$Experiment == 'END'] <- 'Field 2'
metadata$Experiment[metadata$Experiment == 'GH'] <- 'Greenhouse'

# Replace all _ with space


#  + ggtitle("Weighted Unifrac and Genetic Background") 
w1 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% filter(Inbred_or_Hybrid == "Inbred"|
                                   Inbred_or_Hybrid == "Hybrid"|
                                   Inbred_or_Hybrid == "Open_Pollinated") %>%
  ggplot( aes(x=PC1, y=PC2, color=`Inbred_or_Hybrid`)) +
  geom_point(alpha=0.5, size = 5) + scale_color_manual(values = c("Hybrid" = "firebrick",
                                                                  "Inbred" = "royalblue3",
                                                                  "Open_Pollinated" = "orange"),
                                                       labels = c("Hybrid","Inbred","Open Pollinated")) + 
  theme(legend.position="bottom", legend.text = element_text(size = 10)) + theme(axis.text = element_text(size = 20))#alpha controls transparency and helps when points are overlapping

w1$labels$colour <- "Inbred or Hybrid"
#+ ggtitle("Weighted Unifrac and Sample Type")
w2 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% filter(Sample_Type == "Rhizosphere"|
                                   Sample_Type == "Root"|
                                   Sample_Type == "Stalk") %>% ggplot( aes(x=PC1, y=PC2, color=`Sample_Type`)) +
  geom_point(alpha=0.5, size = 5) + scale_color_manual(values = c("Rhizosphere" = "purple",
                                                                  "Root" = "tan4",
                                                                  "Stalk" = "olivedrab")) + 
  theme(legend.position="bottom", legend.text = element_text(size = 10)) + theme(axis.text = element_text(size = 20)) #alpha controls transparency and helps when points are overlapping
w2$labels$colour <- "Sample Type"
#+ ggtitle("Weighted Unifrac and Experiment")

w3 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Experiment`)) +
  geom_point(alpha=0.5, size = 5) +   scale_color_manual(values = c("Field 1" = "cyan3",
                                                                    "Field 2" = "gold3",
                                                                    "Greenhouse" = "green4"
  )) + 
  theme(legend.position="bottom", legend.text = element_text(size = 10)) + theme(axis.text = element_text(size = 20)) #alpha controls transparency and helps when points are overlapping


ggcmb <- ggarrange(w3,w2,w1, nrow = 1, ncol = 3, labels = c("A","B","C"))
ggcmb

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures")
ggsave("Fig2_Beta_D2.png",
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
       ggcmb, device = "png", width = 15, height = 5, dpi = 600)

################### PERMANOVA
######## Load Data
getwd()
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/Phyloseq_Objects")

# Use main data - phy_data

# import 
# metadata <- data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ","))
# 
# otu = otu_table(as.matrix(read.table("phy_data_otu.csv", row.names = 1, header = TRUE, sep = ",")),taxa_are_rows = TRUE)
# taxa = tax_table(as.matrix(read.table("phy_data_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")))
# meta = sample_data(data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ",")), errorIfNULL = FALSE)
# tree = read.tree("phy_data_tree.tree")
# 
# # replace "." with "-" otu dfs
# colnames(otu)
# colnames(otu) <- gsub("\\.", "-", colnames(otu))
# colnames(otu)
# 
# 
# phy_data = phyloseq(otu,taxa, meta, tree)
# phy_data = merge_phyloseq(phy_data, meta, tree)
# phy_data

rbiome_wunifrac <- unifrac(otu_table(phy_data), weighted = TRUE, tree = tree)

jac <- phyloseq::distance(phy_data,"bray")

ctrl <- with(sample_data(phy_data), how(blocks = Inbred_or_Hybrid, nperm = 999))

adonis_results <- adonis2(jac ~ sample_data(phy_data)$Experiment +
                            sample_data(phy_data)$Sample_Type
                          + sample_data(phy_data)$Inbred_or_Hybrid
                          + sample_data(phy_data)$Genotype
                          + sample_data(phy_data)$Location,
                          by = "margin")
adonis_results




adonis_results <- adonis2(rbiome_wunifrac ~ sample_data(phy_data)$Experiment +
                           sample_data(phy_data)$Sample_Type
                          + sample_data(phy_data)$Inbred_or_Hybrid
                          + sample_data(phy_data)$Genotype
                          + sample_data(phy_data)$Location,
                          by = "margin")
adonis_results


# Type II ANOVA adonis
library(RVAideMemoire)
adonisII_results <- adonis.II(jac ~ sample_data(phy_data)$Experiment +
                                sample_data(phy_data)$Sample_Type
                              + sample_data(phy_data)$Inbred_or_Hybrid
                              + sample_data(phy_data)$Genotype
                              + sample_data(phy_data)$Location)

adonisII_results

