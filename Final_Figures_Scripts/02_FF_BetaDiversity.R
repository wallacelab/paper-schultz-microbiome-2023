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
library(ggpubr)

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

#data sets
phyCmbFilt
stalksF
rhizosF
rootsF

##############################################

# Beta Diversity 

##############################################

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

################### Super Filtered
#Export to qiime qza
phyloseq2qiime2(phyCmbFilt)
# phyloseq2qiime2(stalksF)
# phyloseq2qiime2(rhizosF)
# phyloseq2qiime2(rootsF)

library(rbiom)
library(ape)

combo_biom <- read.biom("phyCmbFilt_features-table.biom")
combo_tree <- ape::read.tree("phyCmbFilt_tree-rooted.newick")

combo_w_unifrac = rbiom::unifrac(combo_biom, weighted = TRUE, tree = combo_tree)
combo_un_unifrac = rbiom::unifrac(combo_biom, weighted = FALSE, tree = combo_tree)

#remotes::install_github("cmmr/rbiom")
plot(combo_biom, Unifrac ~ 'Inbred_or_Hybrid', unifrac, weighted = TRUE)

library(RAM)

plot(combo_biom)

# rarify at:      need to include this for qiime alpha diversity that we wont use
sample_sums(phyCmbFilt)

#import
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
library(qiime2R)

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/FinalFigs_qza/combined-core-metrics-results")

metadata = read.csv("Combined_Metadata.csv", header = TRUE, sep = "\t")

wunifrac_cmb <- read_qza("weighted_unifrac_pcoa_results.qza")

unwunifrac_cmb <- read_qza("unweighted_unifrac_pcoa_results.qza")

# Change Experiment Names
levels(metadata$Experiment) <- c(levels(metadata$Experiment), "Field_1")
levels(metadata$Experiment) <- c(levels(metadata$Experiment), "Field_2")
levels(metadata$Experiment) <- c(levels(metadata$Experiment), "Greenhouse")

metadata$Experiment[metadata$Experiment == 'MMH'] <- ('Field_1')
metadata$Experiment[metadata$Experiment == 'END'] <- 'Field_2'
metadata$Experiment[metadata$Experiment == 'GH'] <- 'Greenhouse'

#  + ggtitle("Weighted Unifrac and Genetic Background") 
w1 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% filter(Inbred_or_Hybrid == "Inbred"|
                                   Inbred_or_Hybrid == "Hybrid"|
                                   Inbred_or_Hybrid == "Open_Pollinated") %>%
  ggplot( aes(x=PC1, y=PC2, color=`Inbred_or_Hybrid`)) +
  geom_point(alpha=0.5, size = 3) + scale_color_manual(values = c("Hybrid" = "firebrick",
                                                          "Inbred" = "royalblue3",
                                                          "Open_Pollinated" = "orange")) + 
  theme(legend.position="bottom", legend.text = element_text(size = 10)) #alpha controls transparency and helps when points are overlapping

#+ ggtitle("Weighted Unifrac and Sample Type")
w2 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% filter(Sample_Type == "Rhizosphere"|
                                   Sample_Type == "Root"|
                                   Sample_Type == "Stalk") %>% ggplot( aes(x=PC1, y=PC2, color=`Sample_Type`)) +
  geom_point(alpha=0.5, size = 3) + scale_color_manual(values = c("Rhizosphere" = "purple",
                                                          "Root" = "tan4",
                                                          "Stalk" = "olivedrab")) + 
  theme(legend.position="bottom", legend.text = element_text(size = 10)) #alpha controls transparency and helps when points are overlapping

#+ ggtitle("Weighted Unifrac and Experiment")

w3 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Experiment`)) +
  geom_point(alpha=0.5, size = 3) +   scale_color_manual(values = c("Field_1" = "cyan3",
                                                          "Field_2" = "gold3",
                                                          "Greenhouse" = "green4"
                                                          )) + 
  theme(legend.position="bottom", legend.text = element_text(size = 10)) #alpha controls transparency and helps when points are overlapping


ggcmb <- ggarrange(w3,w2,w1, nrow = 1, ncol = 3, labels = c("A","B","C"))
ggcmb

ggsave("Fig2_Beta.png", 
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/PaperFigures",
       ggcmb, device = "png", width = 15, height = 5, dpi = 600)









u1 <- unwunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Inbred_or_Hybrid`)) +
  geom_point(alpha=0.5) + ggtitle("Unweighted Unifrac and Genetic Background") #alpha controls transparency and helps when points are overlapping

u2 <- unwunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Sample_Type`)) +
  geom_point(alpha=0.5) + ggtitle("Unweighted Unifrac and Sample Type") 

u3 <- unwunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Experiment`)) +
  geom_point(alpha=0.5) + ggtitle("Unweighted Unifrac and Experiment") 



#######################################################################################################################################################




ggarrange(u1,u2,u3, nrow = 1, ncol = 3, labels = c("A","B","C"))

######################## Less Filtered - Use more filtered Weighted for paper
phy2
#Export to qiime qza
phyloseq2qiime2(phy2)
# phyloseq2qiime2(stalksF)
# phyloseq2qiime2(rhizosF)
# phyloseq2qiime2(rootsF)

library(rbiom)
library(ape)

combo_biom <- read.biom("phy2_features-table.biom")
combo_tree <- ape::read.tree("phy2_tree-rooted.newick")

# Now run 02_b and 02_c .sh scripts

combo_w_unifrac = rbiom::unifrac(combo_biom, weighted = TRUE, tree = combo_tree)
combo_un_unifrac = rbiom::unifrac(combo_biom, weighted = FALSE, tree = combo_tree)

#remotes::install_github("cmmr/rbiom")
# plot(combo_biom, Unifrac ~ 'Inbred_or_Hybrid', unifrac, weighted = TRUE)

library(RAM)

# plot(combo_biom)

# rarify at:      need to include this for qiime alpha diversity that we wont use
sample_sums(phy2)

#import
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
library(qiime2R)

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_qza_files/FinalFigs_qza/phy2_lessfiltered-core-metrics-results")

metadata = read.csv("Combined_Metadata.csv", header = TRUE, sep = "\t")

wunifrac_cmb <- read_qza("weighted_unifrac_pcoa_results.qza")

unwunifrac_cmb <- read_qza("unweighted_unifrac_pcoa_results.qza")


w1 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Inbred_or_Hybrid`)) +
  geom_point(alpha=0.5) + ggtitle("Weighted Unifrac and Genetic Background") #alpha controls transparency and helps when points are overlapping

w2 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Sample_Type`)) +
  geom_point(alpha=0.5) + ggtitle("Weighted Unifrac and Sample Type") #alpha controls transparency and helps when points are overlapping

w3 <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Experiment`)) +
  geom_point(alpha=0.5) + ggtitle("Weighted Unifrac and Experiment") #alpha controls transparency and helps when points are overlapping

u1 <- unwunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Inbred_or_Hybrid`)) +
  geom_point(alpha=0.5) + ggtitle("Unweighted Unifrac and Genetic Background") #alpha controls transparency and helps when points are overlapping

u2 <- unwunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Sample_Type`)) +
  geom_point(alpha=0.5) + ggtitle("Unweighted Unifrac and Sample Type") 

u3 <- unwunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Experiment`)) +
  geom_point(alpha=0.5) + ggtitle("Unweighted Unifrac and Experiment") 


ggarrange(w1,w2,w3, nrow = 1, ncol = 3, labels = c("A","B","C"))

ggarrange(u1,u2,u3, nrow = 1, ncol = 3, labels = c("A","B","C"))

# Figures for PAG

w1_c <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Inbred_or_Hybrid`)) +
  geom_point(alpha=0.5) + theme(legend.position = "none")  #alpha controls transparency and helps when points are overlapping

w2_c <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Sample_Type`)) +
  geom_point(alpha=0.5) + theme(legend.position = "none")  #alpha controls transparency and helps when points are overlapping

w3_c <- wunifrac_cmb$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>% ggplot( aes(x=PC1, y=PC2, color=`Experiment`)) +
  geom_point(alpha=0.5) + theme(legend.position = "none") #alpha controls transparency and helps when points are overlapping

ggsave("PAG_beta_background", w1_c, device = "png", width = 5, height = 5, dpi = 300)

ggsave("PAG_beta_tissue", w2_c, device = "png", width = 5, height = 5, dpi = 300)

ggsave("PAG_beta_location", w3_c, device = "png", width = 5, height = 5, dpi = 300)
