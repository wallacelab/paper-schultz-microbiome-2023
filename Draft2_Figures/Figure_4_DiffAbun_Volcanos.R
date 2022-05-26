# Volcano Plots
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
library(DESeq2)
library(ggrepel)


setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS")

# Read in the data

metadata <- read.csv("Combined_Key.csv", header = TRUE, sep = "\t")
metadata$Group <- paste(metadata$Sample_Type_Blanks_differentiated, metadata$Experiment, sep = "_")

newmeta <- read.csv("Combined_Key_withGroup.tsv", header = TRUE, sep = "\t")
head(newmeta)

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

# Class, Genus, or ASV level
phyClass <- tax_glom(phyCmbFilt, taxrank = "Class")
phyGenus <- tax_glom(phyCmbFilt, taxrank = "Genus")
phyHigh <- phyCmbFilt

# Use ASV level for Volcano Plot
phyHigh

#Before filtering out by experiments how many taxa do I have?
phyClassPreFilt <- tax_glom(phy2, taxrank = "Class")
###############
stalks <- subset_samples(phyHigh, Sample_Type_Blanks_differentiated=="Stalk")
rhizos <- subset_samples(phyHigh, Sample_Type_Blanks_differentiated=="Rhizosphere")
roots <- subset_samples(phyHigh, Sample_Type_Blanks_differentiated=="Root")

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
phyHigh
stalksF
rhizosF
rootsF
#### Differential Abundance Functions

diff_abund <- function(phyloseq_obj, metadata_category, cat_quote, A, B, alpha_value,color){
  compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
  phydesq = phyloseq_to_deseq2(compartment_t, design = metadata_category)
  deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
  
  res = results(deseq_obj, cooksCutoff = FALSE, contrast = c(cat_quote, A, B))
  alpha = alpha_value
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
  print(dim(sigtab))
  return(sigtab)
}

plot_diff_abund <- function(phyloseq_obj, metadata_category, cat_quote, A, B, alpha_value){
    compartment_t = transform_sample_counts(phyloseq_obj, function(OTU) OTU +1)
    phydesq = phyloseq_to_deseq2(compartment_t, design = metadata_category)
    deseq_obj = DESeq(phydesq, test = "Wald", fitType = "parametric")
    res = results(deseq_obj, cooksCutoff = FALSE, contrast = c(cat_quote, A, B))
    alpha = alpha_value
    res = cbind(as(res, "data.frame"), as(tax_table(phyloseq_obj)[rownames(res), ], "matrix"))
    #plot <- ggplot(res, aes(x = Class, y = log2FoldChange, fill = which(res$padj < alpha))) + 
    #geom_point(size=6, aes(color = the_color)) + 
    # theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ggtitle(title)
  
  return(res)
}
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/Diff_Abun_ASVs_D2")

### Important: if there are no diff abundant taxa then error = Error in dimnames(x) <- dn : 
# length of 'dimnames' [1] not equal to array extent

# Lets Begin - Inbred vs Hybrid
#########################################################################################
### All tissues
All_Tis_df <- diff_abund(phyHigh,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
write.csv(as.data.frame(All_Tis_df), file="AllTisues_InbredvsHybrid_ASV.csv")
All_Tis_IO <- diff_abund(phyHigh,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Open_Pollinated", .001)
write.csv(as.data.frame(All_Tis_IO), file="AllTisues_InbredvsOP_ASV.csv")
All_Tis_HO <- diff_abund(phyHigh,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Hybrid", "Open_Pollinated", .001)
write.csv(as.data.frame(All_Tis_HO), file="AllTisues_HybridvsOP_ASV.csv")
# Stalks
Stalk_df <- diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
write.csv(as.data.frame(Stalk_df), file="Stalks_IvsH_ASV.csv")
Stalk_IO <- diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Open_Pollinated", .001)
write.csv(as.data.frame(Stalk_IO), file="Stalks_IvsO_ASV.csv")
Stalk_HO <- diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Hybrid", "Open_Pollinated", .001)
write.csv(as.data.frame(Stalk_HO), file="Stalks_HvsO_ASV.csv")
# Rhizos
Rhiz_df <- diff_abund(rhizosF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
write.csv(as.data.frame(Rhiz_df), file="Rhizos_IvsH_ASV.csv")
Rhiz_IO <- diff_abund(rhizosF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Open_Pollinated", .001)
write.csv(as.data.frame(Rhiz_IO), file="Rhizos_IvsO_ASV.csv")
Rhiz_HO <- diff_abund(rhizosF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Hybrid", "Open_Pollinated", .001)
write.csv(as.data.frame(Rhiz_df), file="Rhizos_HvsO_ASV.csv")
# Roots
Root_df <- diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
write.csv(as.data.frame(Root_df), file="Roots_IvsH_ASV.csv")
Root_IO <- diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Open_Pollinated", .001)
write.csv(as.data.frame(Root_IO), file="Roots_IvsO_ASV.csv")
Root_HO <- diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Hybrid", "Open_Pollinated", .001)
write.csv(as.data.frame(Root_HO), file="Roots_HvsO_ASV.csv")
### Compare Tissues No I vs H
stalk_vs_root <- diff_abund(phyHigh,~Sample_Type_Blanks_differentiated,
                            "Sample_Type_Blanks_differentiated", "Stalk", "Root", .001)
write.csv(as.data.frame(stalk_vs_root), file="StalksvsRoots.csv")

stalk_vs_rhizo <- diff_abund(phyHigh,~Sample_Type_Blanks_differentiated,
                            "Sample_Type_Blanks_differentiated", "Stalk", "Rhizosphere", .001)
write.csv(as.data.frame(stalk_vs_rhizo), file="StalksvsRhizos.csv")

root_vs_rhizo <- diff_abund(phyHigh,~Sample_Type_Blanks_differentiated,
                             "Sample_Type_Blanks_differentiated", "Root", "Rhizosphere", .001)
write.csv(as.data.frame(root_vs_rhizo), file="RootvsRhizos.csv")
#############################################################################
# plot inbred vs hybrid into separate graphs
#Stalks
stalk_df_all<- plot_diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
plot_df <- stalk_df_all


plot_df$Significance <- ifelse(plot_df$padj > 0.001 | is.na(plot_df$padj), 'NotSig', ifelse(plot_df$log2FoldChange > 0,'Hybrid','Inbred') )
sih <- ggplot(data=plot_df, aes(x=log2FoldChange, y=log10(pvalue))) + 
  geom_point(aes(colour = Significance), size = 6) + scale_y_reverse(limits = c(0, -15)) + scale_x_continuous(limits = c(-4, 4)) + 
  scale_color_manual(values = c("Inbred" = "royalblue3",
                                "Hybrid" = "firebrick",
                                "NotSig" = "grey"))+ theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) + geom_vline(xintercept = 0, colour = "blue") + 
  geom_hline(yintercept = 0, colour = "blue") + ggtitle("Stalks") +
  theme(legend.position="none", axis.title.x = element_blank()) + theme(plot.title = element_text(colour = "olivedrab", size = 16, face = "bold")) +
  geom_text_repel(
  data = subset(plot_df, padj < 0.001),
  aes(label = Family),
  size = 5)

#Roots
root_df_all<- plot_diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
plot_df <- root_df_all

plot_df$Significance <- ifelse(plot_df$padj > 0.001 | is.na(plot_df$padj), 'NotSig', ifelse(plot_df$log2FoldChange > 0,'Hybrid','Inbred') )
rootih <- ggplot(data=plot_df, aes(x=log2FoldChange, y=log10(pvalue))) + 
  geom_point(aes(colour = Significance), size = 6) + scale_y_reverse(limits = c(0, -15)) + scale_x_continuous(limits = c(-4, 4)) +
  scale_color_manual(values = c("Inbred" = "royalblue3",
                                "Hybrid" = "firebrick",
                                "NotSig" = "grey"))+ theme(axis.text=element_text(size=18),
                                                           axis.title=element_text(size=18),legend.position="none") + geom_vline(xintercept = 0, colour = "blue") + 
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, colour = "blue") + ggtitle("Roots") + theme(plot.title = element_text(colour = "tan3", size = 16, face = "bold")) +
  geom_text_repel(
    data = subset(plot_df, padj < 0.001),
    aes(label = Family),
    size = 5)

#Rhizos
rhizos_df_all<- plot_diff_abund(rhizosF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
plot_df <- rhizos_df_all

plot_df$Significance <- ifelse(plot_df$padj > 0.001 | is.na(plot_df$padj), 'NotSig', ifelse(plot_df$log2FoldChange > 0,'Hybrid','Inbred') )
rhizih <- ggplot(data=plot_df, aes(x=log2FoldChange, y=log10(pvalue))) + 
  geom_point(aes(colour = Significance), size = 6) + scale_y_reverse(limits = c(0, -15)) + scale_x_continuous(limits = c(-4, 4)) +
  scale_color_manual(values = c("Inbred" = "royalblue3",
                                "Hybrid" = "firebrick",
                                "NotSig" = "grey"),
                     labels = c("Inbred","Hybrid", "Not Significant"))+ theme(axis.text=element_text(size=18),
                                                           axis.title=element_text(size=18)) + geom_vline(xintercept = 0, colour = "blue") + 
  geom_hline(yintercept = 0, colour = "blue") + ggtitle("Rhizos") + theme(plot.title = element_text(colour = "purple", size = 16, face = "bold")) + theme(legend.position="bottom",
                                                                          legend.title = element_text(size = 16),
                                                                          legend.text = element_text(size = 16)) +
  geom_text_repel(
    data = subset(plot_df, padj < 0.001),
    aes(label = Family),
    size = 5)
rhizih$labels$colour <- 'More Abundant In:'

IvH_Fig <- ggarrange(sih,rootih,rhizih, nrow = 3, ncol = 1)
IvH_Fig

ggsave("Fig4_IvsH_D2.png",
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
       IvH_Fig, device = "png", width = 10, height = 12, dpi = 600)
################################################
################################################
### Tissue Comparison

#Stalk vs Root
svroot_all<- plot_diff_abund(phyHigh,~Sample_Type_Blanks_differentiated,
                             "Sample_Type_Blanks_differentiated", "Stalk", "Root", .001)
plot_df <- svroot_all

plot_df$Significance <- ifelse(plot_df$padj > 0.001, 'NotSig', ifelse(plot_df$log2FoldChange > 0,'Stalk','Root') )
svroot_fig <- ggplot(data=plot_df, aes(x=log2FoldChange, y=log10(pvalue))) + 
  geom_point(aes(colour = Significance), size = 6) + coord_flip() + scale_y_reverse() +
  scale_color_manual(values = c("Stalk" = "olivedrab",
                                "Root" = "tan3",
                                "NotSig" = "grey"))+ theme(axis.text=element_text(size=18),
                                                           axis.title=element_text(size=18)) + geom_vline(xintercept = 0, colour = "blue") + 
  geom_hline(yintercept = 0, colour = "blue") + ggtitle("Stalk vs Roots")

# Stalk vs Rhizos
svrhizo_all<- plot_diff_abund(phyHigh,~Sample_Type_Blanks_differentiated,
                             "Sample_Type_Blanks_differentiated", "Stalk", "Rhizosphere", .001)
plot_df <- svrhizo_all

plot_df$Significance <- ifelse(plot_df$padj > 0.001, 'NotSig', ifelse(plot_df$log2FoldChange > 0,'Stalk','Rhizo') )
svrhiz_fig <- ggplot(data=plot_df, aes(x=log2FoldChange, y=log10(pvalue))) + 
  geom_point(aes(colour = Significance), size = 6) + coord_flip() + scale_y_reverse() +
  scale_color_manual(values = c("Stalk" = "olivedrab",
                                "Rhizo" = "purple",
                                "NotSig" = "grey"))+ theme(axis.text=element_text(size=18),
                                                           axis.title=element_text(size=18)) + geom_vline(xintercept = 0, colour = "blue") + 
  geom_hline(yintercept = 0, colour = "blue") + ggtitle("Stalk vs Rhizos")

# Roots vs Rhizos
rvrhizo_all<- plot_diff_abund(phyHigh,~Sample_Type_Blanks_differentiated,
                              "Sample_Type_Blanks_differentiated", "Root", "Rhizosphere", .001)
plot_df <- rvrhizo_all

plot_df$Significance <- ifelse(plot_df$padj > 0.001, 'NotSig', ifelse(plot_df$log2FoldChange > 0,'Root','Rhizo') )
rvrhize_fig <- ggplot(data=plot_df, aes(x=log2FoldChange, y=log10(pvalue))) + 
  geom_point(aes(colour = Significance), size = 6) + coord_flip() + scale_y_reverse() +
  scale_color_manual(values = c("Root" = "tan3",
                                "Rhizo" = "purple",
                                "NotSig" = "grey"))+ theme(axis.text=element_text(size=18),
                                                           axis.title=element_text(size=18)) + geom_vline(xintercept = 0, colour = "blue") + 
  geom_hline(yintercept = 0, colour = "blue") + ggtitle("Root vs Rhizos")

Tissue_Fig <- ggarrange(svroot_fig,svrhiz_fig,rvrhize_fig, nrow = 3, ncol = 1)

ggsave("Fig4_Tissue_D2.png",
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
       Tissue_Fig, device = "png", width = 10, height = 8, dpi = 600)















