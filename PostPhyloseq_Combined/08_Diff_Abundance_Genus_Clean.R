# Take pre proccessed qiime objects, filter with phyloseq, then insert back into qiime for diversity metrics!

library(tidyverse)
library(devtools)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(tidyr)

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

# Aglomerate Higher taxonomic levels - or not depending on commented out line
# Phylum = 18
# Class = 38
# Order = 89

#phyHigh <- tax_glom(phyCmbFilt, taxrank = "Species")
phyHigh <- phyCmbFilt

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
##############################################

# Plot the differences ############################################
library("DESeq2")
packageVersion("DESeq2")
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

#### Differential Abundance Functions


diff_abund <- function(phyloseq_obj, metadata_category, cat_quote, A, B, alpha_value){
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

plot_diff_abund <- function(df, title){
  plot <- ggplot(df, aes(x = Class, y = log2FoldChange, fill = Phylum)) + 
    geom_point(size=6, aes(color = Phylum)) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ggtitle(title)
  
  return(plot)
}

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/Diff_Abundance_Clean_NoAglom")

# Test function
da_df <- diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
plot_diff_abund(da_df,"test title")
write.csv(as.data.frame(da_df), file="TheTest.csv")

### Important: if there are no diff abundant taxa then error = Error in dimnames(x) <- dn : 
# length of 'dimnames' [1] not equal to array extent

# Lets Begin
#########################################################################################

# All Tissues
da_df <- diff_abund(phyHigh,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
write.csv(as.data.frame(da_df), file="AllTisues_InbredvsHybrid_ASV.csv")

da_df <- diff_abund(phyHigh,~Location,"Location", "IH", "GH", .001)
write.csv(as.data.frame(da_df), file="AllTissues_IHvsGH.csv")


# Stalks
da_df <- diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)
write.csv(as.data.frame(da_df), file="Stalks_IvsH.csv")

da_df <- diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Open_Pollinated", .001)
write.csv(as.data.frame(da_df), file="Stalks_IvsOp.csv")

da_df <- diff_abund(stalksF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Hybrid", "Open_Pollinated", .001)

da_df <- diff_abund(stalksF,~Location,"Location", "IH", "GH", .001)

# Rhizos
da_df <- diff_abund(rhizosF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(rhizosF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Open_Pollinated", .001)

da_df <- diff_abund(rhizosF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Hybrid", "Open_Pollinated", .001)

da_df <- diff_abund(rhizosF,~Location,"Location", "IH", "GH", .001)

# Roots
da_df <- diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Open_Pollinated", .001)

da_df <- diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Hybrid", "Open_Pollinated", .001)
write.csv(as.data.frame(da_df), file="Roots_HvsOp.csv")

da_df <- diff_abund(rootsF,~Location,"Location", "IH", "GH", .001)
write.csv(as.data.frame(da_df), file="Roots_IHvsGH.csv")
plot_diff_abund(da_df,"Roots Diff Abundance: IH Farm vs GH")


### By Experiment

# Stalks
stalks_GH = subset_samples(stalksF, Experiment=="GH")
stalks_END = subset_samples(stalksF, Experiment=="END")
stalks_MMH = subset_samples(stalksF, Experiment=="MMH")

da_df <- diff_abund(stalks_GH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(stalks_END,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(stalks_MMH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

# Rhizos
rhizos_GH = subset_samples(rhizosF, Experiment=="GH")
rhizos_END = subset_samples(rhizosF, Experiment=="END")
rhizos_MMH = subset_samples(rhizosF, Experiment=="MMH")

da_df <- diff_abund(rhizos_GH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(rhizos_END,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(rhizos_MMH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

# Roots
roots_GH = subset_samples(rootsF, Experiment=="GH")
roots_END = subset_samples(rootsF, Experiment=="END")

da_df <- diff_abund(roots_GH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(roots_END,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

# Compare GH and END
da_df <- diff_abund(stalksF,~Experiment,"Experiment", "GH", "END", .001)
da_df <- diff_abund(rhizosF,~Experiment,"Experiment", "GH", "END", .001)
da_df <- diff_abund(rootsF,~Experiment,"Experiment", "GH", "END", .001)



#### Compare within Experiment without filtering OTUs out that are not in every experiment
phy2

phyHigh <- tax_glom(phy2, taxrank = "Genus")
phyHigh
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

### By Experiment

# Stalks
stalks_GH = subset_samples(stalksF, Experiment=="GH")
stalks_END = subset_samples(stalksF, Experiment=="END")
stalks_MMH = subset_samples(stalksF, Experiment=="MMH")

da_df <- diff_abund(stalks_GH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(stalks_END,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(stalks_MMH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

# Rhizos
rhizos_GH = subset_samples(rhizosF, Experiment=="GH")
rhizos_END = subset_samples(rhizosF, Experiment=="END")
rhizos_MMH = subset_samples(rhizosF, Experiment=="MMH")

da_df <- diff_abund(rhizos_GH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(rhizos_END,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(rhizos_MMH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

# Roots
roots_GH = subset_samples(rootsF, Experiment=="GH")
roots_END = subset_samples(rootsF, Experiment=="END")

da_df <- diff_abund(roots_GH,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

da_df <- diff_abund(roots_END,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred", "Hybrid", .001)

# Compare GH and END
da_df <- diff_abund(stalksF,~Experiment,"Experiment", "END", "GH", .001)
da_df <- diff_abund(rhizosF,~Experiment,"Experiment", "END", "GH", .001)
da_df <- diff_abund(rootsF,~Experiment,"Experiment", "END", "GH", .001)



###Plots for paper:

da_df <- diff_abund(rootsF,~Inbred_or_Hybrid,"Inbred_or_Hybrid", "Inbred","Hybrid", .001)
dim(da_df)
plot_diff_abund(da_df,"All Roots Diff Abundance: Inbred vs Hybrid")

# Lets make a volcano plot
#BiocManager::install('EnhancedVolcano')










