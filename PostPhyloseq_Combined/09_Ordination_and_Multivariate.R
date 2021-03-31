library(tidyverse)
library(devtools)
devtools::install_github("jbisanz/qiime2R")
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

# Stacked bar plot

combined_phylum <- phyCmbFilt %>%                             
  tax_glom(taxrank = "Phylum") %>%                            # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%        # Transform to rel. abundance     
  psmelt() %>%                                                # Melt to long format
  arrange(Phylum)                                             # Sort data frame alphabetically by phylum

# Plot it
#barplot_taxa <- ggplot(combined_phylum, aes(x = Experiment, y = Abundance, fill = Phylum)) +
#  geom_bar(stat = "identity")

barplot_taxa <- plot_bar(phyCmbFilt, x= "Sample_Type", fill = "Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# This needs more work but moving on - Can always plot by tissue with Experiment on x Axis - That would look nice

########################################### https://micca.readthedocs.io/en/latest/phyloseq.html
# Alpha Diversity - not rarefied - can just rarefy and plug that phyloseq object in for CmbFlt

alpha_plot <- plot_richness(phyCmbFilt, x="Experiment", measures=c("Observed", "Shannon", "ACE"), 
                            color = "Inbred_or_Hybrid", 
                            shape = "Sample_Type")


# Unconstrained Ordinations - use tissue specific phyloseq objects for actual pictures and labeling

#Unweighted UniFrac
wunifrac_dist = phyloseq::distance(phyCmbFilt, method = "unifrac", weighted=F)
wuni_ordination = ordinate(phyCmbFilt, method = "PCoA", distance = wunifrac_dist )
plot_ordination(phyCmbFilt, wuni_ordination, color = "Experiment", shape = "Inbred_or_Hybrid")

#Bray Curtis
bc_ordination = ordinate(phyCmbFilt, method = "PCoA", distance = "bray")
plot_ordination(phyCmbFilt, bc_ordination, color = "Experiment", shape = "Inbred_or_Hybrid")

#nMDS - Why is one stalk sample all the way over there lol?
set.seed(18)

# MMH-2GH-42K has an nmds1 score of 44000
combo_nmds <- ordinate(physeq = phyCmbFilt, method = "NMDS", distance = "bray")

plot_ordination(physeq = phyCmbFilt, ordination = combo_nmds,
                color = "Experiment", shape = "Sample_Type")

# Remove MMH-2GH-42K because it is super weird and needs more investigation
combined_new = subset_samples(phyCmbFilt, sample_names(phyCmbFilt) != "MMH-2GH-42K")
combo_nmds2 <- ordinate(physeq = combined_new, method = "NMDS", distance = "bray")
data.scores = as.data.frame(scores(combo_nmds2))
plot_ordination(physeq = combined_new, ordination = combo_nmds2,
                color = "Experiment", shape = "Sample_Type")
# Need to look back at fastqC and see if there is something wrong with this sample???????? - The OTU table is all zeros
rhizos_nmds <- ordinate(physeq = rhizosF, method = "NMDS", distance = "bray")

plot_ordination(physeq = rhizosF, ordination = rhizos_nmds,
                color = "Experiment", shape = "Inbred_or_Hybrid")

weird = subset_samples(phyCmbFilt, sample_names(phyCmbFilt) == "MMH-2GH-42K")
###### Notes about this sample - low %Duplication (most are above 80%) Has like no reads, its OTU table is empty

stalk_nmds <- ordinate(physeq = stalksF, method = "NMDS", distance = "bray")

plot_ordination(physeq = stalksF, ordination = stalk_nmds,
                color = "Experiment", shape = "Inbred_or_Hybrid")

# PERMANOVA
adonist_results <- adonis(wunifrac_dist ~ sample_data(phyCmbFilt)$Experiment 
                          + sample_data(phyCmbFilt)$Sample_Type
                          + sample_data(phyCmbFilt)$Inbred_or_Hybrid
                          + sample_data(phyCmbFilt)$Genotype
                          + sample_data(phyCmbFilt)$Location)
adonist_results
 
# Why do we use sample_data call from phyloseq?   https://micca.readthedocs.io/en/latest/phyloseq.html same as http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#unconstrained_ordinations
# That acesses our data in the phyloseq class

BC.dist = phyloseq::distance(phyCmbFilt, method = "bray")
sampledf = data.frame(sample_data(phyCmbFilt))

bc_adonis <- adonis(BC.dist ~ Sample_Type + 
                      Inbred_or_Hybrid + Inbred_or_Hybrid*Sample_Type + 
                      Inbred_or_Hybrid*Genotype +Genotype, data = sampledf, strata = sampledf$Experiment)

# you can do nesting like A:B:C so I could have A + B + A:B + A:B:C
# A question for doctor wallace

# Constrained Ordination
cap_ord <- ordinate(
  physeq = phyCmbFilt,
  method = "CAP",
  distance = BC.dist,
  formula = BC.dist ~ Experiment + Sample_Type + Inbred_or_Hybrid + Genotype + Location + Block
  #BC.dist ~ Experiment + Sample_Type + Inbred_or_Hybrid + Genotype + Location + Block
)

cap_plot <- plot_ordination(
  physeq = phyCmbFilt,
  ordination = cap_ord,
  color = "Experiment",
  shape = "Sample_Type"
)

arrowmat <- vegan::scores(cap_ord, display = "bp")

#arrowmat <- arrowmat[0:9]
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Only keep things that might be helpful for us
arrowdf <- arrowdf[0:9,]

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

# Next step: break into tissues and do all these plots!




