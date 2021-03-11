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

# Aglomerate Higher taxonomic levels
# Phylum = 18
# Class = 38
# Order = 89
  phyHigh <- tax_glom(phyCmbFilt, taxrank = "Class")
  phyHigh

#Before filtering out by experiments how many taxa do I have?
phyClassPreFilt <- tax_glom(phy2, taxrank = "Class")
###############
#Split by samples then remove blank OTUs and do phyloseq to qiime

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
stalksF
rhizosF
rootsF
##############################################

# Plot the differences ############################################
library("DESeq2")
packageVersion("DESeq2")
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

#stalk          !!! Every gene contains at least one zero
stalks1 = transform_sample_counts(stalksF, function(OTU) OTU + 1)
deStalk = phyloseq_to_deseq2(stalks1, ~ Inbred_or_Hybrid)
deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")

res = results(deStalk, cooksCutoff = FALSE)
alpha = 0.001
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rootsF)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
ggplot(sigtab, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ggtitle("Differential Abundacance Inbred vs Hybrid in Combined: Stalks") 


#rhizos
rhizos1 = transform_sample_counts(rhizosF, function(OTU) OTU + 1)
deRhizos = phyloseq_to_deseq2(rhizos1, ~ Inbred_or_Hybrid)
deRhizos = DESeq(deRhizos, test="Wald", fitType = "parametric")

res = results(deRhizos, cooksCutoff = FALSE)
alpha = 0.001
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab), ], "matrix"))
head(sigtab)
# There are no differentially abundant rhizosphere microbes
dim(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Class), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Differential Abundacance Inbred vs Hybrid in Combined: Rhizos")



#roots
roots1 = transform_sample_counts(rootsF, function(OTU) OTU + 1)
deRoots = phyloseq_to_deseq2(roots1, ~ Inbred_or_Hybrid)
deRoots = DESeq(deRoots, test="Wald", fitType = "parametric")

res = results(deRoots, cooksCutoff = FALSE)
alpha = 0.001
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rootsF)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Differential Abundacance Inbred vs Hybrid in Combined: Roots")

########################    Make Tables
setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/Dif_Abundance_ClassTaxa")

### Stalks
stalks1 = transform_sample_counts(stalksF, function(OTU) OTU + 1)
deStalk = phyloseq_to_deseq2(stalks1, ~ Inbred_or_Hybrid)
deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")

res <- results(deStalk, cooksCutoff = FALSE)
res

res.H.I <- results(deStalk, contrast=c("Inbred_or_Hybrid", "Hybrid", "Inbred"))
res.H.O <- results(deStalk, contrast=c("Inbred_or_Hybrid", "Hybrid", "Open_Pollinated"))
res.I.O <- results(deStalk, contrast=c("Inbred_or_Hybrid", "Inbred", "Open_Pollinated"))

res.H.I

#plotMA uses limma instead of deseq
DESeq2::plotMA(res)
DESeq2::plotMA(res.H.I)
DESeq2::plotMA(res.H.O)
DESeq2::plotMA(res.I.O)

# Diff Abundance
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(stalksF)[rownames(sigtab), ], "matrix"))
dim(sigtab)
write.csv(as.data.frame(sigtab), file="All_stalk_difresults_all.csv")

sigtab1 = res.H.I[which(res.H.I$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(stalksF)[rownames(sigtab1), ], "matrix"))
dim(sigtab1)
write.csv(as.data.frame(sigtab1), file="All_stalk_difresults_HybridvsInbred.csv")

sigtab2 = res.H.O[which(res.H.O$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(stalksF)[rownames(sigtab2), ], "matrix"))
dim(sigtab2)
write.csv(as.data.frame(sigtab2), file="All_stalk_difresults_HybridvsOpenP.csv")

sigtab3 = res.I.O[which(res.I.O$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(stalksF)[rownames(sigtab3), ], "matrix"))
dim(sigtab3)
write.csv(as.data.frame(sigtab3), file="All_stalk_difresults_InbredvsOpenP.csv")

plotDispEsts(deStalk)

###################################################
# Rhizosphere run

rhizos1 = transform_sample_counts(rhizosF, function(OTU) OTU + 1)
dStalk = phyloseq_to_deseq2(rhizos1, ~ Inbred_or_Hybrid)
dStalk = DESeq(dStalk, test="Wald")

res <- results(dStalk, cooksCutoff = FALSE)
res

res.H.I <- results(dStalk, contrast=c("Inbred_or_Hybrid", "Hybrid", "Inbred"))
res.H.O <- results(dStalk, contrast=c("Inbred_or_Hybrid", "Hybrid", "Open_Pollinated"))
res.I.O <- results(dStalk, contrast=c("Inbred_or_Hybrid", "Inbred", "Open_Pollinated"))

res.H.I

#plotMA uses limma instead of deseq
DESeq2::plotMA(res)
DESeq2::plotMA(res.H.I)
DESeq2::plotMA(res.H.O)
DESeq2::plotMA(res.I.O)

alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab), ], "matrix"))
dim(sigtab)
write.csv(as.data.frame(sigtab), file="All_rhizo_difresults_all.csv")

sigtab1 = res.H.I[which(res.H.I$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab1), ], "matrix"))
dim(sigtab1)
write.csv(as.data.frame(sigtab1), file="All_rhizo_difresults_HybridvsInbred.csv")

sigtab2 = res.H.O[which(res.H.O$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab2), ], "matrix"))
dim(sigtab2)
write.csv(as.data.frame(sigtab2), file="All_rhizo_difresults_HybridvsOpenP.csv")

sigtab3 = res.I.O[which(res.I.O$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab3), ], "matrix"))
dim(sigtab3)
write.csv(as.data.frame(sigtab3), file="All_rhizo_difresults_InbredvsOpenP.csv")


#######################################################################
# Roots run
roots1 = transform_sample_counts(rootsF, function(OTU) OTU + 1)
dStalk = phyloseq_to_deseq2(roots1, ~ Inbred_or_Hybrid)
dStalk = DESeq(dStalk, test="Wald")

res <- results(dStalk, cooksCutoff = FALSE)
res

res.H.I <- results(dStalk, contrast=c("Inbred_or_Hybrid", "Hybrid", "Inbred"))
res.H.O <- results(dStalk, contrast=c("Inbred_or_Hybrid", "Hybrid", "Open_Pollinated"))
res.I.O <- results(dStalk, contrast=c("Inbred_or_Hybrid", "Inbred", "Open_Pollinated"))

res.H.I

#plotMA uses limma instead of deseq
DESeq2::plotMA(res)
DESeq2::plotMA(res.H.I)
DESeq2::plotMA(res.H.O)
DESeq2::plotMA(res.I.O)

alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rootsF)[rownames(sigtab), ], "matrix"))
dim(sigtab)
write.csv(as.data.frame(sigtab), file="All_root_difresults_all.csv")

sigtab1 = res.H.I[which(res.H.I$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(rootsF)[rownames(sigtab1), ], "matrix"))
dim(sigtab1)
write.csv(as.data.frame(sigtab1), file="All_root_difresults_HybridvsInbred.csv")

sigtab2 = res.H.O[which(res.H.O$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(rootsF)[rownames(sigtab2), ], "matrix"))
dim(sigtab2)
write.csv(as.data.frame(sigtab2), file="All_root_difresults_HybridvsOpenP.csv")

sigtab3 = res.I.O[which(res.I.O$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(rootsF)[rownames(sigtab3), ], "matrix"))
dim(sigtab3)
write.csv(as.data.frame(sigtab3), file="All_root_difresults_InbredvsOpenP.csv")

############################################ Now look at experiments instead of Inbred and Hybrid
########################    Make Tables


# Stalks
stalks1 = transform_sample_counts(stalksF, function(OTU) OTU + 1)
deStalk = phyloseq_to_deseq2(stalks1, ~ Experiment)
deStalk = DESeq(deStalk, test="Wald", fitType = "parametric")

res <- results(deStalk, cooksCutoff = FALSE)
res

res.M.E <- results(deStalk, contrast=c("Experiment", "MMH", "END"))
res.M.G <- results(deStalk, contrast=c("Experiment", "MMH", "GH"))
res.E.G <- results(deStalk, contrast=c("Experiment", "END", "GH"))

#plotMA uses limma instead of deseq

# Diff Abundance
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(stalksF)[rownames(sigtab), ], "matrix"))
dim(sigtab)
write.csv(as.data.frame(sigtab), file="Combined_Experiments_Stalk.csv")

sigtab1 = res.M.E[which(res.M.E$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(stalksF)[rownames(sigtab1), ], "matrix"))
dim(sigtab1)
write.csv(as.data.frame(sigtab1), file="MMHvsEND_Stalk.csv")

sigtab2 = res.M.G[which(res.M.G$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(stalksF)[rownames(sigtab2), ], "matrix"))
dim(sigtab2)
write.csv(as.data.frame(sigtab2), file="MMHvsGH_Stalk.csv")

sigtab3 = res.E.G[which(res.E.G$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(stalksF)[rownames(sigtab3), ], "matrix"))
dim(sigtab3)
write.csv(as.data.frame(sigtab3), file="ENDvsGH_Stalk.csv")


###################################################
# Rhizosphere run

rhizos1 = transform_sample_counts(rhizosF, function(OTU) OTU + 1)
dStalk = phyloseq_to_deseq2(rhizos1, ~ Experiment)
deStalk = DESeq(dStalk, test="Wald")

res <- results(deStalk, cooksCutoff = FALSE)
res

res.M.E <- results(deStalk, contrast=c("Experiment", "MMH", "END"))
res.M.G <- results(deStalk, contrast=c("Experiment", "MMH", "GH"))
res.E.G <- results(deStalk, contrast=c("Experiment", "END", "GH"))

# Diff Abundance
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab), ], "matrix"))
dim(sigtab)
write.csv(as.data.frame(sigtab), file="Combined_Experiments_Rhizos.csv")

sigtab1 = res.M.E[which(res.M.E$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab1), ], "matrix"))
dim(sigtab1)
write.csv(as.data.frame(sigtab1), file="MMHvsEND_Rhizos.csv")

sigtab2 = res.M.G[which(res.M.G$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab2), ], "matrix"))
dim(sigtab2)
write.csv(as.data.frame(sigtab2), file="MMHvsGH_Rhizos.csv")

sigtab3 = res.E.G[which(res.E.G$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(rhizosF)[rownames(sigtab3), ], "matrix"))
dim(sigtab3)
write.csv(as.data.frame(sigtab3), file="ENDvsGH_Rhizos.csv")

#######################################################################
# Roots run
roots1 = transform_sample_counts(rootsF, function(OTU) OTU + 1)
dStalk = phyloseq_to_deseq2(roots1, ~ Experiment)
deStalk = DESeq(dStalk, test="Wald")

res <- results(deStalk, cooksCutoff = FALSE)
res

# No MMH Roots
#res.M.E <- results(deStalk, contrast=c("Experiment", "MMH", "END"))
#res.M.G <- results(deStalk, contrast=c("Experiment", "MMH", "GH"))
res.E.G <- results(deStalk, contrast=c("Experiment", "END", "GH"))

# Diff Abundance
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rootsF)[rownames(sigtab), ], "matrix"))
dim(sigtab)
write.csv(as.data.frame(sigtab), file="Combined_Experiments_Rootss.csv")

# sigtab1 = res.M.E[which(res.M.E$padj < alpha), ]
# sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(rootsF)[rownames(sigtab1), ], "matrix"))
# dim(sigtab1)
# write.csv(as.data.frame(sigtab1), file="MMHvsEND_Roots.csv")
# 
# sigtab2 = res.M.G[which(res.M.G$padj < alpha), ]
# sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(rootsF)[rownames(sigtab2), ], "matrix"))
# dim(sigtab2)
# write.csv(as.data.frame(sigtab2), file="MMHvsGH_Roots.csv")

sigtab3 = res.E.G[which(res.E.G$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(rootsF)[rownames(sigtab3), ], "matrix"))
dim(sigtab3)
write.csv(as.data.frame(sigtab3), file="ENDvsGH_Rootss.csv")


######################################### Start upset plot

library(HMP2Data)
library(phyloseq)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(ggplot2)
library(UpSetR)
library(tidyr)
library(data.table)

phyHigh

mergedGroups <- merge_samples(phyHigh, "Group", fun = sum)


upset_object <- as.data.frame(t(otu_table(mergedGroups)))
binary_object <- sapply(upset_object, function(x) ifelse(x > 0, 1, 0),
                        USE.NAMES = T)
rownames(binary_object) <- rownames(upset_object)
binary_object <- as.data.frame(binary_object)

upset_order <- colnames(binary_object)

noblanks = subset(binary_object[c(9:13,15:17)])

shared_ASV_plot <- upset(noblanks,
                        nsets = 9,
                        nintersects = NA,
                        keep.order = TRUE,
                        #empty.intersections = "on"
                         )
shared_ASV_plot



