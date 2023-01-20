# Alpha Diversity Figure
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


getwd()
setwd("D:/Manual_Backup/May_2022/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/Phyloseq_Objects")

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

head(otu_table(phy_data_rel))
##############################################

# Alpha Diversity 

##############################################
#All in one and hill numbers - q2 Functional diversity

phy_data_rare <- subset_samples(phy_data_rare, Sample_Type!="Soil")

fig1 <- plot_richness(phy_data_rare, x="Sample_Type",
                           measures=c("Observed", "Shannon"), 
              color = "Inbred_or_Hybrid",
              shape = "Experiment", 
               nrow = 3) + 
  geom_point(size=5, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("Hybrid" = "firebrick",
                                "Inbred" = "royalblue3",
                                "Open_Pollinated" = "orange"),
                     labels = c("Hybrid","Inbred","Open Pollinated")) + 
  theme(strip.text = element_text(size = 12)) + scale_shape_manual(values=c(1,2,3),
                                                                   breaks = c("MMH","END","GH"),
                                                                   labels = c("Field 1", "Field 2", "Greenhouse"),
                                                                   guide = guide_legend(reverse = TRUE)) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.text.y  = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(angle = 0,face = "bold",size = 12, hjust = .5)) + xlab("Sample Type") + 
  labs(color = "Sample Type") + guides(color = guide_legend(order = 1),
                                       shape = guide_legend(order = 2)) 



alpha_df <- fig1$data
# # Subset and append to df for 
# 
q1_df <- alpha_df[grep("Shannon", alpha_df$variable), ]
# q2_df <- alpha_df[grep("Simpson", alpha_df$variable), ]
# 
q1_df$value <- ifelse(q1_df$variable == "Shannon",
                         (exp(q1_df$value)), q1_df$value)
# q2_df$value <- ifelse(q2_df$variable == "Simpson",
#                          (1/(q2_df$value)), q2_df$value)
# 
q1_df$variable <- gsub("Shannon", "Hills q1", q1_df$variable)
# q2_df$variable <- gsub("Simpson", "Hills q2", q2_df$variable)
# 
q1_df
# q2_df
# 
alpha_df <- rbind(alpha_df,q1_df)
# alpha_df <- rbind(alpha_df,q2_df)
alpha_df
# 
fig1$data <- alpha_df

fig1$layers <- fig1$layers[-1]
fig1$data %>% mutate(Experiment = fct_relevel(Experiment, 
                                                           "MMH","END","GH"))
fig1
fig1 +facet_grid(rows = vars(Sample_Type))
fig1


### Create new figs based on alpha_df
alpha_df$Experiment <- factor(alpha_df$Experiment, levels = c("MMH","END","GH"))

obsfig <- ggplot(alpha_df %>% filter(variable == "Observed"), aes(x = Experiment,y = value, shape = Experiment, color = Inbred_or_Hybrid)) +
  geom_point(size=5, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("Hybrid" = "firebrick",
                                "Inbred" = "royalblue3",
                                "Open_Pollinated" = "orange"),
                     guide = "none",
                     labels = c("Hybrid","Inbred","Open Pollinated")) + 
  theme(strip.text = element_text(size = 12)) + scale_shape_manual(values=c(1,2,3),
                                                                   breaks = c("MMH","END","GH"),
                                                                   labels = c("Field 1", "Field 2", "Greenhouse"),
                                                                   guide = "none") + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y  = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  labs(color = "Sample Type") + guides(color = guide_legend(order = 1),
                                       shape = guide_legend(order = 2)) + facet_grid(cols = vars(Sample_Type)) +
  ggtitle("Observed")

shannonfig <- ggplot(alpha_df %>% filter(variable == "Shannon"), aes(x = Experiment,y = value, shape = Experiment, color = Inbred_or_Hybrid)) +
  geom_point(size=5, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("Hybrid" = "firebrick",
                                "Inbred" = "royalblue3",
                                "Open_Pollinated" = "orange"),
                     labels = c("Hybrid","Inbred","Open Pollinated")) + 
  theme(strip.text = element_text(size = 12)) + scale_shape_manual(values=c(1,2,3),
                                                                   breaks = c("MMH","END","GH"),
                                                                   labels = c("Field 1", "Field 2", "Greenhouse"),
                                                                   guide = guide_legend(reverse = TRUE)) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(axis.text.y  = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + ylab("Alpha Diversity Values") + 
  labs(color = "Sample Type") + guides(color = guide_legend(order = 1),
                                       shape = guide_legend(order = 2)) + facet_grid(cols = vars(Sample_Type)) +
  ggtitle("Shannon")

q1fig <- ggplot(alpha_df %>% filter(variable == "Hills q1"), aes(x = Experiment,y = value, shape = Experiment, color = Inbred_or_Hybrid)) +
  geom_point(size=5, position = position_dodge(width = 1)) +
  scale_color_manual(values = c("Hybrid" = "firebrick",
                                "Inbred" = "royalblue3",
                                "Open_Pollinated" = "orange"),
                     guide = "none",
                     labels = c("Hybrid","Inbred","Open Pollinated")) + 
  theme(strip.text = element_text(size = 12)) + scale_shape_manual(values=c(1,2,3),
                                                                   breaks = c("MMH","END","GH"),
                                                                   labels = c("Field 1", "Field 2", "Greenhouse"),
                                                                   guide = "none") + 
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y  = element_text(size = 12, face="bold")) +
  theme(axis.text.x = element_text(angle = 45,face = "bold",size = 12, hjust = .5, vjust = .5)) + 
  xlab("Sample Type") +
  labs(color = "Sample Type") + guides(color = guide_legend(order = 1),
                                       shape = guide_legend(order = 2)) + facet_grid(cols = vars(Sample_Type)) +
  ggtitle("Hills q1") + scale_x_discrete(breaks = c("MMH","END","GH"),
                                         labels = c("Field 1", "Field 2", "Greenhouse"))

library(grid)
#stacked_fig <- ggarrange(obsfig,shannonfig,q1fig, nrow = 3)

gt <- ggplot_gtable(ggplot_build(q1fig))
gt$layout$clip[grep("panel-2-\\d+", q1fig$layout$name)] <- "off"
grid.draw(gt)

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 3, position = c("right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

stacked_fig <- grid_arrange_shared_legend(obsfig,shannonfig,q1fig, ncol = 1, nrow = 3)

#annotate_figure(stacked_fig, bottom = textGrob("Field 1 Field 2 Greenhouse", rot = 45, gp = gpar(cex = 1.3)))
# Try an easier way 
library(cowplot)
thelegend <- get_legend(shannonfig)

rm_legend <- function(p){p + theme(legend.position = "none")}

#plots <- ggarrange(rm_legend(obsfig),rm_legend(shannonfig),rm_legend(q1fig), ncol = 1, nrow = 3, align = "v")

library(patchwork)
plots2 <- rm_legend(obsfig) + rm_legend(shannonfig) + rm_legend(q1fig) + plot_layout(ncol = 1, heights = c(1,1,1))

stacked_fig <- ggarrange(plots2, thelegend, widths = c(1, .15))


ggsave("Fig1_Alpha_D3_rare.png", plot = stacked_fig, device = "png", width = 9, height = 9, 
       units = c("in"), dpi = 750, path = "D:/Manual_Backup/May_2022/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/Draft3_Figures")

################################################################ 9 Panel Figure
# DATA = phy_data_rel
# 
# plot_richness(DATA, x="Experiment", measures=c("Observed", "Shannon", "Simpson"), 
#               color = "Inbred_or_Hybrid",
#               shape = "Sample_Type", title = "Alpha Diversity: Combined Experiments") + geom_point(alpha = .05)
# 
# stalks <- subset_samples(DATA, Sample_Type_Blanks_differentiated=="Stalk")
# rhizos <- subset_samples(DATA, Sample_Type_Blanks_differentiated=="Rhizosphere")
# roots <- subset_samples(DATA, Sample_Type_Blanks_differentiated=="Root")
# 
# stalk_alpha <- plot_richness(stalks, x="Experiment", measures=c("Observed", "Shannon", "Simpson"), 
#                              color = "Inbred_or_Hybrid", 
#                              title = "Alpha Diversity: Combined Experiments - Stalks") + geom_point(position = position_dodge(width = .5))
# stalk_alpha$layers <- stalk_alpha$layers[-1]
# #stalk_alpha
# 
# rhizos_alpha <- plot_richness(rhizos, "Experiment", measures=c("Observed", "Shannon", "Simpson"), 
#                               color = "Inbred_or_Hybrid",
#                               title = "Alpha Diversity: Combined Experiments - Rhizos") + geom_point(position = position_dodge(width = .5))
# rhizos_alpha$layers <- rhizos_alpha$layers[-1]
# #rhizos_alpha
# 
# roots_alpha <- plot_richness(roots, x="Experiment", measures=c("Observed", "Shannon", "Simpson"), 
#                              color = "Inbred_or_Hybrid", 
#                              title = "Alpha Diversity: Combined Experiments - Roots") + geom_point(position = position_dodge(width = .5))
# roots_alpha$layers <- roots_alpha$layers[-1]
# #roots_alpha
# 
# 
# alpha_diversity <- ggarrange(stalk_alpha, roots_alpha, rhizos_alpha, ncol = 1, nrow = 3, labels = c("A","B","C"))
# alpha_diversity

###########################################################################
# Compare with Parametric and Nonparametric tools
############################################################################
metadata <- data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ","))
metadata <- tibble::rownames_to_column(metadata, "SampleID")
phy_data_rel

### ALL SAMPLES
# build sample alpha values plus metadata
library(data.table)
# Write a function for alpha_table with metadata
alpha_tabulate <- function(phyloseq_object, metadata_file){
  all_rich <- estimate_richness(phyloseq_object, measures = c("Observed", "Shannon", "Simpson"))
  #covert rownames to columns
  all_rich <- cbind(rownames(all_rich), data.frame(all_rich, row.names = NULL))
  colnames(all_rich)[1] <- "SampleID"
  all_rich
  # replace . with -
  all_rich$SampleID <- gsub('.', '-', all_rich$SampleID, fixed = TRUE)
  all_rich
  all_rich_sort <- all_rich[order(match(all_rich$SampleID,metadata_file$SampleID)), ]
  all_rich_sort
  # Now that it's sorted cbind meta data labels for ANOVA. 
  big_alpha <- data.frame(data.table::setDT(all_rich_sort)[setDT(metadata_file), on = c("SampleID")])
  # drop blanks
  big_alpha <- big_alpha[!grepl("Blank", big_alpha$Inbred_or_Hybrid),]
  # drop soil
  big_alpha <- big_alpha[!grepl("Soil", big_alpha$Inbred_or_Hybrid),]
}

big_alpha <- alpha_tabulate(phy_data_rel,metadata)

big_a <- big_alpha

big_a <- mutate(big_a, Hillq1 = exp(Shannon))


#################################################################
## Use Non Parametric Tests
stalks_alpha <- big_a %>% filter(Sample_Type == "Stalk")
roots_alpha <- big_a %>% filter(Sample_Type == "Root")
rhizos_alpha <- big_a %>% filter(Sample_Type == "Rhizosphere")

#### Use the pairwise wilcox test and dunnTest - Non Parametric multiple comparison
#install.packages("FSA")
library("FSA")

### All Samples
# Observed
dunnTest(Observed ~ Inbred_or_Hybrid, data = big_a, 
         method = "bh")
dunnTest(Observed ~ Experiment, data = big_a, 
         method = "bh")
dunnTest(Observed ~ Sample_Type, data = big_a, 
         method = "bh")
# Shannon
dunnTest(Shannon ~ Inbred_or_Hybrid, data = big_a, 
         method = "bh")
dunnTest(Shannon ~ Experiment, data = big_a, 
         method = "bh")
dunnTest(Shannon ~ Sample_Type, data = big_a, 
         method = "bh")
# Hills q1
dunnTest(Hillq1 ~ Inbred_or_Hybrid, data = big_a, 
         method = "bh")
dunnTest(Hillq1 ~ Experiment, data = big_a, 
         method = "bh")
dunnTest(Hillq1 ~ Sample_Type, data = big_a, 
         method = "bh")


########################## Stalks
# Observed
dunnTest(Observed ~ Inbred_or_Hybrid, data = stalks_alpha, 
         method = "bh")
dunnTest(Observed ~ Experiment, data = stalks_alpha, 
         method = "bh")
# Shannon
dunnTest(Shannon ~ Inbred_or_Hybrid, data = stalks_alpha, 
         method = "bh")
dunnTest(Shannon ~ Experiment, data = stalks_alpha, 
         method = "bh")
# Hills q1
dunnTest(Hillq1 ~ Inbred_or_Hybrid, data = stalks_alpha, 
         method = "bh")
dunnTest(Hillq1 ~ Experiment, data = stalks_alpha, 
         method = "bh")

########################## Roots
# Observed
dunnTest(Observed ~ Inbred_or_Hybrid, data = roots_alpha, 
         method = "bh")
dunnTest(Observed ~ Experiment, data = roots_alpha, 
         method = "bh")
# Shannon
dunnTest(Shannon ~ Inbred_or_Hybrid, data = roots_alpha, 
         method = "bh")
dunnTest(Shannon ~ Experiment, data = roots_alpha, 
         method = "bh")
# Hills q1
dunnTest(Hillq1 ~ Inbred_or_Hybrid, data = roots_alpha, 
         method = "bh")
dunnTest(Hillq1 ~ Experiment, data = roots_alpha, 
         method = "bh")

########################## Rhizosphere
# Observed
dunnTest(Observed ~ Inbred_or_Hybrid, data = rhizos_alpha, 
         method = "bh")
dunnTest(Observed ~ Experiment, data = rhizos_alpha, 
         method = "bh")
# Shannon
dunnTest(Shannon ~ Inbred_or_Hybrid, data = rhizos_alpha, 
         method = "bh")
dunnTest(Shannon ~ Experiment, data = rhizos_alpha, 
         method = "bh")
# Hills q1
dunnTest(Hillq1 ~ Inbred_or_Hybrid, data = rhizos_alpha, 
         method = "bh")
dunnTest(Hillq1 ~ Experiment, data = rhizos_alpha, 
         method = "bh")


