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
# 
# otu = otu_table(as.matrix(read.table("phy_data_otu.csv", row.names = 1, header = TRUE, sep = ",")),taxa_are_rows = TRUE)
# taxa = tax_table(as.matrix(read.table("phy_data_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")))
# meta = sample_data(data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ",")), errorIfNULL = FALSE)
# tree = read.tree("phy_data_tree.tree")


otu_table = read.table("phy_data_otu.csv", row.names = 1, header = TRUE, sep = ",")
taxa_table = read.table("phy_data_taxonomy.csv", row.names = 1, header = TRUE, sep = ",")
meta = sample_data(data.frame(read.table("phy_data_sampdat.csv", row.names = 1, header = TRUE, sep = ",")), errorIfNULL = FALSE)
tree = read.tree("phy_data_tree.tree")

# replace "." with "-" otu dfs
colnames(otu_table)
colnames(otu_table) <- gsub("\\.", "-", colnames(otu_table))
colnames(otu_table)

df <- as.data.frame(t(otu_table))

dfm <- merge(metadata,df, by = "row.names")

#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
dflong <- dfm %>% make_long(Sample_Type)


# Get taxa nums from stuff
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


stalksum <- as.data.frame(taxa_sums(subset_samples(phy_data, Sample_Type=="Stalk")))
rootsum <- as.data.frame(taxa_sums(subset_samples(phy_data, Sample_Type=="Root")))
rhizosum <- as.data.frame(taxa_sums(subset_samples(phy_data, Sample_Type=="Rhizosphere")))


tiss_tax <- cbind(rhizosum,rootsum)
tiss_tax<- cbind(tiss_tax,stalksum)
colnames(tiss_tax) <- c("Rhizosphere","Root","Stalk")

# Next - color by Phylum?
tiss_tax

tiss_tax2 <- cbind(tiss_tax,taxa_table$Phylum)
colnames(tiss_tax2) <- c("Rhizosphere","Root","Stalk","Phylum")
tiss_tax2$Stalk <- ifelse(tiss_tax2$Stalk > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Root <- ifelse(tiss_tax2$Root > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Rhizosphere <- ifelse(tiss_tax2$Rhizosphere > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2

taxa_long <- tiss_tax2 %>% make_long(Rhizosphere,Root,Stalk)
taxa_filt <- filter(taxa_long, node != 0)

# Try with geom_sankey
gs <- ggplot(taxa_filt, aes(x = x
                            , next_x = next_x
                            , node = node
                            , next_node = next_node
                            , fill = factor(node)
                            , label = node)
) + geom_sankey(flow.alpha = 0.5
                , node.color = "black"
                ,show.legend = TRUE) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18)) + theme(axis.title = element_blank()
                                                       , axis.text.y = element_blank()
                                                       , axis.ticks = element_blank()  
                                                       , panel.grid = element_blank(),
                                                       panel.background = element_blank()) + 
   guides(fill=guide_legend(title="Phylum")) + geom_sankey_label(size = 4,
  color = "black",
   fill= "white",
   hjust = -0.5) 


gs
ggsave("Fig5_ASV_Flow_Label_D2.png",
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
       gs, device = "png", width = 12, height = 6, dpi = 600)





# 
# sum(tiss_tax$Rhizosphere)
# sum(tiss_tax$Root)
# sum(tiss_tax$Stalk)
# library(tidyr)
# 
# # Labels and Math
# sank_long <- dflong <- tiss_tax %>% make_long(Rhizosphere,Root,Stalk)
# sank_long$node2 <- ifelse(sank_long$node > 0, "Found in Tissue","Not Found")
# sank_long$next_node2 <- ifelse(sank_long$next_node > 0, "Found in Tissue","Not Found")
# 
# sank_filt <- filter(sank_long, node2 == "Found in Tissue")
# 
# # Try with geom_sankey
# gs <- ggplot(sank_filt, aes(x = x
#                             , next_x = next_x
#                             , node = node2
#                             , next_node = next_node2
#                             , fill = factor(x)
#                             , label = node2)
# ) + geom_sankey(flow.alpha = 0.5
#                 , node.color = "black"
#                 ,show.legend = FALSE) 
# # + geom_sankey_label(size = 4,
# #  color = "black",
# # fill= "white",
# # hjust = -0.5) 
# gs <- gs+ scale_fill_manual(values = c("Rhizosphere" = "purple",
#                                        "Root" = "tan3",
#                                        "Stalk" = "olivedrab")) + 
#   theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 18)) +
#   theme_bw() + theme(legend.position = "none") + theme(axis.title = element_blank()
#                                                        , axis.text.y = element_blank()
#                                                        , axis.ticks = element_blank()  
#                                                        , panel.grid = element_blank())
# 
# gs


# dflong <- tiss_tax %>% make_long(Rhizosphere,Root,Stalk)
# 
# library(networkD3)
# library(dplyr)
# 
# nodes <- data.frame(
#   name=c(as.character(dflong$x), 
#          as.character(dflong$next_x)) %>% unique()
# )
# dflong$IDsource <- match(dflong$x, nodes$name)-1 
# dflong$IDtarget <- match(dflong$next_x, nodes$name)-1
# 
# 
# my_color <- 'd3.scaleOrdinal() .domain(["Rhizosphere","Root","Stalk"]) .range(["purple", "tan", "green"])'
# 
# # Figure
# p <- sankeyNetwork(Links = dflong, Nodes = nodes,
#                    Source = "IDsource", Target = "IDtarget",
#                    Value = "node", NodeID = "name", 
#                    sinksRight=FALSE, fontSize = 12,colourScale=my_color)
# p
