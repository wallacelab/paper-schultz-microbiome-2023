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

my_labels <- c("Verrucomicrobia","Proteobacteria","Planctomycetes",
               "Chloroflexi","Actinobacteria","Acidobacteria", 
               "Bacteroidetes", "Cyanobacteria", "Firmicutes",
               "Gemmatimonadetes")

taxa_filt$labels <- ifelse((taxa_filt$node %in% my_labels & taxa_filt$x=='Rhizosphere'), taxa_filt$node, NA)

# Try with geom_sankey
gs <- ggplot(taxa_filt, aes(x = x
                            , next_x = next_x
                            , node = node
                            , next_node = next_node
                            , fill = factor(node)
                            , label = labels)
) + geom_sankey(flow.alpha = 0.5
                , node.color = "black"
                ,show.legend = TRUE) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18)) + theme(axis.title.x = element_blank()
                                                       , axis.title.y = element_text(size = 16, face = "bold")
                                                       , axis.text.y = element_blank()
                                                       , axis.ticks = element_blank()  
                                                       , panel.grid = element_blank(),
                                                       panel.background = element_blank()) + 
  guides(fill="none") + geom_sankey_label(size = 4,
                        color = "black",
                        fill= "white") + ylab("Flow of ASVs")


gs
# ggsave("Fig5_ASV_Flow_Label_D2.png",
#        path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
#        gs, device = "png", width = 12, height = 6, dpi = 600)



phy_data
phy.c <- transform_sample_counts(phy_data, function(x) x/sum(x))
phy.c <- subset_samples(phy.c, Sample_Type != "Soil")
phy.melt <- psmelt(phy.c)

relbar <- ggplot(data = phy.melt, mapping = aes_string(x = "Sample_Type",y = "Abundance")) +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill") 

relbar <- relbar + theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = "white",color="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank() ) + ylab("Relative Abundance of Reads")
  

relbar


flow_and_relabund <- ggarrange(gs,relbar, labels = c("A","B"))

ggsave("Fig5_Flow_and_Bar_D2.png",
       path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
       flow_and_relabund, device = "png", width = 12, height = 6, dpi = 600)
########## Compare with node size at 100%
tiss_tax

tiss_rel <- tiss_tax

# Normalize counts to 1,000,000
tiss_rel[1:3] <- lapply(tiss_rel[1:3], function(x) floor(x/sum(x) * 10000))
sum(tiss_rel$Rhizosphere)
sum(tiss_rel$Root)
sum(tiss_rel$Stalk)
library(tidyverse)
huge_df <- data.frame(ASV=NA,Root=NA,Rhizosphere=NA,Stalk=NA)

# Turn into 3 million rows
for(c in 1:ncol(tiss_rel)){
  for(r in 1:nrow(tiss_rel)){
    count = tiss_rel[r,c]
    rname = rownames(tiss_rel)[r]
    if(c == 1){
      newrow <- c((rownames(tiss_rel)[r]),1,0,0)
    }
    if(c == 2){
      newrow <- c((rownames(tiss_rel)[r]),0,1,0)
    }
    if(c == 3){
      newrow <- c((rownames(tiss_rel)[r]),0,1,0)
    }
    if(count > 0){
    for(i in 1:count){
      huge_df <- rbind(huge_df, rep(newrow))}
    }
    
    #rerun(count,huge_df <- rbind(huge_df, newrow))
    #huge_df <- rbind(huge_df, rep(newrow), count)
    #do.call("rbind", replicate(count, newrow, simplify = FALSE))
  }

}
huge_df <- huge_df[-1,]
print(nrow(huge_df))
head(huge_df)

taxa_table <- tibble::rownames_to_column(taxa_table, "ASV")

tiss_tax2 <- left_join(huge_df,taxa_table,by = "ASV")


tiss_tax2$Stalk <- ifelse(tiss_tax2$Stalk > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Root <- ifelse(tiss_tax2$Root > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2$Rhizosphere <- ifelse(tiss_tax2$Rhizosphere > 0, as.character(tiss_tax2$Phylum),0)
tiss_tax2

taxa_long <- tiss_tax2 %>% make_long(Rhizosphere,Root,Stalk)
taxa_filt <- filter(taxa_long, node != 0)

my_labels <- c("Verrucomicrobia","Proteobacteria","Planctomycetes",
               "Chloroflexi","Actinobacteria","Acidobacteria", 
               "Bacteroidetes", "Cyanobacteria", "Firmicutes",
               "Gemmatimonadetes")

taxa_filt$labels <- ifelse((taxa_filt$node %in% my_labels & taxa_filt$x=='Rhizosphere'), taxa_filt$node, NA)

# Try with geom_sankey
gs <- ggplot(taxa_filt, aes(x = x
                            , next_x = next_x
                            , node = node
                            , next_node = next_node
                            , fill = factor(node)
                            , label = labels)
) + geom_sankey(flow.alpha = 0.5
                , node.color = "black"
                ,show.legend = TRUE) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18)) + theme(axis.title = element_blank()
                                                       , axis.text.y = element_blank()
                                                       , axis.ticks = element_blank()  
                                                       , panel.grid = element_blank(),
                                                       panel.background = element_blank()) + 
  guides(fill="none") + geom_sankey_label(size = 4,
                                          color = "black",
                                          fill= "white")


gs

# Regular Relative Abundance graph
phy_data
phy.c <- transform_sample_counts(phy_data, function(x) x/sum(x))
phy.c <- subset_samples(phy.c, Sample_Type != "Soil")
phy.melt <- psmelt(phy.c)

ggplot(data = phy.melt, mapping = aes_string(x = "Sample_Type",y = "Abundance")) +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill")
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
