Read me for Supplementary Info 

Alpha Diversity
- AlphaD_Supp contains a textfile of R output for Dunn's Test using Kruskal-Wallis multiple comparison.
  p-values are adjusted with the Benjamini-Hochberg method. The comparisons are split by All Tissues, Stalk, Root, and Rhizosphere. 

Beta Diversity
- BetaD_Supp contains an image of a PERMANOVA using weighted Unifrac measures and a type II ANOVA. 

Flow Diagram and Relative Abundance
These figure display our phyloseq data, so no supplementary output exists. 

Common ASVs
Shared_Taxa_Table contains the Taxa_Intersections.tsv file, which shows prescense absence of genus across groups. These intersections are plotted in the upset plot. 

Differential Abundance
The Diff_Abun_ASVs_D2 folder contains csv files of differentially abundant ASVs with DESeq2 output including pvalues. Each file is a complete comparison. Fig4_Tissue_D2.png Shows volcano plots for differential abundance in tissues. 

PICRUST
PICRUST2_tables contain DESEq output for all comparisons. Tables ending in _raw means PICRUST descriptors were not aglomerated by functional pathways, while _agl tables were. Picrust_Differences_Table_Number.csv contains the sum of all significant differences for the comparisons. 

MiniMaize Innoculation
InoculatedMM_DriedMasses.csv contains raw data for the inoculation experiment, with dried mass for above and below ground biomass in grams. 

Other: 
Phyloseq_Objects
This folder contains files needed to create phyloseq objects in R: 
Combined_.qza files are directly from qiime2 after preprocessing (raw). 

phy2_ files are raw data after prevelance filtering, and filtering of Mitochondria and Chloroplasts. 

phy_data files are phy2 where blank samples and taxa have been filtered out. This is the data set used throughout most of the analisys. 

phy_shared files are phy_data files but taxa has been pruned to only those ASVs that appear in all three experiments (Field 1,2, and Greenhouse). This data set was used for Differential Abundance, and PICRUST analysis. 

Krona Plots
HTML Files of Krona Plots are included to interactivly compare taxa. 

