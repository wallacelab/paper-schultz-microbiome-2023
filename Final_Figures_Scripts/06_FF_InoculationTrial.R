library(ggplot2)
library(data.table)

setwd("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Inoculation_Trial")

MM_data <- read.csv("InoculatedMM_DriedMasses.csv", header = TRUE, sep = "\t")

ggplot(data = MM_data, aes(x =Inoculation, y = Above, fill = Inoculation )) + geom_point(aes(colour = Inoculation)) + ylab("Dried Mass (g)") +
  xlab("Inoculation Microbes") + ggtitle("Inoculated MiniMaize Dried Above Ground Mass")

ggplot(data = MM_data, aes(x =Inoculation, y = Below, fill = Inoculation )) + geom_point(aes(colour = Inoculation)) + ylab("Dried Mass (g)") +
  xlab("Inoculation Microbes") + ggtitle("Inoculated MiniMaize Dried Below Ground Mass")

f1_data = MM_data[MM_data$Inoculation %like% "F1",]
B73_data = MM_data[MM_data$Inoculation %like% "B73",]
MO17_data = MM_data[MM_data$Inoculation %like% "MO17",]
Con_data = MM_data[MM_data$Inoculation %like% "Control",]

# Above Ground vs Con
t.test(f1_data$Above,Con_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(B73_data$Above,Con_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Above,Con_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)

# Above vs F1
t.test(B73_data$Above,f1_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Above,f1_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)

# Below Ground vs Con
t.test(f1_data$Below,Con_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(B73_data$Below,Con_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Below,Con_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)

# Below vs F1
t.test(B73_data$Below,f1_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Below,f1_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)


##################### Remove 3 outliers for roots

New_MM_data <- read.csv("InoculatedMM_DriedMasses_RmvOut.csv", header = TRUE, sep = "\t")

ggplot(data = New_MM_data, aes(x =Inoculation, y = Above, fill = Inoculation )) + geom_point(aes(colour = Inoculation)) + ylab("Dried Mass (g)") +
  xlab("Inoculation Microbes") + ggtitle("Inoculated MiniMaize Dried Above Ground Mass")

ggplot(data = New_MM_data, aes(x =Inoculation, y = Below, fill = Inoculation )) + geom_point(aes(colour = Inoculation)) + ylab("Dried Mass (g)") +
  xlab("Inoculation Microbes") + ggtitle("Inoculated MiniMaize Dried Below Ground Mass")

f1_data = New_MM_data[New_MM_data$Inoculation %like% "F1",]
B73_data = New_MM_data[New_MM_data$Inoculation %like% "B73",]
MO17_data = New_MM_data[New_MM_data$Inoculation %like% "MO17",]
Con_data = New_MM_data[New_MM_data$Inoculation %like% "Control",]

# Above Ground vs Con
t.test(f1_data$Above,Con_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(B73_data$Above,Con_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Above,Con_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)

# Above vs F1
t.test(B73_data$Above,f1_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Above,f1_data$Above, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)

# Below Ground vs Con
t.test(f1_data$Below,Con_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(B73_data$Below,Con_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Below,Con_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)

# Below vs F1
t.test(B73_data$Below,f1_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)
t.test(MO17_data$Below,f1_data$Below, mu=0, alt="two.sided", conf=0.95, var.eq = F, paired=F)


# PAG FIGURES
library(ggpubr)
mycompares = list(c("F1", "Control"),c("MO17", "Control"),c("B73", "Control"))

New_MM_data$Inoculation <- factor(New_MM_data$Inoculation, levels = c("B73","MO17","F1","Control"))

abv_f <- ggplot(data = New_MM_data, aes(x =Inoculation, y = Above, fill = Inoculation )) + geom_point(aes(colour = Inoculation), size = 5) + ylab("Dried Mass (g)") +
  xlab("Inoculation Source") + ylab("Dried Above Ground Mass (g)") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) + 
  stat_compare_means(comparisons = mycompares, method = "t.test",size = 8) + 
  stat_compare_means(comparisons = mycompares, method = "t.test",size = 8, 
                     label = "p.signif", hide.ns=T, bracket.size = 0, vjust = -1) +
  theme(legend.position="none") + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
#+ annotate(geom = "text",x=3, y=75,label = "*")

blw_f <- ggplot(data = New_MM_data, aes(x =Inoculation, y = Below, fill = Inoculation )) + geom_point(aes(colour = Inoculation), size = 5) + ylab("Dried Mass (g)") +
  xlab("Inoculation Source") + ylab("Dried Below Ground Mass (g)")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) + 
  stat_compare_means(comparisons = mycompares, method = "t.test",size = 8) + 
  stat_compare_means(comparisons = mycompares, method = "t.test",size = 8, 
                     label = "p.signif", hide.ns=T, bracket.size = 0, vjust = -1) +
  theme(legend.position="none") + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

MM_inoc <- ggarrange(abv_f,blw_f, ncol = 2, nrow = 1)
MM_inoc

ggsave("Figure6_MM_D2.png", MM_inoc, path = "/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Scripts/Draft2_Figures/D2_Figures",
       width = 10, height = 8, dpi = 600, device = "png")
