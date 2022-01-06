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

abv_f <- ggplot(data = New_MM_data, aes(x =Inoculation, y = Above, fill = Inoculation )) + geom_point(aes(colour = Inoculation), size = 4) + ylab("Dried Mass (g)") +
  xlab("Inoculation Microbes") + ggtitle("Dried Above Ground Mass") +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) + 
  stat_compare_means(comparisons = mycompares, method = "t.test")

blw_f <- ggplot(data = New_MM_data, aes(x =Inoculation, y = Below, fill = Inoculation )) + geom_point(aes(colour = Inoculation), size = 4) + ylab("Dried Mass (g)") +
  xlab("Inoculation Microbes") + ggtitle("Dried Below Ground Mass")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18)) + 
  stat_compare_means(comparisons = mycompares, method = "t.test")

MM_inoc <- ggarrange(abv_f,blw_f, ncol = 2, nrow = 1)

ggsave("MM_inoc_PAG", MM_inoc, width = 10, height = 5, dpi = 300, device = "png")
