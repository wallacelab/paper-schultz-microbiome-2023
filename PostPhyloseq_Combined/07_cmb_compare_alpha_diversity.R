library(ggplot2)

#Read in FaithPD Files

Key_File = Key_File = read.csv("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied/Combined_Key.csv", sep = "\t")

rhizos_fpd = read.csv("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied/Cmb_rhizosphere-core-metrics-rarefied-results/Cmb_rhizosphere_faithPD/alpha-diversity.tsv", sep = "\t")
roots_fpd = read.csv("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied/Cmb_roots-core-metrics-rarefied-results/Cmb_roots_faithPD/alpha-diversity.tsv", sep = "\t")
stalk_fpd = read.csv("/home/coreyschultz/1.Projects/2.Heterosis.Microbiome/Maize_Het_Microbiome_CS/Combined_CS/Combined_Results/CmbFiltandRarefied/Cmb_stalks-core-metrics-rarefied-results/Cmb_stalks_faithPD/alpha-diversity.tsv", sep = "\t")

names(rhizos_fpd)=c("SampleID", "faith_pd")
names(roots_fpd)=c("SampleID", "faith_pd")
names(stalk_fpd)=c("SampleID", "faith_pd")

#Which groups are more diverse?
rhizos_data=merge(rhizos_fpd,Key_File, all.x = TRUE)
roots_data=merge(roots_fpd,Key_File, all.x = TRUE)
stalk_data=merge(stalk_fpd,Key_File, all.x = TRUE)

aggregate(rhizos_data$faith_pd, by=list(Inbred_or_Hybrid=rhizos_data$Inbred_or_Hybrid), FUN = sum)
aggregate(roots_data$faith_pd, by=list(Inbred_or_Hybrid=roots_data$Inbred_or_Hybrid), FUN = sum)
aggregate(stalk_data$faith_pd, by=list(Inbred_or_Hybrid=stalk_data$Inbred_or_Hybrid), FUN = sum)


#Plotting Inbred or Hybrid

Rhizosphere_plots=ggplot(rhizos_data,aes(x=Inbred_or_Hybrid,y=faith_pd))+
  geom_violin(aes(fill=Inbred_or_Hybrid ))+
  ggtitle("Combined Rhizosphere Faith PD")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=(element_text(angle=90,hjust=1)))+
  theme(legend.position="none")+
  labs(x='Inbred_or_Hybrid',fill='Inbred_or_Hybrid')

Root_plots=ggplot(roots_data,aes(x=Inbred_or_Hybrid,y=faith_pd))+
  geom_violin(aes(fill=Inbred_or_Hybrid ))+
  ggtitle("Combined Root Faith PD")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=(element_text(angle=90,hjust=1)))+
  theme(legend.position="none")+
  labs(x='Inbred_or_Hybrid',fill='Inbred_or_Hybrid')

Stalk_plots=ggplot(stalk_data,aes(x=Inbred_or_Hybrid,y=faith_pd))+
  geom_violin(aes(fill=Inbred_or_Hybrid ))+
  ggtitle("Combined Stalk Faith PD")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=(element_text(angle=90,hjust=1)))+
  theme(legend.position="none")+
  labs(x='Inbred_or_Hybrid',fill='Inbred_or_Hybrid')

Rhizosphere_plots
Root_plots
Stalk_plots


#Plotting - By Experiments

aggregate(rhizos_data$faith_pd, by=list(Experiment=rhizos_data$Experiment), FUN = sum)
aggregate(roots_data$faith_pd, by=list(Experiment=roots_data$Experiment), FUN = sum)
aggregate(stalk_data$faith_pd, by=list(Experiment=stalk_data$Experiment), FUN = sum)


Rhizosphere_plots=ggplot(rhizos_data,aes(x=Experiment,y=faith_pd))+
  geom_violin(aes(fill=Experiment ))+
  ggtitle("Combined Rhizosphere Faith PD")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=(element_text(angle=90,hjust=1)))+
  theme(legend.position="none")+
  labs(x='Experiment',fill='Experiment')

Root_plots=ggplot(roots_data,aes(x=Experiment,y=faith_pd))+
  geom_violin(aes(fill=Experiment ))+
  ggtitle("Combined Root Faith PD")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=(element_text(angle=90,hjust=1)))+
  theme(legend.position="none")+
  labs(x='Experiment',fill='Experiment')

Stalk_plots=ggplot(stalk_data,aes(x=Experiment,y=faith_pd))+
  geom_violin(aes(fill=Experiment ))+
  ggtitle("Combined Stalk Faith PD")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=(element_text(angle=90,hjust=1)))+
  theme(legend.position="none")+
  labs(x='Experiment',fill='Experiment')

Rhizosphere_plots
Root_plots
Stalk_plots






