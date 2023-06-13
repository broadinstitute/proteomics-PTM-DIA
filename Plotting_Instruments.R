require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)
library(UpSetR)
library(ggvenn)
library(readxl)
library(VennDiagram)
install.packages("pacman")
install.packages("VennDiagram")
library(RColorBrewer)

workflow_order <- c("Spectronaut_DDALibrary","Spectronaut_directDIA","DIANN_DDALibrary","DIANN_PredictedLibrary")

###Time###
time <- read.table(file = "clipboard", sep = '\t', header = TRUE)
View(time)

ggplot(time, aes(fill = Workflow,x=Workflow, y=Time)) +
  geom_bar(position='dodge', stat='identity') +
  xlab("Workflow")+
  ylab("Time (hours)")+
  coord_flip()

### Quantification Completeness ###
completeness <- read.table(file = "clipboard", sep = '\t', header = TRUE)
View(completeness)
ggplot(completeness, aes(x=factor(Stringency, level= c('Unlocalized','Localized')),y=Perc, group = Quantified, color = as.factor(Quantified))) + 
  geom_line()+
  geom_point() +
  facet_wrap(~Instrument, nrow= 1) +  #Change to instrument as needed
  ylim(0,100) +
  ylab('Percent Data >= Quantification Threshold') +
  xlab('Localization Stringency') +
  labs(color = 'Quantification Threshold')+
  scale_fill_brewer(palette = 'Set3')+
  theme(legend.position="top")+
  theme(text = element_text(size = 20)) +
  theme(strip.text.x = element_text(size = 13.5))+
  theme(legend.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))
  


### Phosphopeptides and Phosphosites ###
phosphopeptides <- read.table(file = "clipboard", sep = '\t', header = TRUE)
View(phosphopeptides)

ggplot(phosphopeptides, aes(x=factor(Localization, level = c('Unlocalized','Localized')),fill = Instrument,y=Phosphopeptides)) +   #Can change to instrument
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(ymin = Phosphopeptides - Peptide_SD, ymax= Phosphopeptides + Peptide_SD), width = 0.2, position = position_dodge(0.9))+
  xlab('Localization Stringency') +
  ylab('Mean Phosphopeptide Counts') +
  theme(legend.position="none")+
  scale_fill_brewer(palette = 'Set3')+
  theme(text = element_text(size = 20))
  

ggplot(phosphopeptides, aes(x=factor(Localization, level = c('Unlocalized','Localized')),fill = Instrument,y=Phosphosites)) +      #Can change to instrument
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(ymin = Phosphosites - Site_SD, ymax= Phosphosites + Site_SD), width = 0.2, position = position_dodge(0.9))+
  xlab('Localization Stringency') +
  ylab('Mean Phosphosite Counts') +
  scale_fill_brewer(palette = 'Set3')+
  theme(legend.position = 'none')+
  theme(text = element_text(size = 20))


### Site Localization ###
localization <- read.table(file = "clipboard", sep = '\t', header = TRUE) 
localization <- localization %>% drop_na(Perc)

ggplot(localization, aes(x= as.factor(Spike),y=Perc, fill = Instrument)) +
  geom_bar(position = position_dodge2(preserve = "total"), stat='identity') +
  facet_wrap('Instrument', scales = "free", nrow =1) +
  xlab('Spike') +
  ylab('Percent Phosphosites')+
  ylim(0,100)+
  scale_fill_brewer(palette = 'Set3')+
  theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position = 'none')+
  theme(axis.title = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 13.5))+
  theme(text = element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




