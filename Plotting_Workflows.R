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

workflow_order <- c("Spectronaut_directDIA","Spectronaut_DDALibrary","DIANN_PredictedLibrary","DIANN_DDALibrary")

### Time ###
time <- read.table(file = "clipboard", sep = '\t', header = TRUE)

ggplot(time, aes(fill = Workflow, x=Workflow, y=Time)) +
  geom_bar(position='dodge', stat='identity') +
  xlab("Workflow")+
  ylab("Time (hours)")+
  scale_fill_brewer(palette = 'Dark2')+
  coord_flip()+
  theme(text = element_text(size = 20))
  

### Quantification Completeness ###
completeness <- read.table(file = "clipboard", sep = '\t', header = TRUE)
View(completeness)
ggplot(completeness, aes(x=factor(Stringency, level= c('Unlocalized','Localized')),y=Perc, group = Quantified, color = as.factor(Quantified))) + 
  geom_line()+
  geom_point() +
  facet_wrap(~Workflow, nrow = 1) +  #Change to instrument as needed
  ylim(0,100) +
  ylab('Percent Data >= Quantification Threshold') +
  xlab('Localization Stringency') +
  labs(color = 'Quantification Threshold')+
  theme(legend.position="top")+
  theme(text = element_text(size = 20)) +
  theme(strip.text.x = element_text(size = 13.5))+
  theme(legend.title = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))
  

  
  
### Phosphopeptides and Phosphosites ###
phosphopeptides <- read.table(file = "clipboard", sep = '\t', header = TRUE)
View(phosphopeptides)

# phosphopeptides <- phosphopeptides %>% filter(Localization == 'Unfiltered')
ggplot(phosphopeptides, aes(x= factor(Localization, level = c('Unlocalized','Localized')),fill = Workflow,y=Phosphopeptides)) +   
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(ymin = Phosphopeptides - Peptide_SD, ymax= Phosphopeptides + Peptide_SD), width = 0.2, position = position_dodge(0.9))+
  xlab('Localization Stringency') +
  ylab('Mean Phosphopeptide Counts') +
  scale_fill_brewer(palette = 'Dark2')+
  theme(legend.position="top")+
  theme(legend.position =  'none')+
  theme(text = element_text(size = 20))

  

ggplot(phosphopeptides, aes(x=factor(Localization, level = c('Unlocalized','Localized')),fill = Workflow,y=Phosphosites)) +      
  geom_bar(position='dodge', stat='identity') +
  geom_errorbar(aes(ymin = Phosphosites - Site_SD, ymax= Phosphosites + Site_SD), width = 0.2, position = position_dodge(0.9))+
  xlab('Localization Stringency') +
  ylab('Mean Phosphosite Counts') +
  scale_fill_brewer(palette = 'Dark2')+
  theme(legend.position =  'none')+
  theme(text = element_text(size = 14))


### Site Localization ###
localization <- read.table(file = "clipboard", sep = '\t', header = TRUE) 
localization <- localization %>% drop_na(Perc)

ggplot(localization, aes(x= as.factor(Spike),y=Perc, fill = Workflow)) +
  geom_bar(position = position_dodge2(preserve = "total"), stat='identity') +
  facet_wrap('Workflow', scales = "free", nrow = 1) +
  xlab('Spike') +
  ylab('Percent Phosphosites')+
  ylim(0,100) +
  scale_fill_brewer(palette = 'Dark2')+
  theme(legend.position = 'none')+
  theme(axis.title = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 13.5))+
  theme(text = element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  




