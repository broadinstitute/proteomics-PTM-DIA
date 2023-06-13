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
library("RColorBrewer")

#Point to the directory that ends with AllSpikes_Info
setwd('Z:\\Helium_Tan\\R02_PTMDIA\\ASMS\\CVs')



list_files = list.files(pattern = "_output",recursive = T)


LL <- lapply(list_files, FUN = function(L){
  df = fread(L)
})


LL2 <- bind_rows(LL)



LL2 <- LL2 %>% filter(Workflow != 'sn_HybridLibrary' & Workflow != 'diann_Hybrid')
LL2 <- LL2%>%
  mutate(across('Workflow', str_replace, 'sn_DDALibrary', 'Spectronaut_DDALibrary')) %>%
  mutate(across('Workflow', str_replace, 'sn_directDIA', 'Spectronaut_directDIA')) %>%
  mutate(across('Workflow', str_replace, 'diann_DDALibrary', 'DIANN_DDALibrary'))%>%
  mutate(across('Workflow', str_replace, 'diann_Predicted', 'DIANN_PredictedLibrary'))
  

ggplot(data=subset(LL2, !is.na(`Percent CV`)), aes(x=as.factor(Spike), y=`Percent CV`,fill = Instrument), na.rm = TRUE) +
  geom_boxplot(position = position_dodge2(preserve = "single"))+
  ylim(0, 100)+
  xlab('Spike')+
  labs(fill = 'Workflow')+
  scale_fill_brewer(palette = 'Set3')+
  # scale_fill_brewer(palette = 'Dark2')+
  theme(text = element_text(size = 14))
  


comb <- ggplot(data = subset(Combined, !is.na(`Percent CV`)), aes(x= as.factor(Spike), y = `Percent CV`, fill = factor(Filter, levels = c('Unfiltered', 'Stringent'))), na.rm = TRUE)+
  geom_boxplot(position = position_dodge2(preserve = "single"))+
  ylim(0, 100)+
  xlab('Spike')+
  facet_wrap(~Workflow, nrow = 1)+
  labs(fill = 'PTM Site Localization Filter')
  
comb
