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



myCol <- brewer.pal(5, "Pastel2")


setwd('Z:/Helium_Tan/R02_PTMDIA/Pro')

#Load DIA reports from all 6 workflows
sn_directDIA <- read.table('Spectronaut/directDIA/20230227_104026_PhosphoDIA_R02_Report_directDIA.tsv', sep = '\t', header = TRUE, quote = '') %>% filter(., 'EG.PTMAssayProbability' >= 0.99)
sn_DDAlib <- read.table('Spectronaut/DDALibrary/20230308_091202_R02_PhosphoDIA_Pro_DDALibrary_Report.tsv', sep = '\t', header = TRUE, quote = '') %>% filter(., 'EG.PTMAssayProbability' >= 0.99)
sn_hybrid <- read.table('Spectronaut/Hybrid/20230313_103651_R02_PhosphoDIA_Pro_HybridLibrary_Report.tsv', sep = '\t', header = TRUE, quote = '') %>% filter(., 'EG.PTMAssayProbability' >= 0.99)

diann_predicted <- read.table('DIANN/Predicted_Library/dia_nn/out/report.tsv', sep = '\t', header = TRUE, quote = '') %>% filter(., 'PTM.Site.Confidence' >= 0.51)
diann_DDAlib <- read.table('DIANN/DDA_Library/Three_varmods/dia_nn/out/report.tsv', sep = '\t', header = TRUE, quote = '') %>% filter(., 'PTM.Site.Confidence' >= 0.51)
diann_hybrid <- read.table('DIANN/Combined_DDA_Predicted_Library/dia_nn/out/report.tsv', sep = '\t', header = TRUE, quote = '') %>% filter(., 'PTM.Site.Confidence' >= 0.51)



#Convert Spectronaut modifications to DIA-NN modifications
convert <- function(sequence) {
  
    new <- sequence %>% gsub("\\[\\+80\\]", "\\(UniMod:21\\)",.) %>% gsub("_\\[\\+42\\]","\\(UniMod:1\\)",.) %>% gsub("\\[\\+57\\]","\\(UniMod:4\\)",.) %>%
      gsub("\\[\\+16\\]","\\(UniMod:35)",.) %>% gsub("\\[\\+8\\]","\\(UniMod:259)",.) %>% gsub("\\[\\+10\\]","\\(UniMod:267)",.) %>% gsub("_","",.)
  
    return(new)
    
}
  
  
sn_directDIA['Modified.Sequence'] <- lapply(sn_directDIA['FG.IntMID'],convert)
sn_DDAlib['Modified.Sequence'] <- lapply(sn_DDAlib['FG.IntMID'],convert)


overlap <- list(sn_directDIA = sn_directDIA$Modified.Sequence, sn_DDALibrary = sn_DDAlib$Modified.Sequence, diann_predicted = diann_predicted$Modified.Sequence, diann_DDAlib = diann_DDAlib$Modified.Sequence, diann_hybrid = diann_hybrid$Modified.Sequence)

# upset <- upset(fromList(overlap))                
# upset


v1 <- venn.diagram((overlap), filename = NULL, fill = myCol)
grid.draw(v1)







          
                             
                             
                             
                             
          
