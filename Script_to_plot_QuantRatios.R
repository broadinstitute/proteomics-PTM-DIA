require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library("RColorBrewer")
library(ggpubr)


# setwd("Z:\\Helium_Tan\\R02_PTMDIA\\Exploris\\DIANN\\DDALibrary\\dia_nn\\out\\Exploris_AllSpikes_Replicates")
setwd("Z:\\Helium_Tan\\R02_PTMDIA\\FAIMS\\DIANN\\DDALibrary\\dia_nn\\out\\Exploris_FAIMS_AllSpikes_Replicates")
# setwd("Z:\\Helium_Tan\\R02_PTMDIA\\Pro\\DIANN\\DDA_Library\\Three_varmods\\dia_nn\\out\\timsTOF_Pro_AllSpikes_Replicates")
# setwd("Z:\\Helium_Tan\\R02_PTMDIA\\SCP\\DIANN\\DDALibrary\\dia_nn\\out\\TimsTOF_SCP_AllSpikes_Replicates")

# setwd('Z:\\Helium_Tan\\R02_PTMDIA\\Pro\\DIANN\\Predicted_Library\\dia_nn\\out\\Pro_AllSpikes_Replicates')
# setwd('Z:\\Helium_Tan\\R02_PTMDIA\\Pro\\Spectronaut\\DDALibrary\\Pro_AllSpikes_Replicates')
# setwd('Z:\\Helium_Tan\\R02_PTMDIA\\Pro\\Spectronaut\\directDIA\\Pro_AllSpikes_Replicates')


list_files = list.files(pattern = "replicates_output",recursive = T)

#Bind rows into one large dataframe
LL <- lapply(list_files, FUN = function(L){
  df = fread(L)
})

LL2 <- bind_rows(LL)

#Create new value specifically for expected ratios
ratios_for_plotting <- LL2 %>%
  select(Spike, `Expected Ratio`) %>%
  distinct() %>%
  arrange(Spike)


ratios_for_plotting$Spike  <- factor(ratios_for_plotting$Spike, levels = ratios_for_plotting$Spike)
LL2$Spike <- factor(LL2$Spike, levels = ratios_for_plotting$Spike)
LL2$logRatio <- log10(1/LL2$`Actual Ratio`)
ratios_for_plotting$logExpRatio <- log10 (1/ratios_for_plotting$`Expected Ratio`)

#Plot quant ratios as bar plots for each spike level
ggplot(LL2 , mapping = aes(y = logRatio, x= Spike)) +
  geom_boxplot() +
  geom_jitter(aes(color=Peptide), cex = 0.5)+
  #scale_y_log10()+
  theme_bw()+
  geom_point(data = ratios_for_plotting, aes(x = Spike, y = logExpRatio),shape = 95, size =10,color= "black", alpha = 1)+
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position = "none")+
  theme(text = element_text(size = 20))



#Find the outliers at each spike level
outliers <- LL2 %>%
  group_by(Spike) %>%
  identify_outliers("logRatio")



#Find and plot correlation between expected ratio and mean actual ratio at each spike level, across the titration curve. the closer the correlation to 1, the more quantitatively accurate/linear.
corr <- LL2 %>% group_by(`Expected Ratio`) %>% summarise('Mean Actual Ratio' = mean(1/`Actual Ratio`))
corr$`Expected Ratio` <- 1/corr$`Expected Ratio`

perfect <- LL2 %>% group_by(`Expected Ratio`) %>% summarise('Perfect Ratio' = mean(`Expected Ratio`))

ggplot(corr, mapping= aes(x = `Expected Ratio`, y =`Mean Actual Ratio`)) +
  geom_line(color = "#7CAE00")+
  geom_point(color = "#7CAE00") +
  geom_line(perfect, mapping = aes(x = `Expected Ratio`, y = `Perfect Ratio`), linetype = "dotted")

model <- lm(formula = `Mean Actual Ratio` ~ `Expected Ratio`, data = corr)
summary(model)
print(cor(corr$`Expected Ratio`,corr$`Mean Actual Ratio`))  #This is the actually published value



#Find correlation between spike level and ratio (is the relationship linear?)
y <- LL2$logRatio
x <- as.numeric (as.character(LL2$Spike))
print (cor (x, y))


#Enumerate peptide replicates at each spike level
LL2 %>%
  group_by(Spike) %>%
  tally()
