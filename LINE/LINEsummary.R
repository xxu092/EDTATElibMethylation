library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(forcats)
library(stringr)
library(hrbrthemes)
library(viridis)
library(ggridges)
setwd("/bigdata/stajichlab/xxu092/M.cicadina/EDTATElibMethylation/")

TEsummary <- read_tsv("TEmethylsummary.tsv", col_names = TRUE)
LINE <- TEsummary %>% filter(TEsummary$TE_family == "LINE/unknown")

#this mutates CGmethyl to 1 if per is >= 0.5, if not return 0
LINE <- LINE %>% 
  mutate (CGmethyl = ifelse(perc >= 0.25, 1, 0))
#replace NA with 0
LINE$CGmethyl[is.na(LINE$CGmethyl)] <- 0

LINEfam<- LINE %>% group_by(TE_type) %>%
  summarize(counts = n(), 
            methylated_repeats = sum(CGmethyl),
            min_length = min(length),
            max_length = max(length),
            median_length = median(length),
            mean_length = mean(length))

LINEfam <- LINEfam %>%
  mutate(methylated = methylated_repeats / counts)
LINEfam <- LINEfam %>% 
  mutate(family=str_extract(TE_type,"(?<=TE_\\d{4})\\d+"))

write_tsv(LINEfam, "LINE/LINEsummary.tsv")

pdf(file = "LINE/LINEgraph.pdf")


ggplot(LINE, aes( x = length, y = TE_type, fill = TE_type)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  xlab("length") +
  ylab("families")  
  
  ggplot(LINE, aes( x = perc, y = TE_type, fill = TE_type)) +
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position = "none") +
    xlab("methylated C%") +
    ylab("families")  
  
  ggplot(LINE, aes( x = length, y = TE_type)) +
    geom_boxplot(fill="slateblue", alpha = 0.2) +
    xlab("length") +
    ylab("families")
  
  ggplot(LINE, aes( x = perc, y = TE_type)) +
    geom_boxplot(fill="slateblue", alpha = 0.2) +
    xlab("methylated C%") +
    ylab("families")
  
  ggplot(data = LINE, aes(x = perc, color = TE_type)) +
    geom_bar() +
    facet_wrap(~TE_type)
  
  ggplot(LINE, aes(x=length, y=perc, color=TE_type)) +
    geom_point() +
    xlab("length")+
    ylab("methylated C%") 
    
   
  ggplot(LINEfam, aes(x=counts, y=methylated,color=family)) +
    geom_point() +
    geom_text(
      label=LINEfam$family,
      nudge_x = 0.05, nudge_y = 0.05,
      check_overlap = T
    ) +
    xlab("counts") +
    ylab("methylated element %")
  

  ggplot(LINEfam, aes(x=mean_length, y=methylated, color=family)) +
    geom_point() + 
    geom_text(
      label=LINEfam$family,
      nudge_x = 0.05, nudge_y = 0.05,
      check_overlap = T
    ) +
    xlab("mean_length") +
    ylab("methylated element %")

  dev.print(pdf,"Rgraph.pdf")  
  dev.off()
  