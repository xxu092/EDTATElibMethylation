setwd("/bigdata/stajichlab/xxu092/M.cicadina/EDTATElibMethylation/")
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr) #for spearman correlation

#read in file-TE intersected with all CpG sites
all <- read_tsv("TElib_con_5mc_c.tsv",
  col_names = c("CHROM","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE","GROUP","No.CpG"))


dim(all)
head(all)
#read in file -TE intersected with only positively methylated CpG sites
positive <- read_tsv("TElib_con_5mc_20.tsv", 
  col_names = c("CHROM","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE","GROUP","methylC"))
dim(positive)
head(positive)
#join all table and positive table
TE_5mc <- left_join(all, positive, by=c("CHROM","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE","GROUP"))
dim(TE_5mc)
head(TE_5mc)

#calculate percentage of Cs CpG in one element that's methylated
TE_5mc <- mutate(TE_5mc,perc = methylC/ No.CpG)
head(TE_5mc)

#make a new table with length information
new_TE_5mc <- TE_5mc %>% mutate(TE_type = str_replace(GROUP, "^.+Name=([^;]+);Classification=.+$","\\1")) %>%
  mutate(TE_family = str_replace(GROUP,"^.+Classification=([^;]+);Sequence_ontology=.+$","\\1")) %>% 
  mutate(length = END - START + 1) %>%
  select(-c(GROUP,SOURCE,TYPE,PHASE,SCORE))
head(new_TE_5mc)
write_tsv(new_TE_5mc,"TEmethylsummary.tsv")

