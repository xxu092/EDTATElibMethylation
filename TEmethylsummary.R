setwd("/bigdata/stajichlab/xxu092/M.cicadina/EDTATElibMethylation/")
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr) #for spearman correlation

all <- read_tsv("TElib_con_5mc_c.tsv",col_names = c("CHROM","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE","GROUP","No.CpG"))


dim(all)
head(all)
positive <- read_tsv("TElib_con_5mc_20.tsv",col_names = c("CHROM","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE","GROUP","methylC"))
dim(positive)
head(positive)
TE_5mc <- left_join(all, positive, by=c("CHROM","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE","GROUP"))
dim(TE_5mc)
head(TE_5mc)
TE_5mc <- mutate(TE_5mc,perc = methylC/ No.CpG)
head(TE_5mc)

new_TE_5mc <- TE_5mc %>% mutate(TE_type = str_replace(GROUP, "^.+Name=([^;]+);Classification=.+$","\\1")) %>%
  mutate(TE_family = str_replace(GROUP,"^.+Classification=([^;]+);Sequence_ontology=.+$","\\1")) %>% 
  mutate(length = END - START + 1) %>%
  select(-c(GROUP,SOURCE,TYPE,PHASE,SCORE))
head(new_TE_5mc)
write_tsv(new_TE_5mc,"TEmethylsummary.tsv")

## continue from here add length data 
positive_LINE <- filter(new_LINE_5mc, perc >= 0.25)
positive_LINE <- mutate(positive_LINE,halfmethyl=1)
negative_LINE <- filter(new_LINE_5mc, perc < 0.25)
negative_LINE <- mutate(negative_LINE,halfmethyl=0)
Na_LINE <- filter(new_LINE_5mc, perc == "NaN")
Na_LINE <- mutate(Na_LINE,halfmethyl=0)
head(positive_LINE)
head(negative_LINE)
head(Na_LINE)

final_LINE <- full_join(positive_LINE,negative_LINE)
final_LINE <- full_join(final_LINE,Na_LINE)
head(final_LINE)
dim(final_LINE)
final_LINE <- final_LINE %>% 
  group_by(TE_type) %>%
  summarize(n(),sum(halfmethyl))
colnames(final_LINE) <- c("TE_type", "sum_5mc", "positive_5mc")
head(final_LINE)
final_LINE <- mutate(final_LINE,methylated_repeats=final_LINE$positive_5mc/final_LINE$sum_5mc)
print(final_LINE)
families <- read_tsv("/bigdata/stajichlab/xxu092/M.cicadina/LINE/EDTA_LINE/LINE_families_EDTA.tsv", col_names = TRUE)
head(families)
LINE_methyl <- full_join(families,final_LINE, by = c("NAME"="TE_type"))

write_tsv(LINE_methyl,"LINE_methylation_summary.tsv")

ggplot(LINE_methyl, aes(x=counts, y=methylated_repeats)) +
  geom_point()

ggplot(LINE_methyl, aes(x=methylated_repeats, y=counts)) +
  geom_point()

ggplot(LINE_methyl, aes(x=mean_length, y=methylated_repeats)) +
  geom_point()

#calculate spearman correlation
x=LINE_methyl$counts
y=LINE_methyl$methylated_repeats
corr <- cor.test(x, y, method = "spearman")
head(corr)

ggscatter(LINE_methyl, x='x', y='y', cor.coef = TRUE, cor.method = "spearman")
# there is a bug with ggscatter 