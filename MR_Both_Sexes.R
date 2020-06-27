#Intracranial Aneurysm MR Analysis 

#load data table package and set working directory
library(data.table)
library(tidyr)
setwd('~/Documents/MR_Code')

#load the data file
df_both_sexes_1 = fread('I61.gwas.imputed_v3.both_sexes.tsv')

#remove all rows where low_confidence_variant = true
df_both_sexes_2 <- subset(df_both_sexes_1, low_confidence_variant == 'FALSE')

#separate first colunm using ":" as delimiter
#give columns actual names 
df_both_sexes_3 <- separate(df_both_sexes_2, col = 'variant', into = c('c1','c2','c3','c4'), sep = ':')