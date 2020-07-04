#Intracranial Aneurysm MR Analysis 

#load data table package and set working directory
library(data.table)
library(tidyr)
library(dplyr)
setwd('~/Documents/MR_Code')

##load the data file
df_both_sexes_1 = fread('I61.gwas.imputed_v3.both_sexes.tsv')


##remove all rows where low_confidence_variant = true
df_both_sexes_2 <- subset(df_both_sexes_1, low_confidence_variant == 'FALSE')

##separate first column using ":" as delimiter
##give columns actual names 
df_both_sexes_3 <- separate(df_both_sexes_2, col = 'variant', into = c('chr','pos','Reference variant','Alternative variant'), sep = ':')

##rename columns
#effect_allele.outcome --> replaced "minor_allele"
#other_allele.outcome --> replaced "Alternative variant"??? 
#eaf.outcome --> replaced "expected_case_minor_AC"
#samplesize.outcome --> replaced "n_complete_samples"
#beta.outcome --> replaced "beta"
colnames(df_both_sexes_3)[colnames(df_both_sexes_3) %in% c("Alternative variant","minor_allele","expected_case_minor_AC","n_complete_samples","beta")] <- c("other_allele.outcome","effect_allele.outcome","eaf.outcome","samplesize.outcome","beta.outcome")

##creating SNP column 
#delete "minor allele"
