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
df_both_sexes_3 <- separate(df_both_sexes_2, col = 'variant', into = c('chr','pos','other_allele.outcome','effect_allele.outcome'), sep = ':')

##rename columns
colnames(df_both_sexes_3)[colnames(df_both_sexes_3) %in% c("minor_allele","expected_case_minor_AC","n_complete_samples","beta", "se", "pval")] <- c("other_allele.outcome","effect_allele.outcome","eaf.outcome","samplesize.outcome","beta.outcome","se.outcome", "pval.outcome")
df_both_sexes_3$outcome = "ICH"
df_both_sexes_3$originalname.outcome = "Intracranial Hemmorahge"
df_both_sexes_3$mr_keep.outcome = "TRUE"
df_both_sexes_3$data_source.outcome = "Neale"

##Create the subsets correlating with the exposure snp
#systolic bp
exposure_sysbp_id = "ukb-b-6503"
exposure_sysbp_snps = extract_instruments(exposure_sysbp_id)
exposure_sysbp_snps = clump_data(exposure_sysbp_snps)
df_sysbp <- subset(df_both_sexes_3, exposure_sysbp_snps$chr.exposure == df_both_sexes_3$chr & exposure_sysbp_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#diastolic bp
exposure_diabp_id = "ukb-b-7992"
exposure_diabp_snps = extract_instruments(exposure_diabp_id)
exposure_diabp_snps = clump_data(exposure_diabp_snps)
df_diabp <- subset(df_both_sexes_3, exposure_diabp_snps$chr.exposure == df_both_sexes_3$chr & exposure_diabp_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#bmi
exposure_bmi_id = "ukb-a-248"
exposure_bmi_snps = extract_instruments(exposure_bmi_id)
exposure_bmi_snps = clump_data(exposure_bmi_snps)
df_bmi <- subset(df_both_sexes_3, exposure_bmi_snps$chr.exposure == df_both_sexes_3$chr & exposure_bmi_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#platelet count
exposure_pc_id = "ukb-d-30080_irnt"
exposure_pc_snps = extract_instruments(exposure_pc_id)
exposure_pc_snps = clump_data(exposure_pc_snps)
df_pc <- subset(df_both_sexes_3, exposure_pc_snps$chr.exposure == df_both_sexes_3$chr & exposure_pc_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#red blood cell count
exposure_rbccount_id = "ieu-a-275"
exposure_rbccount_snps = extract_instruments(exposure_rbccount_id)
exposure_rbccount_snps = clump_data(exposure_rbccount_snps)
df_rbccount <- subset(df_both_sexes_3, exposure_rbccount_snps$chr.exposure == df_both_sexes_3$chr & exposure_rbccount_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#platelet volume
exposure_plateletvol_id = "ieu-a-1006"
exposure_plateletvol_snps = extract_instruments(exposure_plateletvol_id)
exposure_plateletvol_snps = clump_data(exposure_plateletvol_snps)
df_platletvol <- subset(df_both_sexes_3, exposure_plateletvol_snps$chr.exposure == df_both_sexes_3$chr & exposure_plateletvol_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#red blood cell distribution width (indicator of nutritional definciency)
exposure_rbcdw_id = "ukb-d-30070_irnt"
exposure_rbcdw_snps = extract_instruments(exposure_rbcdw_id)
exposure_rbcdw_snps = clump_data(exposure_rbcdw_snps)
df_rbcdw <- subset(df_both_sexes_3, exposure_rbcdw_snps$chr.exposure == df_both_sexes_3$chr & exposure_rbcdw_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#past tobacco smoking
exposure_tobacco_past_id = "ukb-b-2134"
exposure_tobacco_past_snps = extract_instruments(exposure_tobacco_past_id)
exposure_tobacco_past_snps = clump_data(exposure_tobacco_past_snps)
df_tobacco_past <- subset(df_both_sexes_3, exposure_tobacco_past_snps$chr.exposure == df_both_sexes_3$chr & exposure_tobacco_past_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)

#current tobacco smoking
exposure_tobacco_current_id = "ukb-b-223"
exposure_tobacco_current_snps = extract_instruments(exposure_tobacco_current_id)
exposure_tobacco_current_snps = clump_data(exposure_tobacco_current_snps)
df_tobacco_current <- subset(df_both_sexes_3, exposure_tobacco_current_snps$chr.exposure == df_both_sexes_3$chr & exposure_tobacco_current_snps$pos.exposure == df_both_sexes_3$pos, select=chr:data_source.outcome)
