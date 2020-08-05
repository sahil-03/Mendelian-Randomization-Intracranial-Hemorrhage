#Intracranial Aneurysm MR Analysis
#load data table package and set working directory
library(data.table)
library(tidyr)
library(dplyr)
library(TwoSampleMR)
setwd('~/Documents/MR_Code/')

##load the downloaded RData object
load("/Users/sahiladhawade/Downloads/rf.rdata")
##load the data file
df_both_sexes_1 = fread('I61.gwas.imputed_v3.both_sexes.tsv')

##remove all rows where low_confidence_variant = true
df_both_sexes_2 <- subset(df_both_sexes_1, low_confidence_variant == 'FALSE')

##separate first column using ":" as delimiter
df_both_sexes_3 <- separate(df_both_sexes_2, col = 'variant', into = c('chr.outcome','pos.outcome','other_allele.outcome','effect_allele.outcome'), sep = ':')

##rename columns
names(df_both_sexes_3)[names(df_both_sexes_3) == "expected_case_minor_AF"] = "eaf.outcome"
names(df_both_sexes_3)[names(df_both_sexes_3) == "n_complete_samples"] = "samplesize.outcome"
names(df_both_sexes_3)[names(df_both_sexes_3) == "beta"] = "beta.outcome"
names(df_both_sexes_3)[names(df_both_sexes_3) == "se"] = "se.outcome"
names(df_both_sexes_3)[names(df_both_sexes_3) == "pval"] = "pval.outcome"
names(df_both_sexes_3)[names(df_both_sexes_3) == "minor_AF"] = "eaf.outcome"
subset = df_both_sexes_3$minor_allele != df_both_sexes_3$effect_allele.outcome
df_both_sexes_3[subset,'eaf.outcome'] = 1 - df_both_sexes_3[subset,'eaf.outcome']

df_both_sexes_3$outcome = "ICH"
df_both_sexes_3$originalname.outcome = "Intracranial Hemmorahge"
df_both_sexes_3$mr_keep.outcome = "TRUE"
df_both_sexes_3$data_source.outcome = "Neale"
df_both_sexes_3$id.outcome = 99999999999999999999
df_both_sexes_3$pos.outcome = as.numeric(df_both_sexes_3$pos.outcome)

#merge with outcome file
merge_exposure_id_with_df_both_sexes_3 = function(df_both_sexes_3, exposure_id, output_file) {
  exposure_sysbp_snps = extract_instruments(exposure_id)
  exposure_sysbp_snps = clump_data(exposure_sysbp_snps)
  
  exposure_sysbp_snps = exposure_sysbp_snps[, c('SNP', 'chr.exposure', 'pos.exposure')]
  df_sysbp = merge(df_both_sexes_3, exposure_sysbp_snps,
                   by.x=c('chr.outcome',  'pos.outcome'),
                   by.y=c('chr.exposure', 'pos.exposure'))
  write.csv(df_sysbp, output_file, row.names=F)
  return(df_sysbp)
}


##Conduct the MR Analysis
exposures_df = data.frame(exposure_id=c("ukb-b-12493", "ukb-b-19953", "ukb-b-7992", "ukb-b-20175", "ebi-a-GCST002223", "ieu-a-781"),
                          output_save_name=c("hypertension_outcome", "BMI_outcome", "diastolic_bp_outcome", "systolic_bp_outcome", "HDL_outcome", "LDL_outcome"))

for (rowi in 1:nrow(exposures_df)) {
  exposure_id = exposures_df$exposure_id[rowi]
  output_save_name = exposures_df$output_save_name[rowi]
  if (file.exists(paste0('scatter_plots/', output_save_name, '.png'))) {
    #next
  } 
  exposure_sysbp_snps = extract_instruments(exposure_id)
  # dataframe of variants from extraction_instruments() and their association statistics with ICH
  df_sysbp = merge_exposure_id_with_df_both_sexes_3(df_both_sexes_3, exposure_id, paste0('dataframes/', output_save_name, '.csv'))
  merged_df = harmonise_data(exposure_sysbp_snps, df_sysbp)
  
  #results csv
  results = mr(merged_df)
  results = generate_odds_ratios(results)
  write.csv(results, paste0('results/', output_save_name, '_results.csv'))
  
  #scatter plots
  p = mr_scatter_plot(results, merged_df)
  png(file=paste0('scatter_plots/', output_save_name, '_scatter.png'))
  print(p)
  
  #funnel plots
  funnel <- mr_singlesnp(merged_df)
  p2 <- mr_funnel_plot(funnel)
  png(file=paste0('funnel_plots/', output_save_name, '_funnel.png'))
  print(p2)
  
  #forest plots
  forest <- mr_singlesnp(merged_df)
  p3 <- mr_forest_plot(forest)
  png(file=paste0('forest_plots/', output_save_name, '_forest.png'))
  print(p3)
  
  #leave-one-out forest plots
  res_loo <- mr_leaveoneout(merged_df)
  p4 <- mr_leaveoneout_plot(res_loo)
  png(file=paste0('leave_one_out_forest_plots/', output_save_name, '_leave_one_out.png'))
  
  #MR_MoE (Machine Learning Approach)
  res <- mr_wrapper(merged_df)
  res_moe <- mr_moe(res, rf)
  expid = paste(exposure_id,sep = "", ".1e+20")
  write.csv(res_moe[[expid]][["estimates"]], paste0('mr_moe/', output_save_name, '_moe.csv'))
  
  #Contamination Mixture Method
  by <- merged_df$beta.outcome
  byse <- merged_df$se.outcome
  bx <- merged_df$beta.exposure
  bxse <- merged_df$se.exposure
  estimate_conmix <- contamination_mixture(bx,bxse,by,byse,exposure_id)
  #write.csv(estimate_conmix, paste0('con_mix/', output_save_name, '_conmix.csv'))
  
  #scatter plots con_mix
  
  
  
  
  dev.off()
}


##Contamination Mixture method
contamination_mixture = function(bx,bxse,by,byse, exposure_id) {
  iters = 2001; theta = seq(from=-1, to=1, by=2/(iters-1))
  # if the causal estimate (and confidence interval) is not expected to lie between -1 and 1
  # then change from and to (and maybe increase iters)
  ratio = by/bx; psi = 1.5*sd(ratio) # note: this is suggested as an initial value for psi
  # sensitivity to this parameter should be considered
  ratio.se = abs(byse/bx); # first-order
  #ratio.se = sqrt(byse^2/bx^2+by^2*bxse^2/bx^2) # second-order, assuming two-sample setting
  
  lik=NULL
  for (j1 in 1:iters) {
    lik.inc = exp(-(theta[j1]-ratio)^2/2/ratio.se^2)/sqrt(2*pi*ratio.se^2)
    lik.exc = exp(-ratio^2/2/(psi^2+ratio.se^2))/(sqrt(2*pi*(psi^2+ratio.se^2)))
    valid = (lik.inc>lik.exc)*1
    lik[j1] = prod(c(lik.inc[valid==1], lik.exc[valid==0]))
    if (which.max(lik)==length(lik)) { valid.best = valid }
  }
  phi = ifelse(sum(valid.best)<1.5, 1,
               max(sqrt(sum(((ratio[valid.best==1]-weighted.mean(ratio[valid.best==1],
                                                                 ratio.se[valid.best==1]^-2))^2*
                                                                 ratio.se[valid.best==1]^-2))/(sum(valid.best)-1)), 1))
  loglik = log(lik)
  
  whichin = which(2*loglik>(2*max(loglik)-qchisq(0.95, df=1)*phi^2))
  estimate = theta[which.max(loglik)]          # estimate
  lower_ci = theta[whichin[1]]                 # lower limit of CI
  upper_ci = theta[whichin[length(whichin)]]   # upper limit of CI
  cont_disj = sum(diff(whichin)!=1)             # this is 0 if the CI is a continuous range of values
  #  and 1+ if the CI consists of multiple disjoint ranges
  
  #calculating p-value 
  se_lower = (lower_ci-estimate)/-1.96
  se_upper = (upper_ci-estimate)/1.96
  se = (se_upper+se_lower)/2
  z_score = abs(estimate/se)
  f <- function(x) {1/sqrt(2*pi)*exp(-x^2/2)}
  p_value = integrate(f, lower = z_score, upper = 9e9999999)
  
  #creating dataframe
  df_conmix <- data.frame(exposure_id = c(exposure_id),
                          estimate = c(estimate),
                          lower_ci = c(lower_ci),
                          upper_ci = c(upper_ci),
                          se = c(se),
                          z_score = c(z_score),
                          p_value = c(p_value),
                          cont_disj = c(cont_disj))
  
  
  return(df_conmix)
}







