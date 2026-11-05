# To-think
1. Make presentation files
2. Covariates for metabolites
3. Models: check!!
4. Which sets of analyses should we run? Potentially, it can be computationally expensive (e.g., for Rhineland) - at the moment analyses are performed sequentially

2026-11-05
* The previous versions - I used functions to run each set of analyses; however, it was hard to debug the code when errors occurred, so I have decided to create a separate script for each analysis and source it.
Examples:

#* [Ancestry-stratified analyses of non-underweight individuals at midlife] # 
str_mediation_test(
  data = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"), 
  str.var = "Ancestry",str.var.name = "ancestry",
  cov_cohort_specific=cov_cohort_specific,
  cov_additional = cov_additional,
  is.family.data = opt$is.family.data,
  logfile.name = "04_ancestry_str_mediation_analysis.log")

#* [Pooled analyses of non-underweight individuals at midlife adjusting for ethnicity] # 
pooled_mediation_test(
  data = df_all %>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"),
  analysis.type = 'primary',# c("primary","primary_all",'adjEthnicity')
  cov_cohort_specific = cov_cohort_specific, 
  cov_additional = cov_additional, 
  logfile.name = '04_pooled_mediation_analysis_non_underweight_midlife.log',
  is.family.data = opt$is.family.data,
  roi.34 = roi.34,
  group = 'pooled')
