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
