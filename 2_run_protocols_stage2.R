#*****************************************************************************#
# Load inputs and set working directories: CHANGE HERE ****
#* [Put the full path of the input file]
cohort_input_file = "~/Library/CloudStorage/OneDrive-Personal/00_Work/00_NeuroCHARGE/00_2025/00_Obesity_BrainAtrophy/00_NeuroCHARGE_obesity_ctxVol-main_feb2026/1_cohort_specific_inputs_stage2.R"
#
#*****************************************************************************#Q

# 00 setup ----
## 1. read cohort specific inputs ----
source(cohort_input_file)

## 2. set working directory
setwd(opt$wd)

## 3. load packages and define custom functions ----
source('00_functions_library.R')
source('00_define_custom_functions.R')

## 4. create output directory and logfile ----
outdir00 = file.path(opt$wd,paste0("stage2_outputs_to_send_aug04_sobel_",opt$cohort_name,opt$input_metabo_platform))
if(opt$is.family.data){
  outdir00 = paste0(outdir00,"_FAM") 
}
dir.create(outdir00)

logfile00 = file.path(outdir00,str_replace(basename(cohort_input_file),".R",".log"))
file.copy(cohort_input_file, outdir00,overwrite = T)
file.rename(file.path(outdir00,basename(cohort_input_file)),logfile00)
session_info(to_file = logfile00)

# 01 prepare data ----
cat("Prep: sourcing functions\n")
source("01_prepare_data_stage2.R")

#--- [REMOVE LATER] 
df_all$waist = df_all$BMI
df_all$WHR = df_all$BMI
df_all$ancestry = sample(c("EUR","EAS"),size = nrow(df_all),replace = T)
table(df_all$ancestry)

# --- ANALYSES ---
#* [Pooled analyses of non-underweight individuals at midlife]
source("A1_pooled_mediation_test.R")

#* [APOE4-stratified analyses of non-underweight individuals at midlife]
source("A2_stratified_mediation_test.R")

#### WORK FROM HERE: THINK - each set of analyses can be performed parallely ####
#* [Age-group-stratified analyses of non-underweight individuals]
# source('04_Age_str_mediation_analysis.R')
str_mediation_test(
  data=df_all %>% filter(is.underweight=="No"),
  str.var = "Age",
  str.var.name =  "age.group",
  cov_cohort_specific=cov_cohort_specific,
  cov_additional = cov_additional,
  is.family.data = opt$is.family.data,
  logfile.name='04_Age_str_mediation_analysis.log'
  )

#* [BMI-stratified analyses of individuals at midlife]
# source('04_BMI_str_mediation_analysis.R')
str_mediation_test(
  data = df_all%>% filter(age.group=="(34,65]_(34,65]"),
  str.var = "BMI",
  str.var.name = "is.underweight",
  cov_cohort_specific=cov_cohort_specific,
  cov_additional = cov_additional,
  is.family.data = opt$is.family.data,
  logfile.name = '04_BMI_str_mediation_analysis.log'
)

#* [Ethnicity-stratified analyses of non-underweight individuals at midlife] # 
str_mediation_test(
  data = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"), 
  str.var = "Ethnicity",str.var.name = "ethnicity",
  cov_cohort_specific=cov_cohort_specific,
  cov_additional = cov_additional,
  is.family.data = opt$is.family.data,
  logfile.name = '04_ethnisity_str_mediation_analysis.log'
  )

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

#* [Pooled analyses of all individuals]
pooled_mediation_test(
  data=df_all,
  cov_cohort_specific = cov_cohort_specific, 
  cov_additional = cov_additional, 
  analysis.type = "primary_all",
  is.family.data = opt$is.family.data,
  logfile.name = '04_pooled_mediation_analysis_all.log',
  roi.34 = roi.34,
  group = 'pooled'
)

#* [Statin-use stratified analyses of non-underweight individuals at midlife] #
str_mediation_test(
  data = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"), 
  #* [NEED TO DEAL WITH THIS] *#
  str.var = "StatinUse",str.var.name = "statin_use",
  ##############################
  cov_cohort_specific=cov_cohort_specific,
  cov_additional = cov_additional,
  is.family.data = opt$is.family.data,
  logfile.name = "04_statin_str_mediation_analysis.log")
