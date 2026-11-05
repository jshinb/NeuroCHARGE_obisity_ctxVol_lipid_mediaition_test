# --- Cohort study information ---

# Analyst/PI information for contact and study information: put double or single quotation
# marks
#
analysis_date = '2025-11-11'#*[YYYY-MM-DD]
analyst = 'Jean Shin'
analyst_email = 'jjshin.research@gmail.com'
PI = c('Zdenka Pausova','Tomas Paus')#*[use c('PI1 name','PI2 name',...) if more than 1]
PI_email =  c('zdenka.pausova@umontreal.ca','tomas.paus@umontreal.ca')#*[use c('PI1 email','PI2 email',...) if more than 1]
cohort_full_name = 'Saguenay Youth Study'#*[e.g., 'Saguenay Youth Study']
cohort_short_name = 'SYS_hypothesis_scaled_res1'#*[e.g., 'SYS_parents']
family_data = FALSE #*[e.g., TRUE if family-based data, otherwise (i.e., if unrelated) FALSE]
ethnicities = "White"#*["White" or use c('White', 'African American') in the presence of multiple ethnicities]
ancestries = "EUR"#*["EUR" or c('AFR', 'AMR', 'EAS', 'EUR', 'SAS',...) in the presence of multiple ancestries]
FreeSurfer_version = '5.3.0'#*[e.g., '5.3.0']
ICV_based_on =  "customized"#*[e.g., FreeSurfer/ customized pipeline (ref:PMID)]
ICV_pipeline_reference = 'PMIDXXX'#*['PMIDXXX']

# --- Data file information ---

# working directory(wd), file names for the brain and non-brain datasets, and missing value code
wd = "~/Library/CloudStorage/OneDrive-Personal/00_Work/00_NeuroCHARGE/00_2025/00_Obesity_BrainAtrophy/00_NeuroCHARGE_obesity_ctxVol-main_feb2026"#*[directory where the R scripts and data files are stored]
input_neuro = 'neuro_SPS.txt'#*[name of data file with brain image variables: e.g., "neuro_SPS.txt"]
input_non.neuro = 'non.neuro_hypothetical_SPS.txt'#*[name of the data file with non-brain variables: e.g., 'non.neuro_SPS.txt']
input_metabo = 'metabodata_adults.txt'#*[replace it with the name of data file with metabolomic data: e,g., 'metabo_Metabolon_SPS.txt']
input_metabo_cov = 'metabo_SYS.txt'#*[replace it with the name of data file with metabolomic data: e,g., 'metabo_Metabolon_SPS.txt']
input_metabo_platform = "Nightingale" #*[e.g., 'MetabolonHD', 'Nightingale', ...]
missing_value_code ='NA' #*[how missing values are coded]
is.family.data = FALSE #TRUE/FALSE if or not family data

# --- Columm Name Information  ---

cov_base = c('age_mri','age_adiposity','age_metabo','sex','ICV','fasting_hours') # - change if names are different
cov_additional = c('current_smoking','hypertension','T2D','statin','other_lipid_low_med')# - change if names are different
cov_genoPC = paste('genoPC',1:4,sep='')# - change this if different
cov_cohort_specific = c() # - change this if any

# --- Codings for categorical variables ---

sex_codes = c('M','F')#
current_smoking_codes = c('yes','no')
hypertension_codes = c('yes','no')
T2D_codes = c('yes','no')
statin_codes = c('yes','no') #on statin vs. not on statin
other_lipid_low_med_codes = c('yes','no') #on med or not on med

opt = list (
  # --- Cohort information ---
analysis_date=analysis_date,
analyst=analyst,
analyst_email=analyst_email,
PI=PI,
PI_email=PI_email,
cohort_full_name=cohort_full_name,
cohort_short_name=cohort_short_name,
family_data=family_data,
ethnicities=ethnicities,
ancestries=ancestries,
FreeSurfer_version=FreeSurfer_version,
ICV_based_on=ICV_based_on,
ICV_pipeline_reference=ICV_pipeline_reference,

# --- Data file information ---
wd=wd,
input_neuro=input_neuro,
input_non.neuro=input_non.neuro,
input_metabo=input_metabo,
input_metabo_cov=input_metabo_cov,
input_metabo_platform=input_metabo_platform,
missing_value_code =missing_value_code ,
is.family.data = is.family.data,

# --- Codings for categorical variables ---
# sex_codes=sex_codes,
# current_smoking_codes = current_smoking_codes,
# hypertension_codes=hypertension_codes,
# T2D_codes=T2D_codes,
# statin_codes=statin_codes,
# other_lipid_low_med_codes=other_lipid_low_med_codes,

# --- Columm Name Information  ---
cov_base=cov_base,
cov_additional=cov_additional,
cov_genoPC=cov_genoPC,
cov_cohort_specific=cov_cohort_specific

)
rm(list=setdiff(ls(), c("opt","cohort_input_file")))
