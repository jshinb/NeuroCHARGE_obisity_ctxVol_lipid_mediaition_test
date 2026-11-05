#! /usr/bin/Rscript

#* Script written by Jean Shin <jjshin.research@gmail.com> (adopted scripts written by Anna Furtjes)
#* Created: Jan 14, 2026
#* Last modifioed: Jan 14, 2026
#* Distributed for the NeuroCHARGE project: Midlife obesity and cortical Atrophy

#*****************************************************************************#
# Record output ----
#*****************************************************************************#
# create directory to hold results
print("Starting analysis.")
print(Sys.time())
time1=Sys.time()

starting.message=paste0("\nData preparation starts here.\n\n")
message(starting.message);cat(starting.message,file=logfile00, append=T)

#*****************************************************************************#
# Sanity checks for input files ----
#*****************************************************************************#
#* Check if files exist ----
if(!file.exists(file.path(opt$wd,opt$input_neuro))){
  stop.msg = paste0("The file you indicated with --input_neuro does not exist:\n", opt$wd, opt$input_neuro,"\n\n")
  cat(stop.msg, file=logfile00, append=T); stop(stop.msg)
}
if(!file.exists(file.path(opt$wd,opt$input_non.neuro))){
  stop.msg = paste0("The file you indicated with --input_non.neuro does not exist: ", opt$wd, opt$input_non.neuro,"\n\n") 
  cat(stop.msg, file=logfile00, append=T); stop(stop.msg)
}
out.msg = paste0("Reading in files:\n1. ", file.path(opt$wd,opt$input_neuro), 
                 "\n2. ", file.path(opt$wd,opt$input_non.neuro),
                 "\n3. ", file.path(opt$wd,opt$input_metabo),"\n\n")
message(out.msg); cat(out.msg,file=logfile00,append=T); rm(out.msg)

#* Read in files ----
neuro_all = fread(file.path(opt$wd,opt$input_neuro))
non.neuro_all = fread(file.path(opt$wd,opt$input_non.neuro))
metabo_all = fread(file.path(opt$wd,opt$input_metabo))
metabo_cov_all = fread(file.path(opt$input_metabo_cov))

#* Examine if column names are correct ----
neuro_all_colnames = fread('neuro_columns.txt',header=F) %>% pull(V1)
non.neuro_required_colnames = fread('non.neuro_columns.txt',header=F) %>% pull(V1)
ind_genoPC = str_detect(names(non.neuro_all),"genoPC")
genoPCs = names(non.neuro_all)[ind_genoPC]
cohort_specific_covs = setdiff(
  names(non.neuro_all),
  c(non.neuro_required_colnames,genoPCs)
  )

if(!all(neuro_all_colnames %in% names(neuro_all))){
  stop.msg =paste0("Your --input_neuro file appears to have incorrect column names, or to miss certain columns.", 
                   "Check documentation for correct and case-sensitive naming and re-run this function.",
                   "\n\n")
  cat(stop.msg,file=logfile00,append=T);stop(stop.msg)
}

if(!all(non.neuro_required_colnames %in% names(non.neuro_all))){
  stop.msg = paste0("Your --input_non.neuro file appears to have incorrect column names, or to miss certain columns. ",
                    "Your column names are: ", paste0(names(non.neuro_all), collapse = ", "),
                    ". Check documentation for correct and case-sensitive naming and re-run this function. ",
                    "It should be (something like) ID, age_mri,age_adiposity, sex, ethnicity, genoPC1, genoPC2, genoPC3, genoPC4, [other study-specific covariates].",
                    "\n\n")
  cat(stop.msg,file=logfile00,append=T);stop(stop.msg)
}

if(!any(names(metabo_all)=="ID")){
  stop.msg = paste0("Your --input_metabo file appears to have 'ID' column in a different name. Your ID column name must be \'ID\'\n. ", 
              "Check the name of your participant ID column and rerun the script.")
  cat(stop.msg,file=logfile00,append=T);stop(stop.msg)
}

ind_genoPC=str_detect(names(non.neuro_all),"genoPC")
if(!any(ind_genoPC)){
  stop.msg =paste0(
    "Your --input_non.neuro file appears to miss genotype PCs. Please provied >=4 leading genetic PCs.",
    "Check documentation for correct and case-sensitive naming and re-run this function. It should be genoPC1, genoPC2, genoPC3, genoPC4,....",
    "\n\n")
  cat(stop.msg,file=logfile00,append=T);stop(stop.msg)
}

#* Check if NA is interpreted correctly ----
neuro_all[neuro_all == opt$NAsymbol] <- NA
non.neuro_all[non.neuro_all == opt$NAsymbol] <- NA
metabo_all[metabo_all == opt$NAsymbol] <- NA
metabo_cov_all[metabo_cov_all == opt$NAsymbol] <- NA

#*****************************************************************************#
# Format brain image and other types of data ----
#*****************************************************************************#
#* Format brain image data (i.e., cortical structure and ICV) ----
#** Make sure all MRI variables are numeric and set ID to be character ----
roi.34 = fread('freesurfer_34rois.txt', header=F) %>% pull(V1)

neurovars <- setdiff(neuro_all_colnames,"ID")
if(!all(unlist( lapply( subset(neuro_all,select=neurovars), is.numeric ) ))){
  stop.msg = paste0("All neuro-imaging variables must be numeric: Please check your data.\n\n")
  cat(stop.msg,file=logfile00,append=T);stop(stop.msg)
}

#** Calculate values across right and left hemispheres ---- 
neuro_all = neuro_all %>% mutate(ID = as.character(ID))
#* calculate the regional and total cortical volumes and format 'neuro' to 
#* include the 34 regional-, total-volumes
#* include the 34 regional-, mean-thickness
#* include the 34 regional-, total-surface area
#* ICV
neuro_volume = get_sums_volume(neuro_all,ID_col = "ID")
neuro_area = get_sums_area(neuro_all,ID_col = "ID")
neuro_thickness = get_means_thickness(neuro_all,ID_col = "ID")
neuro = neuro_volume %>%
  left_join(neuro_thickness, join_by(ID)) %>%
  left_join(neuro_area, join_by(ID)) %>%
  left_join(neuro_all %>% select(ID,ICV), join_by(ID))

#* Format adiposity and covariate data ----
#** ID and FID are characters 
non.neuro_all = non.neuro_all %>% mutate(ID = as.character(ID), FID = as.character(FID))
non.neuro = non.neuro_all %>% 
  left_join(metabo_cov_all,join_by(ID)) %>% #merge with metabo covariates
  mutate(age_diff2 = age_mri - age_metabo) %>%
  mutate(is.underweight = ifelse(BMI<20, "Yes","No")) %>% 
  mutate(age_diff=age_mri-age_adiposity) %>%
  mutate(age.group_mri = cut(age_mri,breaks = c(34,65), right=TRUE) ) %>%
  mutate(age.group_adiposity = cut(age_adiposity,breaks = c(34,65), right=TRUE)) %>%
  mutate(age.group = paste(age.group_mri,age.group_adiposity,sep="_")) %>%
  mutate(APOE4.present = ifelse(APOE %in% c("e2e2","e2e3",'e3e3'),'No','Yes')) %>%
  mutate(APOE4.present = ifelse(APOE == 'e2e4',NA, APOE4.present))

#* Format metabolomic data ----
#* Note: Assuming the metabolomic data has been cleaned (e.g., outliers, missing rates)
# read in data
metabo_all = metabo_all %>% mutate(ID = as.character(ID))
metaboID_ordered = names(metabo_all)[!names(metabo_all) %in% c("ID",'FID')]

#do.I.need <- function(# do I need to merge? - Yes for using the previous script as much as possible
# Assumption - data are cleaned
cov_metabo = c('age_metabo','statin_use','othermed_use','eGFR')#required
cov_genoPC = names(non.neuro_all)[str_detect(names(non.neuro_all),"genoPC")]
cov_additional = c('current_smoking','hypertension','T2D','eGFR')##
cov_cohort_specific = setdiff(names(non.neuro_all),
                              c(non.neuro_required_colnames,cov_genoPC,cov_metabo))

#1. Merge all types of datasets ----
df_all <- non.neuro %>% 
  inner_join(neuro, join_by(ID)) %>%
  inner_join(metabo_all, join_by(ID))

out.msg = paste0("\nData preparation is complete. Final sample size: N = ", nrow(df_all),'.\n')
message(out.msg); cat(out.msg,file=logfile00,append=T); rm(out.msg)

cat("\n **** Please proceed to the next steps. **** \n")
