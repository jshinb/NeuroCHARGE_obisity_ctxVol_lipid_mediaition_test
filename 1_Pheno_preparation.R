#! /usr/bin/Rscript

### Script written by  Jean Shin <jjshin.research@gmail.com> (adopted scripts written by Anna Furtjes)
### Distributed as part of our lifetime brain atrophy GWAS consortium effort
### Last edited Oct 2025

#*****************************************************************************#
# Install packages ----
#*****************************************************************************#
source('00_install_libs.R')


#*****************************************************************************#
# source functions
#*****************************************************************************#
if(!file.exists(file.path(opt$dir, "Rfunctions.R"))){
  stop("Please make sure you save the file called 'Rfunctions.R' in the same directory as 'pheno_preparation.R' as well as your neuroimaging and covariate input files. Currently it is not in ", opt$wd,".")
}
source(file.path(opt$dir, "Rfunctions.R"))

#*****************************************************************************#
# Record output ----
#*****************************************************************************#
# create directory to hold results
print("Starting analysis.")
print(Sys.time())
time1=Sys.time()

#*****************************************************************************#
# Sanity checks for input files ----
#*****************************************************************************#
#* Check if files exist ----
if(!file.exists(file.path(opt$wd,opt$input_neuro))){
  stop(paste0("The file you indicated with --input_neuro does not exist: ", opt$wd, opt$input_neuro))
}
if(!file.exists(file.path(opt$wd,opt$input_non.neuro))){
  stop(paste0("The file you indicated with --input_non.neuro does not exist: ", opt$wd, opt$input_non.neuro))
}

message(
  paste0("Reading in files:\n1. ", file.path(opt$wd,opt$input_neuro), 
         "&\n2. ", file.path(opt$wd,opt$input_non.neuro))
)

#* Read in files ----
neuro_all = fread(file.path(opt$wd,opt$input_neuro))
non.neuro_all = fread(file.path(opt$wd,opt$input_non.neuro))

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
  stop(paste0("Your --input_neuro_all file appears to have incorrect column names, or to miss certain columns.", "Check documentation for correct and case-sensitive naming and re-run this function."))
}

if(!all(non.neuro_required_colnames %in% names(non.neuro_all))){
  stop(paste0("Your --input_non.neuro_all file appears to have incorrect column names, or to miss certain columns. Your column names are: ", paste0(names(non.neuro_all), collapse = ", "),". Check documentation for correct and case-sensitive naming and re-run this function. It should be (something like) ID, age_mri,age_adiposity, sex, ethnicity, genoPC1, genoPC2, genoPC3, genoPC4, [other study-specific covariates]."))
}

ind_genoPC=str_detect(names(non.neuro_all),"genoPC")
if(!any(ind_genoPC)){
  stop(
    paste0(
      "Your --input_non.neuro_all file appears to miss. Please provied >=4 leading genetic PCs.",
      "Check documentation for correct and case-sensitive naming and re-run this function. It should be genoPC1, genoPC2, genoPC3, genoPC4,...."))
}

#* Check if NA is interpreted correctly ----
neuro_all[neuro_all == opt$NAsymbol] <- NA
non.neuro_all[non.neuro_all == opt$NAsymbol] <- NA

#*****************************************************************************#
# Format brain image and other types of data ----
#*****************************************************************************#
#* Format brain image data (i.e., cortical structure and ICV) ----
#** Make sure all MRI variables are numeric and set ID to be character ----
neurovars <- setdiff(neuro_all_colnames,"ID")
if(!all(unlist( lapply( subset(neuro_all,select=neurovars), is.numeric ) ))){
  stop("All neuro-imaging variables must be numeric: Please check your data.\n")
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
  
cov_genoPC = names(non.neuro_all)[str_detect(names(non.neuro_all),"genoPC")]
cov_additional = c('current_smoking','hypertension','T2D')
cov_cohort_specific = setdiff(names(non.neuro_all),
                              c(non.neuro_required_colnames,cov_genoPC))

#* Format other types of data ----
#** ID and FID are characters 
non.neuro_all = non.neuro_all %>% mutate(ID = as.character(ID), FID = as.character(FID))
non.neuro = non.neuro_all %>% 
  mutate(is.underweight = ifelse(BMI<20, "Yes","No")) %>% 
  mutate(age_diff=age_mri-age_adiposity) %>%
  mutate(age.group_mri = cut(age_mri,breaks = c(34,65), right=TRUE) ) %>%
  mutate(age.group_adiposity = cut(age_adiposity,breaks = c(34,65), right=TRUE)) %>%
  mutate(age.group = paste(age.group_mri,age.group_adiposity,sep="_")) %>%
  mutate(APOE4.present = ifelse(APOE %in% c("e2e2","e2e3",'e3e3'),'No','Yes')) %>%
  mutate(APOE4.present = ifelse(APOE == 'e2e4',NA, APOE4.present))

#* Merge neuro and non-neuro datasets ----
df_all = non.neuro %>% inner_join(neuro, join_by(ID)) 

cat("\n **** Please proceed to the next steps. **** \n",sep='')
