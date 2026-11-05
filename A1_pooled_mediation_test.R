#* [Pooled analyses of non-underweight individuals at midlife]

#-------------------------------------------------------------------------@
# change this part ----
analysis.type = 'primary'# c("primary","primary_all",'adjEthnicity')
group = 'pooled'    
logfile.name = '04_pooled_non_underweight_at_midlife.log'
anal.age.group = "(34,65]_(34,65]"
anal.is.underweight = "No"

# does not change but need to define ----
data = df_all %>% #from previous step
  filter(is.underweight==anal.is.underweight,age.group==anal.age.group)
cov_cohort_specific=cov_cohort_specific
cov_additional = cov_additional
is.family.data = opt$is.family.data
roi.34 = roi.34
#-------------------------------------------------------------------------@
# starting ----
analysis_data = data
outdir = file.path(outdir00,paste0(group,"_",analysis.type))
dir.create(outdir)
logfile=file.path(outdir,logfile.name)
msg = paste0('\nMediation analysis starts here.\n'); message(msg);cat(msg,file=logfile)

##* identify columns with all missing or with SD=0 (no variation)
missingCol <- names(analysis_data)[colSums(is.na(analysis_data)) == nrow(analysis_data)]
col.SD = apply(analysis_data %>% select(where(is.numeric)),2,sd,na.rm=T)
mono = apply( analysis_data %>% select(where(is.character)), 2, function(x){ all(length(unique(na.omit(x))) == 1) } )
monoVar = attr(mono,"names")[mono]
rm.var = c(names(col.SD)[is.na(col.SD)|col.SD==0],monoVar,missingCol)
rm.var = unique(rm.var)
analysis_data <- analysis_data %>% select(-all_of(rm.var))
msg =paste0('\nRemoved ',length(rm.var),' variables due to 100% missing or no variations: ', paste0(rm.var,collapse = ", "),".","\n")
message(msg); cat(msg,file=logfile,append=T)

# 2. descriptive stats
if(str_detect(analysis.type, 'primary')){
  do.run = TRUE
}else{
  do.run = length(unique(analysis_data[['ethnicity']])) >1 #run if more than 1 ethnicity
}
if(do.run){
  source("02_get_desc_statistics_stage2.R",local = T)
  # check sample sizes for each categorical variable:
  small.sample_base <- check.sample.size(outdir,cov_cohort_specific)
  if(!is.null(small.sample_base$lessthan_5)){
    # discrete cohort-specific covs exist and can do sex-specific adjustment --> proceed to fitting
    run_base = !small.sample_base$lessthan_5
  }else{#is.null(small.sample_base$lessthan_5
    run_base = TRUE# discrete cohort-specific covs don't exist --> proceed to fitting
  }
  if(run_base){
    source('03_mediation_analysis_baseModel.R', local = T)#need to put is.family
  }else{
    msg = paste0('\nBefore 03: ',group," (n=",nrow(analysis_data),
                 '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                 'we will not run association tests. \n')
    message(msg)
    print(small.sample_base$analysis_data_minFreq)
    cat(msg,file=logfile,append = T)
    capture.output(print(small.sample_base$analysis_data_minFreq),file = logfile,append=T)
  }#else for if(run_base)
  small.sample_full <- check.sample.size(outdir,c(cohort_specific_covs,cov_additional))
  run_full = !small.sample_full$lessthan_5
  if(run_full){
    source('03_mediation_analysis_fullModel.R', local = T)
  }else{
    msg = paste0('\nBefore 3b: ',group," (n=",nrow(analysis_data),
                 '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                 'we will not run association tests. \n')
    message(msg)
    print(small.sample_full$analysis_data_minFreq)
    cat(msg,file=logfile,append = T)
    capture.output(print(small.sample_full$analysis_data_minFreq),file = logfile,append=T)
  }#else for if(run_full)
}else{#!do.run
  msg = paste0("\n There is only one ethnic group, \'",unique(analysis_data$ethnicity),"\' (n=",nrow(analysis_data),"): We will not run ethnicity-adjusted analyses.\n")
  message(msg)
  capture.output(msg,file=logfile,append=T)
}

# clean
rm(data,analysis_data,group,analysis.type)
capture.output(summary(warnings()),file=logfile,append=T)
