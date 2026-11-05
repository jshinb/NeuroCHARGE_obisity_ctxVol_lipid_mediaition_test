#* [APOE4-stratified analyses of non-underweight individuals at midlife]

#-------------------------------------------------------------------------@
# change this part ----
analysis.type = 'primary'# c("primary","primary_all",'adjEthnicity')
str.var = "APOE4"
str.var.name = "APOE4.present"
anal.age.group = "(34,65]_(34,65]"
anal.is.underweight = "No"
logfile.name='04_APOE4_str_mediation_analysis.log'

# does not change but need to define ----
data = df_all %>% #from previous step
  filter(is.underweight==anal.is.underweight,age.group==anal.age.group)
cov_cohort_specific=cov_cohort_specific
cov_additional = cov_additional
is.family.data = opt$is.family.data

#-------------------------------------------------------------------------@
# starting ----
groups = unique(na.omit(data[[str.var.name]]))

roi.34 = fread("freesurfer_34rois.txt",header=F) %>% pull(V1)
if(length(groups)>1){
  # 2. descriptive stats
  for(group in groups){
    outdir = file.path(outdir00,paste0(str.var,"_",str.var.name,"_",group))
    dir.create(outdir)
    logfile=file.path(outdir,logfile.name)
    msg = paste0('\nMediation analysis starts here.\n'); message(msg);cat(msg,file=logfile)
    analysis_data = data %>% filter(.data[[str.var.name]] == group)
    
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
    
    source("02_get_desc_statistics_stage2.R",local = T)
    # check sample sizes for each categorical variable:
    small.sample_base <- check.sample.size(outdir,cov_cohort_specific)
    if(!is.null(small.sample_base$lessthan_5)){
      # discrete cohort-specific covs exist and can do sex-specific adjustment --> proceed to fitting
      run_base = !small.sample_base$lessthan_5
    }else{#is.null(small.sample_base$lessthan_5
      run_base = TRUE# discrete cohort-specific covs don't exist --> proceed to fitting
    }
    if(run_base) {
      
      source("03_mediation_analysis_baseModel.R",local = T)
      
    }else{
      msg = paste0('\nBefore 3a: ',str.var,"_",group," (n=",nrow(analysis_data),
                   '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                   'we will not run association tests. \n')
      message(msg)
      print(small.sample_base$df_minFreq)
      cat(msg,file=logfile,append = T)
      capture.output(print(small.sample_base$df_minFreq),file = logfile,append=T)
    }#else for if(run_base)
    small.sample_full <- check.sample.size(outdir,c(cohort_specific_covs,cov_additional))
    run_full = !small.sample_full$lessthan_5
    if(run_full){
      source("03_mediation_analysis_fullModel.R",local = T)
    }else{
      msg = paste0('\nBefore 3b: ',str.var,"_",group," (n=",nrow(analysis_data),
                   '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                   'we will not run association tests. \n')
      message(msg)
      print(small.sample_full$df_minFreq)
      cat(msg,file=logfile,append = T)
      capture.output(print(small.sample_full$df_minFreq),file = logfile,append=T)
    }#else for if(run_full)
  }#for(group in groups)
}else{#length(groups)=1 or 0
  group = groups[1]
  outdir = file.path(outdir00,paste0(str.var,"_",str.var.name,"_",group))
  dir.create(outdir)
  logfile=file.path(outdir,logfile.name)
  analysis_data = data
  msg = paste0("\n There is only one group:", str.var,"_",group," (n=",nrow(analysis_data),"): No need to run ",str.var,"-stratified analyses.\n")
  message(msg)
  cat(msg,file=logfile,append = T)
}
# clean
rm(data,analysis_data,group)
