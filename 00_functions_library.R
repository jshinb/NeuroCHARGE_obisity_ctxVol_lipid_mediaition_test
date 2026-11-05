# *****************************************************************************
# Consortium Mediation: Library of Functions
# Script written by Anna Furtjes (for Lifetime Brain Atrophy GWAS)
# Adopted/modified by: Jean Shin, jjshin.research@gmail.com for any questions 
# Edited on: 2026-02-06 
# Version 2.1 - Optimized for Parallel Execution and Technical Invariance
# *****************************************************************************

# --- 0. Install and load libraries ---

#!/usr/bin/Rscript
library(psych)
library(lme4)

if(!is.element("pacman",installed.packages()[,1])){
  install.packages("pacman")
}

pacman::p_load(#check if a package is installed, if not, it attempts to install it
  sessioninfo,  #
  tidyverse,    #readr, dplyr, ggplot2, purrr, tidyr
  mgcv,         # for adjusting continuous covariates using non-parametric smooth models
  lme4,lmerTest,# for fitting regression to family-data
  EnvStats,     # for outliers [grubs]
  mediation,    # for mediation analyses (stage 2)
  pbapply,      # for showing progressing for mediation analyses
  future.apply, # for parallel computation
  optparse, psych, data.table, tableone, reshape2, # data manipulation
  ggpubr, cowplot, patchwork, GGally, corrplot)    # for ggplots


# manual installation of hrbrtheme package 
if(!is.element("hrbrthemes",installed.packages()[,1])){
  remotes::install_github("hrbrmstr/hrbrthemes")
};library(hrbrthemes)

# re-define some functions:
select = dplyr::select
describe = psych::describe
mediate = mediation::mediate

# --- 1. Basic Statistical Helpers ---

#' Inverse normal transformation (INT)
#'
#' https://www.biostars.org/p/80597/
#' See the supplement of Yang et al. Nature 2012. 
#'
#' @example
#' x1 <- 5:1
#' inormal(x1)
#'
#' x2 <- c(NA, 5:1, 5:1, NA)
#' inormal(x2)
inormal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

#' Outlier removal (4 Standard Deviations)
rm.sd4 = function(x){
  x.scaled = scale(x)[,1]
  ind = !is.na(x.scaled) & abs(x.scaled) > 4
  x[ind] <- NA
  x
}

# --- 2. Intermediate functions ---
get_sums_volume <- function(data,ID_col){
  #data -> cortical structure data 
  rois = fread('freesurfer_34rois.txt',header=F)$V1
  totalVol = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"volume",sep="_")
    Rh = paste("rh",roi,"volume",sep="_")
    totalVol_roi = (data[[Lh]] + data[[Rh]])
    totalVol = cbind(totalVol,totalVol_roi) #34 rois
  }
  L_Vol = 'lh_CortexVol_volume'
  R_Vol = 'rh_CortexVol_volume'
  global_totalVol = (data[[L_Vol]] + data[[R_Vol]])/2
  totalVol = data.frame(totalVol, global_totalVol)
  names(totalVol) = c(ID_col,paste(c(rois,"total"),"volume",sep="_"))
  # names(totalVol) = c(ID_col,rois,"total")#ID, roi-volume, and total volume
  totalVol
}

get_sums_area = function(data,ID_col){
  #data -> cortical structure
  rois = fread('freesurfer_34rois.txt',header=F)$V1
  total_area = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"area",sep="_")
    Rh = paste("rh",roi,"area",sep="_")
    total_area_roi = (data[[Lh]] + data[[Rh]])
    total_area = cbind(total_area,total_area_roi) #34 rois
  }
  L_area = 'lh_WhiteSurfArea_area'
  R_area = 'rh_WhiteSurfArea_area'
  global_total_area = (data[[L_area]] + data[[R_area]])/2
  total_area = data.frame(total_area,global_total_area)
  names(total_area) = c(ID_col,paste(c(rois,"total"),"area",sep="_")) #ID, roi-area, and total area
  total_area
}

get_means_thickness <- function(data,ID_col){
  #data -> cortical thickness data 
  rois = fread('freesurfer_34rois.txt',header=F)$V1
  avgTH = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"thickness",sep="_")
    Rh = paste("rh",roi,"thickness",sep="_")
    avgTH_roi = (data[[Lh]] + data[[Rh]])/2
    avgTH = cbind(avgTH,avgTH_roi) #34 rois
  }
  L_TH = 'lh_MeanThickness_thickness'
  R_TH = 'rh_MeanThickness_thickness'
  global_avgTH = (data[[L_TH]] + data[[R_TH]])/2
  avgTH = data.frame(avgTH, global_avgTH)
  names(avgTH) = c(ID_col,paste(c(rois,"mean"),"thickness",sep="_"))
  avgTH
}

# --- 3. Data exploration: plots ---
plot_hist <- function(dat = neuro, var, col = "grey"){
  
  # make sure input data is data.frame
  dat = as.data.frame(dat)
  # rename for simplicity
  dat$var = dat[,var]
  
  # calculate mean from original variable
  varOr = str_remove(var, "_stand")
  mean = mean(dat[,varOr], na.rm=T)
  sd = sd(dat[,varOr], na.rm=T)
  
  # get title
  if(var == "total"){
    axisName = "total cortical volume"}else
    {
      axisName=var
    }
  
  # binwidth should differ dependeding on sample size
  if(nrow(dat) < 500){width = 0.5}
  if(nrow(dat) >=500){width = 0.1}
  
  # PLOT
  # different output when there is a "sample" column
  plot = ggplot(dat, aes(x = var))+
    geom_histogram(alpha = 0.5, fill = col) + #, binwidth = width)+
    xlab(axisName)+
    ylab("Count")+
    theme_bw()+
    #scale_x_continuous(sec.axis = sec_axis(trans=~.*sd+mean))+
    theme(panel.border = element_blank())
  
  return(plot)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

plot_heatmap = function(dat, axisNames)
{
  require(ggplot2)
  # define function to obtain lower triangle
  # get correlation matrix for both samples together
  cor = cor(dat, use = "p")
  cor = get_lower_tri(cor)
  # melt matrix
  melted = reshape2::melt(cor)
  # get rounded value
  melted$value_round = round(melted$value, digit = 2)
  melted$distance0 = abs(melted$value)
  
  # plot
  p = ggplot(data = melted)+
    geom_point(aes(x = Var1, y = Var2, shape = value, fill = value, size = distance0), 
               shape = 21, alpha = 0.7, colour = "white") +
    scale_fill_gradient2(low = "#82A0D8", high = "#8DDFCB", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab" ,name="Correlation", guide = "legend")+
    scale_size_continuous(range = c(1, 15), guide = "none")+
    geom_text(aes(Var1, Var2, label = value_round), color = "black", size = 4)+
    xlab("")+
    ylab("")+
    scale_x_discrete(labels = axisNames)+
    scale_y_discrete(labels = axisNames)+
    guides(fill = "none")+
    theme_bw()+
    theme(panel.border = element_blank(),
          axis.text.x = element_text(angle=-45,vjust = 0.5, hjust=0))
  
  return(p)
}

generate_histograms <- function(outdir, metabo_OUT, metabo_info, metabo_columns) {
  dir.create(file.path(outdir, "plots"), showWarnings = FALSE)
  
  # Melt data for facetting
  d_long <- metabo_OUT %>% 
    select(all_of(metabo_columns)) %>% 
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  # Generate plots in batches of 20 for readability
  vars <- unique(d_long$variable)
  batches <- split(vars, ceiling(seq_along(vars)/20))
  
  for(i in seq_along(batches)) {
    p <- ggplot(d_long %>% filter(variable %in% batches[[i]]), aes(x = value)) +
      geom_histogram(fill = "lightsteelblue2", bins = 30) +
      facet_wrap(~variable, scales = "free", ncol = 5) +
      theme_bw() +
      labs(title = paste("Metabolite Distributions - Batch", i),
           caption = "Note: Red dashed lines indicate approximate 4SD boundaries") +
      geom_vline(xintercept = c(-4, 4), linetype = "dashed", color = "red", alpha = 0.5)
    
    ggsave(file.path(outdir, "plots", paste0("metabo_dist_batch_", i, ".pdf")), p, width = 12, height = 8)
  }
}

# --- 4. Brain ROI Adjustment Functions ---

# Adjustment for Volume (Adjusts for ICV)
get_adj_ctxV = function(roii, d, covariate_names, logfile){
  yname = paste(roii, "volume", sep="_")
  d = d %>% mutate(age.c = scale(age_mri)[,1], age.c2 = (scale(age_mri)/2)^2,
                   icv.c = scale(ICV)[,1], icv.c2 = (scale(ICV)/2)^2)
  
  analdati = d %>% select(ID, sex, all_of(covariate_names), age.c, age.c2, icv.c, icv.c2, !!sym(yname))
  analdati[['y']] = inormal(analdati[[yname]])
  analdati[['y']] = rm.sd4(analdati[['y']])
  
  mod_rhs = paste0(c('age.c','age.c2','icv.c','icv.c2', setdiff(covariate_names, c('age_mri',"ICV",'sex'))), collapse = " + ")
  formula_str = paste("y ~", mod_rhs)
  
  fitM = lm(as.formula(formula_str), data=analdati %>% filter(sex=="M"), na.action=na.exclude)
  fitF = lm(as.formula(formula_str), data=analdati %>% filter(sex=="F"), na.action=na.exclude)
  
  resM = analdati %>% filter(sex=="M") %>% mutate(!!paste0(roii,"_volume.adjSex") := resid(fitM)) %>% select(ID, last_col())
  resF = analdati %>% filter(sex=="F") %>% mutate(!!paste0(roii,"_volume.adjSex") := resid(fitF)) %>% select(ID, last_col())
  
  capture.output(paste0("get_adj_ctxV [", roii, "]: ", formula_str), file=logfile, append=T)
  return(bind_rows(resM, resF))
}

# Adjustment for Thickness (Typically NOT adjusted for ICV)
get_adj_ctxTH = function(roii, d, covariate_names, logfile){
  yname = paste(roii, "thickness", sep="_")
  d = d %>% mutate(age.c = scale(age_mri)[,1], age.c2 = (scale(age_mri)/2)^2)
  
  analdati = d %>% select(ID, sex, all_of(covariate_names), age.c, age.c2, !!sym(yname))
  analdati[['y']] = inormal(analdati[[yname]])
  analdati[['y']] = rm.sd4(analdati[['y']])
  
  mod_rhs = paste0(c('age.c','age.c2', setdiff(covariate_names, c('age_mri','sex'))), collapse = " + ")
  formula_str = paste("y ~", mod_rhs)
  
  fitM = lm(as.formula(formula_str), data=analdati %>% filter(sex=="M"), na.action=na.exclude)
  fitF = lm(as.formula(formula_str), data=analdati %>% filter(sex=="F"), na.action=na.exclude)
  
  resM = analdati %>% filter(sex=="M") %>% mutate(!!paste0(roii,"_thickness.adjSex") := resid(fitM)) %>% select(ID, last_col())
  resF = analdati %>% filter(sex=="F") %>% mutate(!!paste0(roii,"_thickness.adjSex") := resid(fitF)) %>% select(ID, last_col())
  
  capture.output(paste0("get_adj_ctxTH [", roii, "]: ", formula_str), file=logfile, append=T)
  return(bind_rows(resM, resF))
}

# Adjustment for Surface Area (Adjusts for ICV)
get_adj_ctxSA = function(roii, d, covariate_names, logfile){
  yname = paste(roii, "area", sep="_")
  d = d %>% mutate(age.c = scale(age_mri)[,1], age.c2 = (scale(age_mri)/2)^2,
                   icv.c = scale(ICV)[,1], icv.c2 = (scale(ICV)/2)^2)
  
  analdati = d %>% select(ID, sex, all_of(covariate_names), age.c, age.c2, icv.c, icv.c2, !!sym(yname))
  analdati[['y']] = inormal(analdati[[yname]])
  analdati[['y']] = rm.sd4(analdati[['y']])
  
  mod_rhs = paste0(c('age.c','age.c2','icv.c','icv.c2', setdiff(covariate_names, c('age_mri',"ICV",'sex'))), collapse = " + ")
  formula_str = paste("y ~", mod_rhs)
  
  fitM = lm(as.formula(formula_str), data=analdati %>% filter(sex=="M"), na.action=na.exclude)
  fitF = lm(as.formula(formula_str), data=analdati %>% filter(sex=="F"), na.action=na.exclude)
  
  resM = analdati %>% filter(sex=="M") %>% mutate(!!paste0(roii,"_area.adjSex") := resid(fitM)) %>% select(ID, last_col())
  resF = analdati %>% filter(sex=="F") %>% mutate(!!paste0(roii,"_area.adjSex") := resid(fitF)) %>% select(ID, last_col())
  
  capture.output(paste0("get_adj_ctxSA [", roii, "]: ", formula_str), file=logfile, append=T)
  return(bind_rows(resM, resF))
}

# --- Adjustment for BMI ---
get_adj_BMI = function(d, covariate_names, logfile){
  d = d %>% mutate(age.c = scale(age_adiposity)[,1], age.c2 = (scale(age_adiposity)/2)^2)
  
  analdati = d %>% select(ID, sex, all_of(covariate_names), age.c, age.c2, BMI)
  analdati[['y']] = inormal(analdati[['BMI']])
  analdati[['y']] = rm.sd4(analdati[['y']])
  
  mod_rhs = paste0(c('age.c','age.c2', setdiff(covariate_names, c('age_adiposity','sex'))), collapse = " + ")
  formula_str = paste("y ~", mod_rhs)
  
  fitM = lm(as.formula(formula_str), data=analdati %>% filter(sex=="M"), na.action=na.exclude)
  fitF = lm(as.formula(formula_str), data=analdati %>% filter(sex=="F"), na.action=na.exclude)
  
  resM = analdati %>% filter(sex=="M") %>% mutate(BMI.adjSex = resid(fitM)) %>% select(ID, BMI.adjSex)
  resF = analdati %>% filter(sex=="F") %>% mutate(BMI.adjSex = resid(fitF)) %>% select(ID, BMI.adjSex)
  
  capture.output(paste0("get_adj_BMI: ", formula_str), file=logfile, append=T)
  return(bind_rows(resM, resF))
}

# --- 5. Adjustment for Waist Circumference (WC) ---

# Note: Usually adjusted for Height to distinguish central adiposity from overall size
get_adj_WC = function(d, covariate_names, logfile){
  d = d %>% mutate(age.c = scale(age_adiposity)[,1], age.c2 = (scale(age_adiposity)/2)^2,
                   height.c = scale(height)[,1])
  
  analdati = d %>% select(ID, sex, all_of(covariate_names), age.c, age.c2, height.c, waist)
  analdati[['y']] = inormal(analdati[['waist']])
  analdati[['y']] = rm.sd4(analdati[['y']])
  
  mod_rhs = paste0(c('age.c','age.c2','height.c', setdiff(covariate_names, c('age_adiposity','sex','height'))), collapse = " + ")
  formula_str = paste("y ~", mod_rhs)
  
  fitM = lm(as.formula(formula_str), data=analdati %>% filter(sex=="M"), na.action=na.exclude)
  fitF = lm(as.formula(formula_str), data=analdati %>% filter(sex=="F"), na.action=na.exclude)
  
  resM = analdati %>% filter(sex=="M") %>% mutate(waist.adjSex = resid(fitM)) %>% select(ID, waist.adjSex)
  resF = analdati %>% filter(sex=="F") %>% mutate(waist.adjSex = resid(fitF)) %>% select(ID, waist.adjSex)
  
  capture.output(paste0("get_adj_WC: ", formula_str), file=logfile, append=T)
  return(bind_rows(resM, resF))
}

# --- 6. Adjustment for Waist-to-Hip Ratio (WHR) ---
get_adj_WHR = function(d, covariate_names, logfile){
  d = d %>% mutate(age.c = scale(age_adiposity)[,1], age.c2 = (scale(age_adiposity)/2)^2)
  
  analdati = d %>% select(ID, sex, all_of(covariate_names), age.c, age.c2, WHR)
  analdati[['y']] = inormal(analdati[['WHR']])
  analdati[['y']] = rm.sd4(analdati[['y']])
  
  mod_rhs = paste0(c('age.c','age.c2', setdiff(covariate_names, c('age_adiposity','sex'))), collapse = " + ")
  formula_str = paste("y ~", mod_rhs)
  
  fitM = lm(as.formula(formula_str), data=analdati %>% filter(sex=="M"), na.action=na.exclude)
  fitF = lm(as.formula(formula_str), data=analdati %>% filter(sex=="F"), na.action=na.exclude)
  
  resM = analdati %>% filter(sex=="M") %>% mutate(WHR.adjSex = resid(fitM)) %>% select(ID, WHR.adjSex)
  resF = analdati %>% filter(sex=="F") %>% mutate(WHR.adjSex = resid(fitF)) %>% select(ID, WHR.adjSex)
  
  capture.output(paste0("get_adj_WHR: ", formula_str), file=logfile, append=T)
  return(bind_rows(resM, resF))
}
# --- 7. Adjustment for Metabolite ---
adjust_metabolites <- function(data, lipids, age_col, sex_col, other_covs=NULL, logfile) {
  message("Adjusting metabolites: Sex-stratified residuals with 4SD cleaning.\n")
  
  adjusted_data <- data
  lipids_to_proc <- intersect(lipids, names(data))
  
  for (lipid in lipids_to_proc) {
    # 1. Prepare data and INT transform the raw lipid
    temp_df <- data %>% select(ID, sex, all_of(c(age_col, sex_col, other_covs, lipid)))
    temp_df$y <- inormal(temp_df[[lipid]])
    
    # 2. Apply 4SD outlier removal (consistent with brain/adiposity logic)
    temp_df$y <- ifelse(abs(scale(temp_df$y)[,1]) > 4, NA, temp_df$y)
    
    # 3. Define formula (excluding sex as it is the strata)
    cov_terms <- c(age_col, other_covs)
    mod_rhs <- paste(cov_terms, collapse = " + ")
    formula_str <- as.formula(paste("y ~", mod_rhs))
    
    # 4. Sex-stratified adjustment
    fitM <- lm(formula_str, data = temp_df %>% filter(sex == "M"), na.action = na.exclude)
    fitF <- lm(formula_str, data = temp_df %>% filter(sex == "F"), na.action = na.exclude)
    
    # 5. Extract residuals
    resM <- temp_df %>% filter(sex == "M") %>% mutate(!!lipid := resid(fitM)) %>% select(ID, !!lipid)
    resF <- temp_df %>% filter(sex == "F") %>% mutate(!!lipid := resid(fitF)) %>% select(ID, !!lipid)
    
    # Join back to main frame
    lipid_resids <- bind_rows(resM, resF)
    adjusted_data[[lipid]] <- lipid_resids[[lipid]][match(adjusted_data$ID, lipid_resids$ID)]
    
    capture.output(paste0("adjust_metabolites [", lipid, "]: ", formula_str), file=logfile, append=T)
  }
  return(adjusted_data)
}

# --- 8. Stratified or pooled analysis functions ---

 str_analysis = function(data,str.var,str.var.name,cov_cohort_specific,cov_additional,is.family.data){
  groups = unique(na.omit(data[[str.var.name]]))
  if(length(groups)>1){
    # 2. descriptive stats
    for(group in groups){
      outdir = file.path(outdir00,paste0(str.var,"_",str.var.name,"_",group))
      dir.create(outdir)
      df = data %>% filter(.data[[str.var.name]] == group)
      source("02_get_desc_statistics_revised.R", local = T)
      capture.output(cat("\n--- Adjustment definitions updated to version 2.1 (Sex-Stratified) ---\n"), 
               file = logfile, append = T)
      # check sample sizes for each categorical variable:
      small.sample_base <- check.sample.size(outdir,cov_cohort_specific)
      if(!is.null(small.sample_base$lessthan_5)){
        # discrete cohort-specific covs exist and can do sex-specific adjustment --> proceed to fitting
        run_base = !small.sample_base$lessthan_5
      }else{#is.null(small.sample_base$lessthan_5
        run_base = TRUE# discrete cohort-specific covs don't exist --> proceed to fitting
      }
      if(run_base) {
        if(is.family.data){
          msg = paste0('\n run family-based association tests under base model.\n')
          message(msg)
          cat(msg,file=logfile,append = T)
          source("3a_roiAssociation_statistics_baseModel_FAM.R",local = T)
        }else{
          source("3a_roiAssociation_statistics_baseModel.R",local = T)
        }
      }else{
        msg = paste0('\nBefore 3a: ',str.var,"_",group," (n=",nrow(df),
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
        
        if(is.family.data){
          msg = paste0('\n run family-based association tests under full model.\n')
          message(msg)
          cat(msg,file=logfile,append = T)
          
          source("3b_roiAssociation_statistics_fullModel_FAM.R",local = T)
        }else{
          source("3b_roiAssociation_statistics_fullModel.R",local = T)
        }
        source("3c_plot_roiAssociation_results_base_vs_fullModels.R",local=T)
      }else{
        msg = paste0('\nBefore 3b: ',str.var,"_",group," (n=",nrow(df),
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
    logfile=file.path(outdir,"strtified_analyses.log")
    df = data
    msg = paste0("\n There is only one group:", str.var,"_",group," (n=",nrow(df),"): No need to run ",str.var,"-stratified analyses.\n")
    message(msg)
    cat(msg,file=logfile,append = T)
  }
}

pooled_analysis = function(
    data,
    cov_cohort_specific,
    cov_additional,
    analysis.type,# c("primary","primary_all",'adjEthnicity')
    is.family.data)
  {
  group = "pooled"
  # 2. descriptive stats
  outdir = file.path(outdir00,paste0(group,"_",analysis.type))
  dir.create(outdir)
  df = data
  if(str_detect(analysis.type, 'primary')){
    do.run = TRUE
  }else{
    do.run = length(unique(df[['ethnicity']])) >1 #run if more than 1 ethnicity
  }
  if(do.run){
      source("02_get_desc_statistics_revised.R", local = T)
      capture.output(cat("\n--- Adjustment definitions updated to version 2.1 (Sex-Stratified) ---\n"), 
               file = logfile, append = T)
    # check sample sizes for each categorical variable:
    small.sample_base <- check.sample.size(outdir,cov_cohort_specific)
    if(!is.null(small.sample_base$lessthan_5)){
      # discrete cohort-specific covs exist and can do sex-specific adjustment --> proceed to fitting
      run_base = !small.sample_base$lessthan_5
    }else{#is.null(small.sample_base$lessthan_5
      run_base = TRUE# discrete cohort-specific covs don't exist --> proceed to fitting
    }
    if(run_base) {
      if(is.family.data){
        source("3a_roiAssociation_statistics_baseModel_FAM.R",local = T)
      }else{
        source("3a_roiAssociation_statistics_baseModel.R",local = T)
      }
    }else{
      msg = paste0('\nBefore 3a: ',group," (n=",nrow(df),
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
      if(is.family.data){
        source("3b_roiAssociation_statistics_fullModel_FAM.R",local = T)
      }else{
        source("3b_roiAssociation_statistics_fullModel.R",local = T)
      }
      source("3c_plot_roiAssociation_results_base_vs_fullModels.R",local=T)
    }else{
      msg = paste0('\nBefore 3b: ',group," (n=",nrow(df),
                   '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                   'we will not run association tests. \n')
      message(msg)
      print(small.sample_full$df_minFreq)
      cat(msg,file=logfile,append = T)
      capture.output(print(small.sample_full$df_minFreq),file = logfile,append=T)
    }#else for if(run_full)
  }else{#!do.run
    logfile=file.path(outdir,"pooled_analyses.log")
    msg = paste0("\n There is only one ethnic group, \'",unique(df$ethnicity),"\' (n=",nrow(df),"): We will not run ethnicity-adjusted analyses.\n")
    message(msg)
    capture.output(msg,file=logfile,append=T)
  }
}

# --- 9. Core Mediation Engines ---

run_single_mediation <- function(data, exp_var, med_var, out_var, covars, sims, X_label, Y_label, is.family.data) {
  
  data_sub <- data %>% select(ID, FID, all_of(c(exp_var, med_var, out_var, covars))) %>% na.omit()
  if (nrow(data_sub) < 30) return(NULL) 
  
  data_sub$X <- inormal(data_sub[[exp_var]])
  data_sub$M <- inormal(data_sub[[med_var]])
  data_sub$Y <- inormal(data_sub[[out_var]])
  
  cov_str <- if(length(covars) > 0) paste(" +", paste(covars, collapse = " + ")) else ""
  
  if(is.family.data) {
    med_mod <- lmer(as.formula(paste("M ~ X", cov_str, "+ (1|FID)")), data = data_sub)
    out_mod <- lmer(as.formula(paste("Y ~ X + M", cov_str, "+ (1|FID)")), data = data_sub)
  } else {
    med_mod <- lm(as.formula(paste("M ~ X", cov_str)), data = data_sub)
    out_mod <- lm(as.formula(paste("Y ~ X + M", cov_str)), data = data_sub)
  }
  
  med_out <- mediate(med_mod, out_mod, treat = "X", mediator = "M", sims = sims)
  
  return(tibble(adiposity = X_label, metabo_ID = med_var, brain = Y_label, 
                ACME_beta = med_out$d0, ACME_pval = med_out$d0.p, Prop_Med = med_out$n0))
}

# --- 4. Pipeline Execution Function ---

run_consortium_mediation <- function(analysis_data, roi_list, metabo_list, cov_genoPC, 
                                     cov_cohort_specific = NULL, cov_additional = NULL, 
                                     sims = 100, is.family.data = FALSE, n_cores = parallel::detectCores() - 1) {
  
  # Filter constant covariates
  all_covs <- c(cov_genoPC, cov_cohort_specific, cov_additional)
  valid_covs <- names(analysis_data)[sapply(analysis_data[all_covs], function(x) length(unique(na.omit(x))) > 1)]
  v_add <- intersect(cov_additional, valid_covs)

  # Parallel Setup
  plan(multisession, workers = n_cores)
  
  # Task setup (Assuming Step 1 Residualization is completed outside or before this)
  xnames <- intersect(c("BMI", "waist", "WHR"), names(analysis_data))
  task_grid <- expand.grid(x = xnames, y = roi_list, m = metabo_list, stringsAsFactors = F)

  message("Starting Parallel Mediation Sweep...")
  results_list <- pblapply(1:nrow(task_grid), function(i) {
    task <- task_grid[i,]
    run_single_mediation(analysis_data, task$x, task$m, task$y, v_add, sims, task$x, task$y, is.family.data)
  }, cl = "future")

  return(bind_rows(results_list))
}