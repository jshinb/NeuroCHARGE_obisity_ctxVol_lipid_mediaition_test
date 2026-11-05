### Functions written by Anna Furtjes, afurtjes@ed.ac.uk for any questions
### Distributed as part of our lifetime brain atrophy GWAS consortium effort
### Last edited May 2025

select = dplyr::select
describe = psych::describe
mediate = mediation::mediate

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

successivelyReduceAge = function(data = Share, ageCut = ageCut){
  # object to store results
  storeNames = c("Age cut-off value", "Cor", "p", "ci_l", "ci_u", "Measure")
  store = as.data.frame(matrix(nrow = length(ageCut)*3, ncol = length(storeNames)))
  names(store) = storeNames
  
  # scale data
  data$ratio = data$ratio*(-1)
  data$resid = data$resid*(-1)
  
  # iterate over each age cut-off and calculate scores
  for(i in ageCut){
    # store which age cut off iteration this is
    loc = which(is.na(store$`Age cut-off value`))[1]
    store[loc:(loc+2),"Age cut-off value"] = i
    
    # cut sample
    Youngdata = data[which(data$age <= i),]
    
    # calculate correlations
    ## difference
    store[loc,"Cor"] =
      with(Youngdata, cor.test(age, diff))$estimate
    
    store[loc,"p"] =
      with(Youngdata, cor.test(age, diff))$p.value
    
    store[loc,"ci_l"] =
      with(Youngdata, cor.test(age, diff))$conf.int[1]
    
    store[loc,"ci_u"] =
      with(Youngdata, cor.test(age, diff))$conf.int[2]
    
    store[loc,"Measure"] = "Difference score"
    
    ## ratio
    store[loc+1,"Cor"] =
      with(Youngdata, cor.test(age, ratio))$estimate
    
    store[loc+1,"p"] =
      with(Youngdata, cor.test(age, ratio))$p.value
    
    store[loc+1,"ci_l"] =
      with(Youngdata, cor.test(age, ratio))$conf.int[1]
    
    store[loc+1,"ci_u"] =
      with(Youngdata, cor.test(age, ratio))$conf.int[2]
    
    store[loc+1,"Measure"] = "Ratio score"
    
    ## residual
    store[loc+2,"Cor"] =
      with(Youngdata, cor.test(age, resid))$estimate
    
    store[loc+2,"p"] =
      with(Youngdata, cor.test(age, resid))$p.value
    
    store[loc+2,"ci_l"] =
      with(Youngdata, cor.test(age, resid))$conf.int[1]
    
    store[loc+2,"ci_u"] =
      with(Youngdata, cor.test(age, resid))$conf.int[2]
    
    store[loc+2,"Measure"] = "Residual score"
    
  }
  return(store)
}

agePlot = function(data = agePlotData){
  # work out nreaks for y-axis
  breaks = round(seq(from = minAge-1, to = maxAge+1, by = 2) )
  # if there are too many and they can't be read, reduce
  if(length(breaks) > 20){
    breaks = breaks[breaks %% 5 ==0]
  }
  
  p=ggplot(data = data)+
    geom_point(aes(x = Cor, y = `Age cut-off value`, colour = Measure), alpha = 0.5)+
    geom_errorbar(aes(y = `Age cut-off value`, xmin = ci_l, xmax = ci_u, colour = Measure), alpha = 0.3)+
    geom_vline(xintercept = 0, color = "grey")+
    xlab("Lifetime brain atrophy\ncorrelation with age")+
    ylab("Maximum age in subset\n(cut-off in age in years)")+
    scale_y_continuous(limits = c(minAge-1, maxAge+1), breaks = breaks)+
    #scale_x_continuous(limits = c(-0.1, 0.2), breaks = seq(from = -0.5, to = 0.35, by = 0.1))+
    scale_color_manual(values = c("#D81B60","#FFC107","#004D40"))+
    ggtitle(paste0(opt$cohort_name, " (N = ", nrow(neuro),")"))+
    theme_bw()+
    theme(panel.border = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  return(p)
}


#* will think about the following for outlier detection/removal:
rm.sd4 = function(x){
  x.scaled = scale(x)[,1]
  print(range(x.scaled),na.rm=T)
  ind = !is.na(x.scaled)
  ind = ind & abs(x.scaled)>4
  x[ind] <- NA
  x
}

#* adjust brain outcomes ----
get_adj_ctxV = function(roii,d,covariate_names,remove.ID=T,logfile){
  #age, sex and ICV-corrected
  yname = paste(roii,"volume",sep="_")
  
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2) %>%
    mutate(icv.c = scale(ICV)[,1],icv.c2 = (scale(ICV)/2)^2)
  
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2','icv.c', 'icv.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2','icv.c','icv.c2',setdiff(covariate_names,c('age_mri',"ICV",'sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_volume.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_volume.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_volume.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_volume.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxV: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_volume.adjSex"),paste0(roii,"_volume.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_ctxTH = function(roii,d,covariate_names,remove.ID,logfile){
  #age and sex-corrected
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2)
  yname = paste(roii,"thickness",sep="_")
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2',setdiff(covariate_names,c('age_mri','sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_thickness.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_thickness.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_thickness.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_thickness.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxTH: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_thickness.adjSex"),paste0(roii,"_thickness.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_ctxSA = function(roii,d,covariate_names,remove.ID=T,logfile){
  #* [age and sex-corrected and **ICV**]
  yname = paste(roii,"area",sep="_")
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2) %>%
    mutate(icv.c = scale(ICV)[,1],icv.c2 = (scale(ICV)/2)^2)
  
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2','icv.c', 'icv.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2','icv.c','icv.c2',setdiff(covariate_names,c('age_mri',"ICV",'sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_area.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_area.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_area.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_area.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxSA: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_area.adjSex"),paste0(roii,"_area.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_ctxSA.woICV = function(roii,d,covariate_names,remove.ID=T,logfile){
  #* [age and sex-corrected - not used at the moment]
  #age and sex-corrected (but...what about ICV??? it is affected by ICV)
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2)
  yname = paste(roii,"area",sep="_")
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2',setdiff(covariate_names,c('age_mri','sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_area.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_area.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_area.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_area.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxSA: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_area.adjSex"),paste0(roii,"_area.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

#* adjust adiposity measures ----
get_adj_BMI = function(d,covariate_names,remove.ID,logfile){
  analdati = d %>%
    mutate(age = age_adiposity) %>%
    mutate(age.c = scale(age)[,1],)
  
  analdati[['x']] = inormal(analdati[["BMI"]])
  # remove outliers based on BMI
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c',setdiff(covariate_names,c('age_adiposity','sex'))),collapse = " + ")
  mod0 = paste("x ~ ",mod0,sep="")
  
  # age, sex and other covariates: 
  fit0M = lm(as.formula(mod0),data=analdati%>% filter(sex=="M"),na.action=na.exclude)
  fit0F = lm(as.formula(mod0),data=analdati%>% filter(sex=="F"),na.action=na.exclude)
  analdatiM = analdati %>% filter(sex=="M") %>% select(ID) %>% 
    mutate(BMI.adjSex = resid(fit0M)) %>%
    mutate(BMI.adj = resid(fit0M)+coef(fit0M)[1])
  analdatiF = analdati %>% filter(sex=="F") %>% select(ID) %>% 
    mutate(BMI.adjSex = resid(fit0F)) %>%
    mutate(BMI.adj = resid(fit0F)+coef(fit0F)[1])
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))

  capture.output(paste0("get_adj_BMI: ",mod0),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste("BMI","adjSex",sep="."),paste("BMI","adj",sep=".")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_WC = function(d,covariate_names,remove.ID=F,logfile){
  analdati = d %>%
    mutate(age = age_adiposity) %>%
    mutate(age.c = scale(age)[,1],)
  
  analdati[['x']] = inormal(analdati[["waist"]])
  # remove outliers based on waist
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c','height',setdiff(covariate_names,c('age_adiposity','sex','height'))),collapse = " + ")
  mod0 = paste("x ~ ",mod0,sep="")
  
  # age, sex and other covariates: 
  fit0M = lm(as.formula(mod0),data=analdati%>% filter(sex=="M"),na.action=na.exclude)
  fit0F = lm(as.formula(mod0),data=analdati%>% filter(sex=="F"),na.action=na.exclude)
  analdatiM = analdati %>% filter(sex=="M") %>% select(ID) %>% 
    mutate(WC.adjSex = resid(fit0M)) %>%
    mutate(WC.adj = resid(fit0M)+coef(fit0M)[1])
  analdatiF = analdati %>% filter(sex=="F") %>% select(ID) %>% 
    mutate(WC.adjSex = resid(fit0F)) %>%
    mutate(WC.adj = resid(fit0F)+coef(fit0F)[1])
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  capture.output(paste0("get_adj_WC: ",mod0),file=logfile,append=T)
  
  ret = subset(analdati,select=c("ID",paste("WC","adjSex",sep="."),paste("WC","adj",sep=".")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_WHR = function(d,covariate_names,remove.ID=F,logfile){
  analdati = d %>%
    mutate(age = age_adiposity) %>%
    mutate(age.c = scale(age)[,1],)
  
  analdati[['x']] = inormal(analdati[["WHR"]])
  # remove outliers based on waist
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c',setdiff(covariate_names,c('age_adiposity','sex'))),collapse = " + ")
  mod0 = paste("x ~ ",mod0,sep="")
  
  # age, sex and other covariates: 
  fit0M = lm(as.formula(mod0),data=analdati%>% filter(sex=="M"),na.action=na.exclude)
  fit0F = lm(as.formula(mod0),data=analdati%>% filter(sex=="F"),na.action=na.exclude)
  
  analdatiM = analdati %>% filter(sex=="M") %>% select(ID) %>% 
    mutate(WHR.adjSex = resid(fit0M)) %>%
    mutate(WHR.adj = resid(fit0M)+coef(fit0M)[1])
  analdatiF = analdati %>% filter(sex=="F") %>% select(ID) %>% 
    mutate(WHR.adjSex = resid(fit0F)) %>%
    mutate(WHR.adj = resid(fit0F)+coef(fit0F)[1])
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  capture.output(paste0("get_adj_WHR: ",mod0),file=logfile,append=T)
  
  ret = subset(analdati,select=c("ID",paste("WHR","adjSex",sep="."),paste("WHR","adj",sep=".")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

#* check sample sizes for discrete covariste (if 1, cannot run regression adjusting for that var)----
check.sample.size <- function(anal_outdir,disc_cov_names){
  load(file.path(anal_outdir,"descriptive_discrete_tables.RData"))
  sm.F = desc_Female_discrete
  sm.M = desc_Male_discrete
  disc_cov_names = intersect(names(sm.F),disc_cov_names)
  
  if(length(disc_cov_names)>0){
    minFreq.F <- minFreq.M <- minLevel.F <- minLevel.M <- c(); for(i in disc_cov_names){
      min.ind = which.min(sm.F[[i]] %>% pull(Count))
      minFreq.F <- c(minFreq.F, sm.F[[i]][["Count"]][min.ind])
      minLevel.F <- c(minLevel.F,as.character(sm.F[[i]][['Value']])[min.ind])
      
      min.ind = which.min(sm.M[[i]] %>% pull(Count))
      minFreq.M <- c(minFreq.M, sm.M[[i]][["Count"]][min.ind])
      minLevel.M <- c(minLevel.M,as.character(sm.M[[i]][['Value']])[min.ind])
    }
    df_minFreq = data.frame(vars=disc_cov_names,minLevel.F,minFreq.F,minLevel.M,minFreq.M)
    lessthan_5 = min(minFreq.F,minFreq.M)<5
  }else{#length(disc_cov_names)=0
    df_minFreq = NULL
    lessthan_5 = NULL
  }
  ret = list(df_minFreq=df_minFreq,lessthan_5 = lessthan_5)
  ret
}
#* stratified or pooled analysis functions ----
str_analysis = function(data,str.var,str.var.name,cov_cohort_specific,cov_additional,is.family.data){
  groups = unique(na.omit(data[[str.var.name]]))
  if(length(groups)>1){
    # 2. descriptive stats
    for(group in groups){
      outdir = file.path(outdir00,paste0(str.var,"_",str.var.name,"_",group))
      dir.create(outdir)
      df = data %>% filter(.data[[str.var.name]] == group)
      source("2_Descriptive_statistics.R",local = T)
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
    is.family.data,
    outdir00,
    group = 'pooled')
  {
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
    source("02_Descriptive_statistics_stage2.R",local = T)
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

generate_histograms <- function(outdir,metabo_OUT,metabo_missing_zero_info,metabo_columns){
  #create a sub-directory to store the pdf files
  dir.create(file.path(outdir,"plots"))
  names(metabo_missing_zero_info)[1] <- "Abbreviation2"
  
  d <- metabo_OUT
  summarydata = psych::describe(d)
  summarydata <- data.frame(Abbreviation2=rownames(summarydata),summarydata)
  summarydata <- metabo_missing_zero_info %>% left_join(summarydata,join_by(Abbreviation2)) 
  summarydata$Abbreviation2 <- factor(summarydata$Abbreviation2,levels=metabo_columns)
  summarydata = summarydata %>% arrange(Abbreviation2)
  
  myvars <- paste(summarydata$Abbreviation2[which(summarydata$n>=2)])
  i=1;for(l in seq(1,length(myvars), by=60)){
    l0 = l
    l1 = min((l+59),length(myvars))
    temp_vars <- myvars[l0:l1]
    temp <- subset(d,select=temp_vars)
    tmp.colnames = colnames(temp)
    temp = data.frame(temp)
    names(temp) = tmp.colnames
    temp_long <- melt(temp)
    temp_info <- summarydata %>% filter(Abbreviation2 %in% temp_vars)
    names(temp_info)[1] <- "variable"
    # generate plots
    myplot <- ggplot()+
      geom_histogram(data=temp_long, aes(x=value), fill="lightsteelblue2", bins=30)+
      geom_text(data=temp_info,
                aes(x=Inf, y=Inf, label=ifelse(mean<100, paste("mean", formatC(mean, format="f", digits=3)), paste("mean", round(mean, digits=0)))),
                hjust=1, vjust=2, size=2.3)+
      geom_text(data=temp_info,
                aes(x=Inf, y=Inf, label=ifelse(median<100, paste("median", formatC(median, format="f", digits=3)), paste("median", round(median, digits=0)))),
                hjust=1, vjust=3.5, size=2.3)+
      geom_text(data=temp_info,
                aes(x=Inf, y=Inf, label=paste("N", n)),
                hjust=1, vjust=5, size=2.3)+
      scale_y_continuous(expand=c(0,0))+
      facet_wrap(~variable, scales="free", ncol=10)+
      theme(strip.background=element_blank(),
            strip.text=element_text(size=8, face="bold"),
            axis.text.y=element_text(size=7),
            axis.text.x=element_text(size=7, angle=30, hjust=1, vjust=1),
            axis.title=element_text(size=8),
            plot.title=element_text(size=9, face="bold"))
    nvars <- length(temp_vars[which(!is.na(temp_vars))])
    h <- ifelse(nvars<11, 2.1, ifelse(nvars>=11 & nvars<21, 3.2, ifelse(nvars>=21 & nvars<31, 4.3, ifelse(nvars>=31 & nvars<41, 5.4, ifelse(nvars>=41 & nvars<51, 6.5, ifelse(nvars>=51, 7.6, NA))))))
    ggsave(file.path(outdir,"plots",paste0("histo",i,".pdf")), myplot, width=15, height=h)
    print(paste("variables",l, "-",(l+59),"done")) ; i=i+1  
  }
}#generate_histograms

adjust_metabolites <- function(data, lipids, age_col, sex_col, other_covs=NULL) {
  message("Adjusting metabolites for Age and Sex.\n")
  
  # just in case 
  adjusted_data <- data
  vars.all.missing = row.names(adjusted_data %>% describe %>% select(n) %>% filter(n==0))
  vars.zero.sd = row.names(adjusted_data %>% describe %>% select(sd) %>% filter(sd==0))
  adjusted_data = adjusted_data %>% select(-all_of(c(vars.all.missing,vars.zero.sd)))
  lipids = intersect(lipids,names(adjusted_data))
  i=0
  for (lipid in lipids) {
    i = i+1
    # 1. Fit linear model: Lipid ~ Age + Sex
    # We use na.exclude so the residuals line up with the original rows
    data[['metabo']] <- inormal(data[[lipid]]) #INT transformation
    mod= paste("metabo", "~", age_col, "*", sex_col)
    if(!is.null(other_covs)){
      other_covs = intersect(other_covs,names(adjusted_data))
      mod = paste0(
        c(paste("metabo", "~", age_col, "*", sex_col),other_covs), 
        collapse = ' + '
      )
    }
    formula_adj <- as.formula(mod)
    fit <- lm(formula_adj, data = data, na.action = na.exclude)
    if(i==1){
      print(mod)
    }
    # 2. Get residuals (The part of the lipid not explained by Age/Sex)
    # We add the mean back or scale it to keep the units interpretable
    adjusted_data[[lipid]] <- residuals(fit)
    
    # 3. Optional: Inverse Normal Rank Transformation (INT)
    # Common in consortiums to handle outliers and non-normality
    adjusted_data[[lipid]] <- qnorm((rank(adjusted_data[[lipid]], na.last="keep") - 0.5) / 
                                      sum(!is.na(adjusted_data[[lipid]])))
  }
  
  return(adjusted_data)
}

run_single_mediation_mediate <- function(data, exp_var, med_var, out_var, covars, sims, X_label = NULL, Y_label = NULL, do.boot=FALSE, is.family.data) {
  inormal <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  }
  
  # 1. Prepare Data Subset & Normalize
  # Using generic names 'x', 'm', 'y' simplifies formulas and ensures 'mediate' finds them
  relevant_vars <- c("ID","FID",exp_var, med_var, out_var, covars)
  data_sub <- na.omit(data %>% select(all_of(relevant_vars)))
  
  # Check if we have enough observations after na.omit
  if (nrow(data_sub) < 10) return(NULL) 
  
  data_sub[["x"]] <- inormal(data_sub[[exp_var]])
  data_sub[["m"]] <- inormal(data_sub[[med_var]])
  data_sub[["y"]] <- inormal(data_sub[[out_var]])
  
  # 2. Construct formulas
  covar_string <- if (length(covars) == 0) "" else paste(" +", paste(covars, collapse = " + "))
  exist.family = any(table(data$FID) > 0)
  run.lm = !is.family.data | !exist.family
  run.lmer = is.family.data & exist.family
  
  if(run.lm){
    f_m   <- paste("m ~ x", covar_string)#mediator model
    f_y   <- paste("y ~ x + m", covar_string)#Direct effect: outcome model
    f_pathb <- paste("y ~ m", covar_string) # Used for coefficient extraction
    f_total   <- paste("y ~ x", covar_string)#Total effect
    
    # 3. Fit Models
    # Using quote(data_sub) is essential for the bootstrap environment
    fit_m   <- do.call("lm", list(formula = as.formula(f_m), data = quote(data_sub)))
    fit_y   <- do.call("lm", list(formula = as.formula(f_y), data = quote(data_sub)))
    fit_pathb <- do.call("lm", list(formula = as.formula(f_pathb), data = quote(data_sub)))
    fit_total <- do.call("lm", list(formula = as.formula(f_total), data = quote(data_sub)))
    
    model.patha_pval = summary(fit_m)$coef['x', "Pr(>|t|)"]
    model.pathb_pval = summary(fit_pathb)$coef['m', "Pr(>|t|)"]
    ADE_pval2 = summary(fit_y)$coef['x', "Pr(>|t|)"]
    Total_pval2 = summary(fit_total)$coef['x',"Pr(>|t|)"]
    
  }else{#is.family.data=T
    f_m   <- paste("m ~ x", covar_string, "+ (1|FID)")
    f_y   <- paste("y ~ x + m", covar_string, "+ (1|FID)")#Direct effect
    f_pathb <- paste("y ~ m", covar_string, "+ (1|FID)") # Used for coefficient extraction
    f_total   <- paste("y ~ x", covar_string, "+ (1|FID)")#Total effect

    # Fit models
    fit_m   <- do.call(lmerTest::lmer, list(formula = as.formula(f_m), data = quote(data_sub)))
    fit_y   <- do.call(lmerTest::lmer, list(formula = as.formula(f_y), data = quote(data_sub)))
    fit_pathb <- do.call(lmerTest::lmer, list(formula = as.formula(f_pathb), data = quote(data_sub)))
    fit_total <- do.call(lmerTest::lmer, list(formula = as.formula(f_total), data = quote(data_sub)))
    
    model.patha_pval = summary(fit_m)$coef['x', "Pr(>|t|)"]
    model.pathb_pval = summary(fit_pathb)$coef['m', "Pr(>|t|)"]
    ADE_pval2 = summary(fit_y)$coef['x', "Pr(>|t|)"]
    Total_pval2 = summary(fit_total)$coef['x',"Pr(>|t|)"]
    
    class(fit_m) <- "lmerMod"
    class(fit_y) <- "lmerMod"
    class(fit_pathb) <- "lmerMod"
    class(fit_total) <- "lmerMod"
  }
  
  # 4. Run Mediation
  med_out <- tryCatch({
    mediation::mediate(
      model.m = fit_m, 
      model.y = fit_y,
      treat = 'x', 
      mediator = 'm', 
      boot = do.boot, #default=FALSE: a quasi-Bayesian approximation
      sims = sims
    )
  }, error = function(e) return(NULL))
  
  if (is.null(med_out)) return(NULL)
  
  # 5. Extract Statistics
  res <- data.frame(
    metabo_ID      = med_var,
    brain_original = Y_label %||% out_var, # Uses original name if label provided
    adiposity_orig = X_label %||% exp_var,
    ACME_est       = med_out$d0,
    ACME_pval      = med_out$d0.p,
    ADE_est        = med_out$z0,
    ADE_pval       = med_out$z0.p,
    Total_est      = med_out$tau.coef,
    Total_pval     = med_out$tau.p,
    Prop_Med       = med_out$n0,
    Lower_CI       = med_out$d0.ci[1],
    Upper_CI       = med_out$d0.ci[2],
    # path a: X->M
    model.patha_est  = summary(fit_m)$coef['x', "Estimate"],
    model.patha_se  = summary(fit_m)$coef['x', "Std. Error"],
    model.patha_pval =  model.patha_pval,
    # path b: M->Y
    model.pathb_est  = summary(fit_pathb)$coef['m', "Estimate"],
    model.pathb_se  = summary(fit_pathb)$coef['m', "Std. Error"],
    model.pathb_pval = model.pathb_pval,
    #direct exposure effect
    ADE_est2  = summary(fit_y)$coef['x', "Estimate"],
    ADE_se  = summary(fit_y)$coef['x', "Std. Error"],
    ADE_pval2 = ADE_pval2,
    # total exposure effect
    Total_est2  = summary(fit_total)$coef['x', "Estimate"],    
    Total_se = summary(fit_total)$coef['x',"Std. Error"],
    Total_pval2 = Total_pval2,
    #
    model.patha = f_m,
    model.pathb = f_pathb,
    model.ADE = f_y,
    model.Total = f_total,
    nsample = nobs(fit_y),
    stringsAsFactors = FALSE
  )
  
  # Trigger Garbage Collection to free RAM
  gc()
  
  return(res)
}

# Helper for the OR logic in labeling
`%||%` <- function(a, b) if (!is.null(a)) a else b

pooled_mediation_test <- function(
    data,# = df_all %>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"),
    is.family.data,#=opt$is.family.data
    cov_cohort_specific,# = cov_cohort_specific, 
    cov_additional,# = cov_additional, 
    analysis.type,# = 'primary',# c("primary","primary_all",'adjEthnicity')
    logfile.name,# logfile name,
    roi.34,# roi.34 = fread("freesurfer_34rois.txt",header=F) %>% pull(V1)
    group# = "pooled"
){
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
}#pooled_mediation_test ends

# str_mediation test
str_mediation_test = function(
    data,# = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]")
    str.var,# = "Ancestry"
    str.var.name,#= "ancestry"
    cov_cohort_specific,#=cov_cohort_specific
    cov_additional,# = cov_additional
    is.family.data,# = opt$is.family.data
    logfile.name
    ){# = "04_ancestry_str_mediation_analysis.log")
  
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
  
}#function ends

run_single_mediation_fast <- function(
    data, exp_var, med_var, out_var, covars,
    sims = 500,
    mode = c("screen", "mediate", "both"),
    X_label = NULL, Y_label = NULL,
    do.boot = FALSE,
    is.family.data = FALSE
) {
  
  inormal <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  }
  
  mode <- match.arg(mode)
  
  # --- 1. Subset (data.table + data.frame safe) ---
  relevant_vars <- c("ID","FID",exp_var, med_var, out_var, covars)
  
  if (inherits(data, "data.table")) {
    data_sub <- data[, ..relevant_vars]
    data_sub <- na.omit(data_sub)
  } else {
    data_sub <- data[relevant_vars]
    data_sub <- data_sub[complete.cases(data_sub), ]
  }
  
  if (nrow(data_sub) < 10) return(NULL)
  
  # --- 2. Normalize ---
  data_sub$x <- inormal(data_sub[[exp_var]])
  data_sub$m <- inormal(data_sub[[med_var]])
  data_sub$y <- inormal(data_sub[[out_var]])
  
  # --- 3. Model setup ---
  covar_string <- if (length(covars)) paste("+", paste(covars, collapse = " + ")) else ""
  exist.family <- any(duplicated(data$FID))
  
  run.lm   <- !is.family.data || !exist.family
  run.lmer <-  is.family.data && exist.family
  
  if (run.lm) {
    f_m <- as.formula(paste("m ~ x", covar_string))
    f_y <- as.formula(paste("y ~ x + m", covar_string))
    
    fit_m <- lm(f_m, data = data_sub)
    fit_y <- lm(f_y, data = data_sub)
    
  } else {
    f_m <- as.formula(paste("m ~ x", covar_string, "+ (1|FID)"))
    f_y <- as.formula(paste("y ~ x + m", covar_string, "+ (1|FID)"))
    
    fit_m <- lmerTest::lmer(f_m, data = data_sub)
    fit_y <- lmerTest::lmer(f_y, data = data_sub)
    
    class(fit_m) <- "lmerMod"
    class(fit_y) <- "lmerMod"
  }
  
  # --- 4. Extract coefficients once ---
  sm <- coef(summary(fit_m))
  sy <- coef(summary(fit_y))
  
  a  <- sm["x", 1]
  sa <- sm["x", 2]
  
  b  <- sy["m", 1]
  sb <- sy["m", 2]
  
  c_prime <- sy["x", 1]
  
  # --- 5. Stage 1 (fast mediation) ---
  ab <- a * b
  
  # Aroian SE
  se_ab <- sqrt(b^2 * sa^2 + a^2 * sb^2 + sa^2 * sb^2)
  
  z_ab <- ab / se_ab
  p_ab <- 2 * pnorm(-abs(z_ab))
  
  # total effect approximation
  c_total <- c_prime + ab
  prop_med <- ab / c_total
  
  res_screen <- data.frame(
    metabo_ID = med_var,
    brain_original = Y_label %||% out_var,
    adiposity_orig = X_label %||% exp_var,
    
    a = a, sa = sa,
    b = b, sb = sb,
    
    #BASED ON SOBEL 
    ACME_est = ab,
    ACME_se  = se_ab,
    ACME_z   = z_ab,
    ACME_pval = p_ab,
    
    c_prime = c_prime,
    c_total = c_total,
    prop_med = prop_med,
    
    nsample = nobs(fit_y),
    stringsAsFactors = FALSE
  )
  
  # --- 6. Stage 2 (mediate) ---
  res_mediate <- NULL
  
  if (mode %in% c("mediate", "both")) {
    
    med_out <- tryCatch(
      mediation::mediate(
        model.m = fit_m,
        model.y = fit_y,
        treat = "x",
        mediator = "m",
        sims = sims,
        boot = do.boot
      ),
      error = function(e) NULL
    )
    
    if (!is.null(med_out)) {
      res_mediate <- data.frame(
        ACME_est_mediate  = med_out$d0,
        ACME_pval_mediate = med_out$d0.p,
        ADE_est           = med_out$z0,
        ADE_pval          = med_out$z0.p,
        Total_est         = med_out$tau.coef,
        Total_pval        = med_out$tau.p,
        Prop_Med_mediate  = med_out$n0,
        Lower_CI          = med_out$d0.ci[1],
        Upper_CI          = med_out$d0.ci[2]
      )
    }
  }
  
  # --- 7. Return ---
  if (mode == "screen") {
    return(res_screen)
  } else if (mode == "mediate") {
    return(res_mediate)
  } else {
    return(cbind(res_screen, res_mediate))
  }
}