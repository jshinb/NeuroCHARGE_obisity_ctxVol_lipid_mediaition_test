#*****************************************************************************#
#
# 03_mediation_analysis_fullModel.R: 
# Perform mediation tests: Adiposity -> Lipids -> Brain
# 
#*****************************************************************************#

## starting ----
op <- options(warn=1)
starting.message=paste0("\n",'03_mediation_analysis_fullModel.R: ',"Started mediation tests. \n\n")
message(starting.message)
cat(starting.message,file=logfile00,append=T)

## remove the identified columns from the analyses ----
## rm.var will be created within the function
##* age_diff: age_mri-age_adiposity; age_diff2 = age_mri - age_metabo
cov_ctxV = setdiff(c('age_mri','age_diff','age_diff2','sex','ICV',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
cov_ctxSA = setdiff(c('age_mri','age_diff','age_diff2','sex','ICV',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
cov_ctxTH = setdiff(c('age_mri','age_diff','age_diff2','sex',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)

# adjust ctx outcomes ----
adj_ctxV <- adj_ctxSA <- adj_ctxTH <-analysis_data %>% select(ID)
for (roii in c(roi.34,'global')){ 
  adj_ctxV = adj_ctxV %>% 
    left_join(get_adj_ctxV(
      roii=roii,
      d=analysis_data%>% rename(global_volume = total_volume),
      covariate_names = cov_ctxV,remove.ID=F,logfile=logfile00),join_by(ID))
  adj_ctxSA = adj_ctxSA %>% 
    left_join(get_adj_ctxSA(
      roii=roii,d=analysis_data%>% rename(global_area = total_area),
      covariate_names = cov_ctxSA,remove.ID=F,logfile=logfile00),join_by(ID))
  adj_ctxTH = adj_ctxTH %>% 
    left_join(get_adj_ctxTH(
      roii=roii,d=analysis_data%>% rename(global_thickness = mean_thickness),
      covariate_names = cov_ctxTH,remove.ID=F,logfile=logfile00),join_by(ID))
}

# adjust adiposity measures ----
if(has_element(names(analysis_data),"BMI")){
  cov_BMI = intersect(
    names(analysis_data),
    c('age_adiposity','age_diff','age_diff2','sex',cov_genoPC,cov_cohort_specific,cov_additional)
  )#setdiff(c('age_adiposity','age_diff','sex',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
  adj_BMI = get_adj_BMI(d=analysis_data,covariate_names = cov_BMI,remove.ID = F,logfile=logfile00)
}else{
  cov_BMI <- adj_BMI <- c()
}
if(has_element(names(analysis_data),"waist")){
  cov_WC = intersect(
    names(analysis_data),
    c('age_adiposity','age_diff','age_diff2','sex','height',cov_genoPC,cov_cohort_specific,cov_additional)
  )#setdiff(c('age_adiposity','age_diff','sex','height',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
  adj_WC = get_adj_WC(d=analysis_data,covariate_names = cov_WC,remove.ID = F,logfile=logfile00)
}else{
  cov_WC <- adj_WC <- c()
}
if(has_element(names(analysis_data),"WHR")){
  cov_WHR = intersect(
    names(analysis_data),
    c('age_adiposity','age_diff','age_diff2','sex',cov_genoPC,cov_cohort_specific,cov_additional)
  )#setdiff(c('age_adiposity','age_diff','sex',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
  adj_WHR = get_adj_WHR(d=analysis_data,covariate_names = cov_WHR,remove.ID = F,logfile=logfile00)
}else{
  cov_WHR <- adj_WHR <- c()
}

# adjust metabo ----
lipids = intersect(names(analysis_data),metaboID_ordered) #filtered
adj_lipid <- adjust_metabolites(
  data = subset(
    analysis_data,
    select=intersect(
      names(analysis_data),c("ID",'age_metabo','age_diff','age_diff2','sex',lipids)
    )
  ), 
  lipids = lipids, 
  age_col = "age_metabo",# 
  sex_col = "sex",
  other_covs = c('age_diff',"age_diff2")
)

# combine adjusted data ----
adjusted_data = analysis_data %>% select(all_of(c("ID","FID",cov_additional)))
if(!is.null(adj_ctxV)){
  adjusted_data = adjusted_data %>% left_join(adj_ctxV,by="ID")
}
if(!is.null(adj_ctxSA)){
  adjusted_data = adjusted_data %>% left_join(adj_ctxSA,by="ID")
}
if(!is.null(adj_ctxTH)){
  adjusted_data = adjusted_data %>% left_join(adj_ctxTH,by="ID")
}
if(!is.null(adj_BMI)){
  adjusted_data = adjusted_data %>% left_join(adj_BMI,by="ID")
}
if(!is.null(adj_WC)){
  adjusted_data = adjusted_data %>% left_join(adj_WC,by="ID")
  names(adjusted_data) = str_replace(names(adjusted_data),"WC","waist")
}
if(!is.null(adj_WHR)){
  adjusted_data = adjusted_data %>% left_join(adj_WHR,by="ID")
}
if(!is.null(adj_lipid)){
  adjusted_data = adjusted_data %>% left_join(adj_lipid,by="ID")
}

# execute in parallel ----
xnames <- intersect(c("BMI","waist","WHR"),str_split(names(adjusted_data),"[.]",simplify = T)[,1])
ynames <- intersect(ctx_vars,str_split(names(adjusted_data),"[.]",simplify = T)[,1])
msg = paste("Starting mediation for", length(lipids), "lipids, ", length(xnames), 'adiposity variable(s), and',
            length(ynames), "brain phenotypes.")
message(msg);cat("\n\n",msg,"\n\n",sep='',file=logfile00, append=T)

## 1. create the grid of all 3 loops ----
# Total tasks = length(xnames) * length(ynames) * length(lipids)
task_grid <- expand.grid(
  x_name = xnames,
  y_name = ynames,
  lipid  = lipids,
  stringsAsFactors = FALSE
)

# Force warnings to stay "buffered" so the handler can catch them properly
options(warn = 0) 

log_file <- file.path(outdir, "mediation_warnings_fullModel.log")
write(paste("Mediation Analysis Log - Started:", Sys.time()), file = log_file)

t0=1;t1=nrow(task_grid)

# ===== CLUSTER SETUP =====
require(parallel)

# Create cluster (adjust number of cores as needed)
### CHANGED to be more robust
# n_cores <- as.numeric(
#   Sys.getenv("SLURM_CPUS_PER_TASK", detectCores() - 1))
n_cores <- 4L
cl <- makeCluster(n_cores, type = "PSOCK")

# Export all required objects to cluster nodes
clusterExport(cl, c(
  "adjusted_data", 
  "task_grid",
  "opt",
  "run_single_mediation_fast",
  "log_file",
  "t0",
  "t1",
  "inormal"
), envir = environment())

# Load required packages on each node
clusterEvalQ(cl, {
  require(dplyr)
  require(mediation)
  require(lmerTest)
  require(lme4)
})

# ===== PARALLEL LOOP =====
final_results_list <- parLapply(cl, c(t0:t1), function(i) {
  task <- task_grid[i, ]
  
  # Use withCallingHandlers for BOTH warnings and messages
  withCallingHandlers(
    {
      run_single_mediation_fast(
        data=adjusted_data, 
        exp_var = paste0(task$x_name, ".adjSex"), 
        med_var = task$lipid, 
        out_var = paste0(task$y_name, ".adjSex"), 
        covars = NULL,
        mode = "screen",
        X_label = task$x_name, 
        Y_label = task$y_name,
        sims = 500,
        do.boot = FALSE,
        is.family.data = opt$is.family.data
      )
    },
    # Handle Warnings
    warning = function(w) {
      log_entry <- sprintf("[%s] Row %d | Lipid: %s | Warning: %s", 
                           Sys.time(), i, task$lipid, conditionMessage(w))
      write(log_entry, file = log_file, append = TRUE)
      if (!is.null(findRestart("muffleWarning"))) invokeRestart("muffleWarning")
    },
    # Handle Messages (Common in lmerTest/lme4)
    message = function(m) {
      log_entry <- sprintf("[%s] Row %d | Lipid: %s | Message: %s", 
                           Sys.time(), i, task$lipid, conditionMessage(m))
      write(log_entry, file = log_file, append = TRUE)
      if (!is.null(findRestart("muffleMessage"))) invokeRestart("muffleMessage")
    }
  )
})

# Stop cluster
stopCluster(cl)

# 4. Clean and Save
final_df <- bind_rows(final_results_list)
resfile = file.path(outdir,"mediation_results_summary_fullModel.txt")
write_tsv(final_df, resfile)

capture.output(summary(warnings()),file=logfile00,append=T)
ending.message=paste0("\n",'03_mediation_analysis_fullModel.R: ',"Analysis complete. Results saved to output directory. \n\n")
message(ending.message); cat(ending.message,"\n",file=logfile00,append=T)
