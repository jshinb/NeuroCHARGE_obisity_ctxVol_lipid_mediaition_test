#*****************************************************************************#
#
# Step 2: Descriptive statistics for data set (including metabodata): table and plots
#*****************************************************************************#
op <- options(warn=1)

starting.message=paste0("\n",'02_get_desc_statistics_stage2.R: ',"Started obtaining descriptive statistics. \n\n")
message(starting.message)
cat(starting.message,file=logfile00,append=T)

# Descriptive statistics ------------------------------------------------------
if(nrow(analysis_data) < 30){
  msg=paste0("There were only ", nrow(analysis_data), " participants left in ", group,". Skipping analysis of this group\n\n.")
  message(msg); 
  cat(msg,file=logfile00,append=T)
  next
}

#* remove variables that are 100% missing in the merged dataset ----
data_missing <- colSums(is.na(analysis_data))
if(any(data_missing == nrow(analysis_data))){
  print(names(analysis_data)[data_missing == nrow(analysis_data)])
  missingCol <- names(analysis_data)[data_missing == nrow(analysis_data)]
  print(missingCol)
  analysis_data[, missingCol] <- NULL
}else{
  missingCol = NULL
}

metabo_columns = names(analysis_data)[names(analysis_data) %in% metaboID_ordered]
metabo_OUT<-subset(analysis_data, select = metabo_columns) #Metabolites only
metabo_OUT_males<-subset(analysis_data, sex=="M", select = metabo_columns) #Male metabolites only
metabo_OUT_females<-subset(analysis_data, sex=="F", select = metabo_columns) #Female metabolites only

metabo_missing_count <- colSums(is.na(metabo_OUT))
metabo_zero_count <- colSums(metabo_OUT == 0, na.rm = TRUE)
metabo_missing_count <- data.frame(metabo_ID=names(metabo_missing_count),n_missing=metabo_missing_count)
metabo_zero_count <- data.frame(metabo_ID=names(metabo_zero_count),n_zero=metabo_zero_count)
metabo_missing_zero_info <- metabo_missing_count %>% 
  left_join(metabo_zero_count,join_by(metabo_ID)) %>% 
  mutate(nsample = nrow(metabo_OUT)-n_missing) 

filename3 = paste("metabo_missing_zero_freqs.tsv", sep="_")
write_tsv(metabo_missing_zero_info,file.path(outdir,filename3))

#* record the numbers of missing values for each column
to_write = data.frame(data_missing) %>% filter(!rownames(.) %in% metaboID_ordered);
to_write = data.frame(vars = rownames(to_write),to_write)
file_to_write=file.path(outdir,paste("table_data_missing.txt",sep=""))
write_tsv(to_write,file_to_write);rm(to_write,file_to_write)

#*****************************************************************************#
# Standardise and remove outliers for cortical volume outside SD of 4 ----
#*****************************************************************************#
#* record number of participants with regional volume with SD >4
ctx_vars.ind = str_detect(names(analysis_data),"_volume");sum(ctx_vars.ind)
ctx_vars.ind = ctx_vars.ind | str_detect(names(analysis_data),"_thickness");sum(ctx_vars.ind)
ctx_vars.ind = ctx_vars.ind | str_detect(names(analysis_data),"_area");sum(ctx_vars.ind)
ctx_vars = names(analysis_data)[ctx_vars.ind]
df_tmp = analysis_data;count_4SD = c();for(xi in ctx_vars) {
  x =  df_tmp[[xi]]
  ind.4SD = abs(scale(x)[,1])>4
  count_4SD = c(count_4SD,sum(ind.4SD, na.rm=T))
  x[!is.na(ind.4SD) & ind.4SD] <- NA
  df_tmp[[xi]] <- x
}

df.count_4SD = data.frame(vars=ctx_vars,count_4sd = count_4SD)
write_tsv(df.count_4SD,file=file.path(outdir,"df.counts_4sd.txt"))
analysis_data= df_tmp;rm(df_tmp)

#*****************************************************************************#
# Summarise neuroimaging data (moved to 2_Descriptive_statistics) ----
#*****************************************************************************#
#* plot histogram
n = nrow(analysis_data)
if(!any("BMI"==missingCol)){
  phist_BMI <- plot_hist(dat = analysis_data, var = "BMI", col = "#D81B60")
}else{
  phist_BMI <- NULL
}
if(!any("waist"==missingCol)){
  phist_waist <- plot_hist(dat = analysis_data, var = "waist", col = "#FFC107")
}else{
  phist_waist <- NULL
}
if(!any("WHR"==missingCol)){
  phist_WHR <- plot_hist(dat = analysis_data, var = "WHR", col = "#004D40")
}else{
  phist_WHR <- NULL
}
#* arrange plots (with safest annotating option)
p <- (phist_BMI + phist_waist + phist_WHR) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = paste0(opt$cohort_name, " (", group, ")\nN = ", n)
  )
p
#* save plot
p_hist_file=file.path(outdir,"histograms.jpg")
ggsave(p_hist_file, 
       bg = "white",
       plot = p, 
       width = 17, 
       height = 7, 
       units = "cm", 
       dpi = 150)

#* plot heat map
#* suppress warnings because a heatmap always has missing values 
suppressWarnings({
  vnames = unlist(str_split('BMI,waist,WHR,height,total_volume,total_area,mean_thickness,ICV,age_adiposity,age_mri',','))
  vnames = setdiff(vnames,missingCol)
  heat = plot_heatmap(dat = analysis_data%>% select(all_of(vnames)),
                      axisNames = vnames)
  
  # annotate
  heat <- annotate_figure(heat, 
                          top = text_grob(paste0(opt$cohort_name, " (", group, ")\nN = ",n),
                                          color = "black", 
                                          face = "bold", 
                                          size = 14))
  
  # save plot
  ggsave(file.path(outdir,"corrplot.jpg"), 
         bg = "white",
         plot = heat, 
         width = 15, 
         height = 15, 
         units = "cm", 
         dpi = 150)
})

#* 
desc_All = analysis_data%>% select(-ID,-FID) %>% psych::describe(IQR = T,quant=c(0.25,0.75))
desc_FemalesMales = analysis_data%>% select(-ID,-FID) %>% 
  mutate(sex = as.character(sex)) %>% 
  psych::describeBy(group="sex",IQR=T,quant=c(0.25,0.75))
save(desc_All,desc_FemalesMales,file=file.path(outdir,"descriptive_tables.RData"))

# categorical variables ----
descriptive_stats_discrete <- function(x,data){
  ret = table(data[[x]])
  ret = data.frame(ret,prop.table(ret))
  ret = subset(ret,select=-Var1.1)
  names(ret) = c('Value',"Count","Proportion")
  ret
}#descriptive_stats_discrete

cat_vars = attr(unlist(lapply(analysis_data, class)),'names')[which(unlist(lapply(analysis_data, class))%in%c("factor","character"))]
cat_vars = setdiff(cat_vars,"ID")

desc_All_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(analysis_data,select = cat_vars))
names(desc_All_discrete) = cat_vars

desc_Female_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(analysis_data,select = cat_vars, sex=="F"))
names(desc_Female_discrete) = cat_vars

desc_Male_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(analysis_data,select = cat_vars, sex=="M"))
names(desc_Male_discrete) = cat_vars

save(desc_All_discrete,desc_Female_discrete,desc_Male_discrete,
     file=file.path(outdir,"descriptive_discrete_tables.RData"))

# correlation matrices and their plots for 'all' variables ----
corm = cor(subset(analysis_data,select=setdiff(names(analysis_data),c("ID",cat_vars))),use="p")
cormF = cor(subset(analysis_data,sex=="F",select=setdiff(names(analysis_data),c("ID",cat_vars))),use="p")
cormM = cor(subset(analysis_data,sex=="M",select=setdiff(names(analysis_data),c("ID",cat_vars))),use="p")
cor_list = list(corm=corm, cormF=cormF, cormM=cormM)
save(cor_list,file=file.path(outdir,'correlation_matrices_continuous_variables.Rdata'))

# scatter plots with correlation ----
dA=subset(analysis_data,select=c( 
  setdiff(intersect(names(non.neuro),names(analysis_data)),c('ID',"FID",missingCol)), "ICV"))

## sex-combined
png(file.path(outdir,"correlation_pairwise_scatter.png"),
    width=22,height=22,units="in",res=300)
pairs.panels(data.frame(dA,stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#BEBEBE33",
             density = TRUE,  # show density plots
             ellipses = TRUE) # show correlation ellipses
dev.off()

## females
png(file.path(outdir,"correlation_pairwise_scatter_Female.png"),
    width=22,height=22,units="in",res=300)
pairs.panels(data.frame(subset(dA,sex=="F"),stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#8B000033",
             density = TRUE,  # show density plots
             ellipses = TRUE) # show correlation ellipses
dev.off()

## males
png(file.path(outdir,"correlation_pairwise_scatter_Male.png"),
    width=22,height=22,units="in",res=300)
pairs.panels(data.frame(subset(dA,sex=="M"),stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#00008B33",
             density = TRUE,  # show density plots
             ellipses = TRUE) # show correlation ellipses
dev.off()

# characteristic tables ----
# 
tab_vars = setdiff(c(names(analysis_data)),c("ID","FID"))
tableOne_All <- CreateTableOne(vars = tab_vars,
                               #strata = "trt", 
                               data = analysis_data, 
                               factorVars = cat_vars)

tableOne_sex_str <- CreateTableOne(vars = tab_vars,
                                   strata = "sex", 
                                   data = analysis_data, 
                                   factorVars = cat_vars)

if(length(table(analysis_data$is.underweight, useNA = "ifany"))>1){
  tableOne_bmi_str <- CreateTableOne(vars = tab_vars,
                                     strata = "is.underweight", 
                                     data = analysis_data, 
                                     factorVars = cat_vars)
  save(tableOne_bmi_str,
       file=file.path(outdir,
                      "characteristic_table_bmi_stratified.RData"))
}else{
  message("\n02_get_desc_statistics_stage2.R: Only one bmi-group exisits.\n")
}

if(length(table(analysis_data$age.group, useNA = "ifany"))>1){
  tableOne_age_str <- CreateTableOne(vars = tab_vars,
                                     strata = "age_group",
                                     data = analysis_data, 
                                     factorVars = cat_vars)
  save(tableOne_age_str,
       file=file.path(outdir,"characteristic_table_age_stratified.RData"))
}else{
  message("\n02_get_desc_statistics_stage2.R: Only one age-group exisits.\n")
}
save(tableOne_All,
     file=file.path(outdir,
                    "characteristic_table.RData"))
save(tableOne_sex_str,
     file=file.path(outdir,
                    "characteristic_table_sex_stratified.RData"))
# generate histograms
generate_histograms(outdir = outdir ,
                    metabo_OUT,
                    metabo_missing_zero_info, 
                    metabo_columns = names(metabo_OUT))

# ending ----
capture.output(summary(warnings()), file=logfile00,append=T)
ending.message = paste0("\n",'02_get_desc_statistics_stage2.R: ',"Finished obtaining descriptive statistics. \n\n")
message(ending.message,sep='')
cat(ending.message,sep='',file=logfile00,append=T)

# clearing ----
options(op)