#!/usr/bin/Rscript

### Script written by Anna Furtjes (for Lifetime Brain Atrophy GWAS)
### Adopted/modified by: Jean Shin, jjshin.research@gmail.com for any questions 
### Edited on: oct2025 
# install libraries

if(!is.element("pacman",installed.packages()[,1])){
  install.packages("pacman")
}

pacman::p_load(#check if a package is installed, if not, it attempts to install it
  sessioninfo, #
  tidyverse, #readr, dplyr, ggplot2, purrr, tidyr
  optparse, psych, data.table, tableone, reshape2,
  ggpubr,cowplot,patchwork, GGally,corrplot,
  mgcv,         # for adjusting continuous covariates using non-parametric smooth models
  lme4,lmerTest,# for fitting regression to family-data
  EnvStats,     # for outliers [grubs]
  mediation,    # for mediation analyses (stage 2)
  pbapply) # for showing progressing for mediation analyses

if(!is.element("hrbrthemes",installed.packages()[,1])){
  remotes::install_github("hrbrmstr/hrbrthemes")
};library(hrbrthemes)
