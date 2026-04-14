#' ---
#' title: Determinants of PFAS exposure 
#' subtitle: Data analysis sensitivity on women
#' author: Julien Riou
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     code_folding: hide
#'     theme: cosmo
#'     highlight: pygments
#'     fig_width: 10
#'     fig_height: 8
#'     fig_caption: true
#' params:
#'   analysis_date:
#' bibliography: misc/bib.bib  
#' ---

#+ results="hide", warnings="false", echo="false"
analysis_date = "2026-02-23"
run_models = TRUE
controls = readRDS(paste0("savepoints/savepoint_",analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds"))

controls$covariates_list <- controls$covariates_list[-3]
controls$covariates_list <- c(controls$covariates_list,"hm_kids","hm_breastfeeding","hm_menopause")
# controls$labs$hm_kids <- "Number of biological children"

# keep only relevant variables and remove men
pei = shesp_3  %>% 
  dplyr::filter(gender=="Females") %>% 
  dplyr::select(did,
                pfas_cluster,
                all_of(paste0(controls$pfas_list_retained,"_nonzero")),
                all_of(paste0(controls$pfas_list_retained,"_quant")),
                all_of(controls$covariates_list))

#' # Description
#' 

fn_105_summary_pfas_cluster(dat_in=pei,
                            cap="PFAS measurements within clusters of overall PFAS exposition.",
                            exclusions=controls$pfas_allzero)
pei  %>%
  dplyr::select(-did,-starts_with("pfas")) %>%
  frq_table(cap="Covariates.",missing="ifany")



#' # Results by PFAS cluster
#' 
#' ## Univariable regression
#' 

if(TRUE) {
  fn_010_clean_dlls() # removes temporary files
  uni_ = fn_207_multinomial(dat_in=pei, 
                            covariates=controls$covariates_list, 
                            multiple_runs=TRUE,
                            multivariable=FALSE,
                            impute=c("age","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label")
  )
  write_rds(uni_,file=file.path(controls$savepoint,"4_univariable_pfas_cluster_women.rds"))
}

uni_ = readRDS(file.path(controls$savepoint,"4_univariable_pfas_cluster_women.rds"))
#+ fig.width=11, fig.height=18
fn_309_fig_regression_cluster(uni_)

#' ## Multivariable regression adjusted on age
#' 

if(TRUE) {
  fn_010_clean_dlls() # removes temporary files
  mul_ = fn_207_multinomial(dat_in=pei, 
                            covariates=controls$covariates_list, 
                            multiple_runs=TRUE,
                            multivariable=TRUE,
                            regularization=FALSE,
                            impute=c("age","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                            adj=c("age")
  )
  write_rds(mul_,file=file.path(controls$savepoint,"4_multivariable_pfas_cluster_women.rds"))
}

mul_ = readRDS(file.path(controls$savepoint,"4_multivariable_pfas_cluster_women.rds"))
#+ fig.width=11, fig.height=18
fn_309_fig_regression_cluster(mul_)


#' ## Multivariable regression adjusted on age, sex, nationality, income, education, center
#' 

if(TRUE) {
  fn_010_clean_dlls() # removes temporary files
  mul2_ = fn_207_multinomial(dat_in=pei, 
                             covariates=controls$covariates_list, 
                             multiple_runs=TRUE,
                             multivariable=TRUE,
                             regularization=FALSE,
                             impute=c("age","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                             adj=c("age","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label")
  )
  write_rds(mul2_,file=file.path(controls$savepoint,"4_multivariable2_pfas_cluster_women.rds"))
}

mul2_ = readRDS(file.path(controls$savepoint,"4_multivariable2_pfas_cluster_women.rds"))
#+ fig.width=11, fig.height=18
fn_309_fig_regression_cluster(mul2_)


#' ## Multivariable regression based on top 10 from BART variable selection
#' 

if(TRUE) {
  controls$covariates_list = controls$covariates_list[!controls$covariates_list %in% c("age_group")]
  sel_ = fn_205_bart_selection(dat_in=pei,
                               dependent="pfas_cluster",
                               covariates=controls$covariates_list,
                               impute=c("age","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
  sel_ = names(sel_$var_importance[1:10])
  bmul_ = fn_207_multinomial(dat_in=pei, 
                             covariates=sel_, 
                             multiple_runs=FALSE,
                             multivariable=TRUE,
                             regularization=FALSE,
                             impute=c("age","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                             adj=sel_
                             
  )
  write_rds(bmul_,file=file.path(controls$savepoint,"4_bart_multivariable_pfas_cluster_women.rds"))
}

bmul_ = readRDS(file.path(controls$savepoint,"4_bart_multivariable_pfas_cluster_women.rds")) %>%
  dplyr::mutate(adjusted="Adjusted, BART selection")
#+ fig.width=11, fig.height=18
fn_309_fig_regression_cluster(bmul_)

#' ## Multivariable penalized regression
#' 


if(TRUE) {
  fn_010_clean_dlls() # removes temporary files
  pmul_ = fn_207_multinomial(dat_in=pei, 
                             covariates=controls$covariates_list, 
                             multiple_runs=TRUE,
                             multivariable=TRUE,
                             regularization=TRUE,
                             impute=c("age","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label")
  )
  write_rds(pmul_,file=file.path(controls$savepoint,"4_pen_multivariable_pfas_cluster_women.rds"))
}

pmul_ = readRDS(file.path(controls$savepoint,"4_pen_multivariable_pfas_cluster_women.rds"))  %>%
  dplyr::mutate(adjusted="Adjusted, penalization")
#+ fig.width=11, fig.height=18
fn_309_fig_regression_cluster(pmul_)

#' # Results by covariate
#' 

all_ = dplyr::bind_rows(uni_,mul_,mul2_,bmul_)

#' ## Socio-demographic
#' 
#' ### Age group

asso = fn_302_table_regression_by_covariate(cov_="age",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Sex

asso = fn_302_table_regression_by_covariate(cov_="gender",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Foreign nationality

asso = fn_302_table_regression_by_covariate(cov_="a1_foreign",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Education

asso = fn_302_table_regression_by_covariate(cov_="a1_education_3classes",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Household monthly income

asso = fn_302_table_regression_by_covariate(cov_="a1_household_monthly_income_3classes",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Study site

asso = fn_302_table_regression_by_covariate(cov_="center_label",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ## Residence
#' 
#' ### Residence in urban or industrial area compared to rural area

asso = fn_302_table_regression_by_covariate(cov_="ex1_residence_urban",coefs_=all_,dat_in=pei,cov_dataname="ex1_residence_urban")
asso[[1]]
asso[[2]]

#' ### Synthetic or natural carpet at home 

asso = fn_302_table_regression_by_covariate(cov_="ex1_carpet",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ## Breastfeeding
#' 

asso = fn_302_table_regression_by_covariate(cov_="hm_breastfeeding",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]


#' ## Food
#' 

#' ### Food group
#' 
asso = fn_302_table_regression_by_covariate(cov_="cluster",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Tuna (frequent consumption of)

asso = fn_302_table_regression_by_covariate(cov_="ex1_tuna",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ###  Salmon (frequent consumption of)

asso = fn_302_table_regression_by_covariate(cov_="ex1_salmon",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Pangasius (frequent consumption of)

asso = fn_302_table_regression_by_covariate(cov_="ex1_pangasius",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Trout (frequent consumption of)

asso = fn_302_table_regression_by_covariate(cov_="ex1_trout",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Shrimp (frequent consumption of)

asso = fn_302_table_regression_by_covariate(cov_="ex1_shrimp",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Local fish from Swiss lakes (frequent consumption of)

asso = fn_302_table_regression_by_covariate(cov_="ex1_swiss_fish",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]


#' ## Cosmetics 
#' 
#' ### Use of body lotion (frequency)

asso = fn_302_table_regression_by_covariate(cov_="ex1_bodylotion_freq",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Use of deodorant roll (frequency)

asso = fn_302_table_regression_by_covariate(cov_="ex1_deo_roll_freq",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Use of deodorant spray (frequency)

asso = fn_302_table_regression_by_covariate(cov_="ex1_deo_spray_freq",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Use of hair spray (frequency)

asso = fn_302_table_regression_by_covariate(cov_="ex1_hairspray_freq",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Use of perfume (frequency)

asso = fn_302_table_regression_by_covariate(cov_="ex1_perfume_freq",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Use of hand cream (frequency)

asso = fn_302_table_regression_by_covariate(cov_="ex1_handcream_freq",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ## Tobacco and alcohol
#' 
#' ### Consumption of smoked tobacco

asso = fn_302_table_regression_by_covariate(cov_="ex1_cig",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]


#' ### Consumption of other tobacco products

asso = fn_302_table_regression_by_covariate(cov_="ex1_other_tobacco",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Consumption of alcoholic products

asso = fn_302_table_regression_by_covariate(cov_="ex1_alcool_weekly",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]


#' ## Other personal exposures
#' 
#' ### Hot meals from disposable package (frequency)

asso = fn_302_table_regression_by_covariate(cov_="ex1_hot_meals_package_freq",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Use of ski fart

asso = fn_302_table_regression_by_covariate(cov_="ex2_skiwax",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Use of water-repellent spray

asso = fn_302_table_regression_by_covariate(cov_="ex2_spray_imp",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]




#' ## Professional expositions
#' 
#' ### Smoke

asso = fn_302_table_regression_by_covariate(cov_="ex2_prof_smoke",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Exhaust gas

asso = fn_302_table_regression_by_covariate(cov_="ex2_prof_exhaust_gas",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Solvent vapours

asso = fn_302_table_regression_by_covariate(cov_="ex2_prof_solvent_vapours",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Cleaning agents

asso = fn_302_table_regression_by_covariate(cov_="ex2_prof_cleaning_agents",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Dust

asso = fn_302_table_regression_by_covariate(cov_="ex2_prof_dust",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Ski wax

asso = fn_302_table_regression_by_covariate(cov_="ex2_prof_skiwax",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]

#' ### Impregnation sprays

asso = fn_302_table_regression_by_covariate(cov_="ex2_prof_spray_imp",coefs_=all_,dat_in=pei)
asso[[1]]
asso[[2]]





