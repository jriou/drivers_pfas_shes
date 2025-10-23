#' ---
#' title: Determinants of PFAS exposure 
#' subtitle: Data analysis 1
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
analysis_date = "2025-09-04"
controls = readRDS(paste0("savepoints/savepoint_",analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds"))

# keep only relevant variables
pei = shesp_3 %>% 
  dplyr::select(did,
                all_of(paste0(controls$pfas_list_retained,"_nonzero")),
                all_of(paste0(controls$pfas_list_retained,"_quant")),
                all_of(controls$covariates_list)) 

#' # Methods
#' 
#' Methods (to be detailed):
#' 
#' - representative sample of `r nrow(pei)` participants from the SHeS-pilot study
#' 
#' - `r length(controls$pfas_list_retained)` different PFAS measurements, some of which are combinations of others, and 
#' each detected in at least 5 participants
#' 
#' - `r length(controls$covariates_list)` potentially relevant covariates based on the literature
#' 
#' - specificity of measurement (detection + quantification) suggests hurdle gamma regression when detection only concerns a subsample of participants,
#' and gamma regression when the substance is detected in all participants.
#' 
#' - outputs of hurdle gamma model $\exp(\beta)$: odds ratio of detection AND relative change in measure (only relative change in measure for gamma regression)
#' 
#' - first univariable (for each PFAS), then multivariable 
#' 
#' - Bayesian inference with `brms` package
#' 
#' - weakly informative prior distributions (Gelman et al 2008) on parameters: normal(0,5) on the log or logit scale
#' 
#' - because of the large number of covariates, we use penalized regression with horseshoe priors (Piironen and Vehtari, 2017). 
#' This parametrization pulls all estimates towards 1 (no effect), while allowing an expected 10% of covariates to be away from 1.
#' 
#' - in several cases we impose a monotonic constraint on the effect of covariates, so that higher exposure translates into a larger effect (https://cran.r-project.org/web/packages/brms/vignettes/brms_monotonic.html)
#' 
#' - imputation of missing data within `brms` (https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html)
#' 

#' # Description
#' 
# 
# fn_100_summary_pfas(dat_in=pei,
#                     cap="PFAS measurements.",
#                     exclusions=controls$pfas_allzero)
# 
# fn_102_hist_pfas(pei,
#                  var="pfas_pfda",
#                  cap="PFDA [335-76-2]",
#                  lower_limit=5)
# 
# pei  %>%
#   dplyr::select(-did,-starts_with("pfas")) %>%
#   frq_table(cap="Covariates.",missing="ifany")



#' # Results by PFAS
#' 
#' ## Univariable regression
#' 

if(TRUE) {
  uni_ = list()
  for(i in 1:length(controls$pfas_list_retained)) {
    fn_010_clean_dlls() # removes temporary files
    dep_ = controls$pfas_list_retained[i]
    if(sum(pei[,paste0(dep_,"_nonzero")]==0)>0) {
      uni_[[i]] = fn_200_hurdle_gamma(dat_in=pei, 
                                      dependent=dep_, 
                                      covariates=controls$covariates_list, 
                                      multiple_runs=TRUE,
                                      multivariable=FALSE,
                                      impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label")
      )
    } else {
      uni_[[i]] = fn_201_gamma(pei, 
                               dependent=dep_, 
                               covariates=controls$covariates_list, 
                               multiple_runs=TRUE,
                               multivariable=FALSE,
                               impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label")
      )
    }
    write_rds(uni_,file=file.path(controls$savepoint,"2_univariable_pfas_tmp.rds"))
  }
  uni_ = uni_ %>% 
    do.call("rbind",.) 
  write_rds(uni_,file=file.path(controls$savepoint,"2_univariable_pfas.rds"))
}
uni_ = readRDS(file.path(controls$savepoint,"2_univariable_pfas.rds"))

#+ fig.width=11, fig.height=18
# fn_300_fig_regression(uni_)

#' ## Multivariable regression adjusted on age, sex
#' 

if(TRUE) {
  mul_ = list()  
  controls$covariates_list = controls$covariates_list[!controls$covariates_list %in% c("age_group")]
  for(i in 1:length(controls$pfas_list_retained)) {
    # suppressWarnings(fn_010_clean_dlls())
    dep_ = controls$pfas_list_retained[i]
    if(sum(pei[,paste0(dep_,"_nonzero")]==0)>0) {
      mul_[[i]] = fn_200_hurdle_gamma(pei, 
                                      dependent=dep_, 
                                      covariates=controls$covariates_list, 
                                      multivariable=TRUE, 
                                      multiple_runs=TRUE,
                                      regularization=FALSE,
                                      impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                                      adj=c("age","gender"))
    } else {
      mul_[[i]] = fn_201_gamma(pei, 
                               dependent=dep_, 
                               covariates=controls$covariates_list, 
                               multivariable=TRUE, 
                               multiple_runs=TRUE,
                               regularization=FALSE,
                               impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                               adj=c("age","gender"))
    }
    write_rds(mul_,file=file.path(controls$savepoint,"2_multivariable_pfas_tmp.rds"))
  }
  mul_ = mul_ %>% 
    do.call("rbind",.) %>% 
    mutate(adjusted="Adjusted")
  write_rds(mul_,file=file.path(controls$savepoint,"2_multivariable_pfas.rds"))
}
# mul_ = readRDS(file.path(controls$savepoint,"2_multivariable_pfas.rds")) %>% 
  # dplyr::mutate(adjusted="Adjusted, simple")

#+ fig.width=11, fig.height=18
# fn_300_fig_regression(mul_)


#' ## Multivariable regression adjusted on age, sex, nationality, income, education, center
#' 

if(TRUE) {
  mul2_ = list()
  for(i in 1:length(controls$pfas_list_retained)) {
    # fn_010_clean_dlls()
    dep_ = controls$pfas_list_retained[i]
    if(sum(pei[,paste0(dep_,"_nonzero")]==0)>0) {
      mul2_[[i]] = fn_200_hurdle_gamma(pei, 
                                       dependent=dep_, 
                                       covariates=controls$covariates_list, 
                                       multivariable=TRUE, 
                                       multiple_runs=TRUE,
                                       regularization=FALSE,
                                       impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                                       adj=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
    } else {
      mul2_[[i]] = fn_201_gamma(pei, 
                                dependent=dep_, 
                                covariates=controls$covariates_list, 
                                multivariable=TRUE, 
                                multiple_runs=TRUE,
                                regularization=FALSE,
                                impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                                adj=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
    }
    write_rds(mul2_,file=file.path(controls$savepoint,"2_multivariable2_pfas_tmp.rds"))
  }
  mul2_ = mul2_ %>% 
    do.call("rbind",.) %>% 
    mutate(adjusted="Adjusted")
  write_rds(mul2_,file=file.path(controls$savepoint,"2_multivariable2_pfas.rds"))
}
mul2_ = readRDS(file.path(controls$savepoint,"2_multivariable2_pfas.rds")) %>% 
  dplyr::mutate(adjusted="Adjusted, extended")

#+ fig.width=11, fig.height=18
# fn_300_fig_regression(mul2_)

#' ## Multivariable regression based on top 10 from BART variable selection
#' 

if(TRUE) {
  bmul_ = list()
  controls$covariates_list = controls$covariates_list[!controls$covariates_list %in% c("age_group")]
  for(i in 1:length(controls$pfas_list_retained)) {
    cat("\n",i,"---------------------------- \n")
    # suppressWarnings(fn_010_clean_dlls())
    dep_ = controls$pfas_list_retained[i]
    sel_ = fn_205_bart_selection(dat_in=pei,
                                 dependent=dep_,
                                 covariates=controls$covariates_list,
                                 impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
    sel_ = names(sel_$var_importance[1:10])
    if(sum(pei[,paste0(dep_,"_nonzero")]==0)>0) {
      bmul_[[i]] = fn_200_hurdle_gamma(pei, 
                                       dependent=dep_, 
                                       covariates=sel_, 
                                       multivariable=TRUE, 
                                       multiple_runs=TRUE,
                                       regularization=FALSE,
                                       adj=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                                       impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
    } else {
      bmul_[[i]] = fn_201_gamma(pei, 
                                dependent=dep_, 
                                covariates=sel_, 
                                multivariable=TRUE, 
                                multiple_runs=TRUE,
                                adj=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"),
                                regularization=FALSE,
                                impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
    }
    write_rds(bmul_,file=file.path(controls$savepoint,"2_bart_multivariable_pfas_tmp.rds"))
  }
  bmul_ = bmul_ %>% 
    do.call("rbind",.) %>% 
    mutate(adjusted="Adjusted penalized")
  write_rds(bmul_,file=file.path(controls$savepoint,"2_bart_multivariable_pfas.rds"))
}
bmul_ = readRDS(file.path(controls$savepoint,"2_bart_multivariable_pfas.rds")) %>% 
  dplyr::mutate(adjusted="Adjusted, BART selection")



#+ fig.width=11, fig.height=18
# fn_300_fig_regression(bmul_)

#' ## Multivariable penalized regression
#' 

if(TRUE) {
  pmul_ = list()
  controls$covariates_list = controls$covariates_list[!controls$covariates_list %in% c("age_group")]
  for(i in 1:length(controls$pfas_list_retained)) {
    cat("\n",i,"---------------------------- \n")
    # suppressWarnings(fn_010_clean_dlls())
    dep_ = controls$pfas_list_retained[i]
    if(sum(pei[,paste0(dep_,"_nonzero")]==0)>0) {
      pmul_[[i]] = fn_200_hurdle_gamma(pei, 
                                       dependent=dep_, 
                                       covariates=controls$covariates_list, 
                                       multivariable=TRUE, 
                                       multiple_runs=FALSE,
                                       regularization=TRUE,
                                       impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
    } else {
      pmul_[[i]] = fn_201_gamma(pei, 
                                dependent=dep_, 
                                covariates=controls$covariates_list, 
                                multivariable=TRUE, 
                                multiple_runs=FALSE,
                                regularization=TRUE,
                                impute=c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"))
    }
    write_rds(pmul_,file=file.path(controls$savepoint,"2_pen_multivariable_pfas_tmp.rds"))
  }
  pmul_ = pmul_ %>% 
    do.call("rbind",.) %>% 
    mutate(adjusted="Adjusted penalized")
  write_rds(pmul_,file=file.path(controls$savepoint,"2_pen_multivariable_pfas.rds"))
}
pmul_ = readRDS(file.path(controls$savepoint,"2_pen_multivariable_pfas.rds"))  %>% 
  dplyr::mutate(adjusted="Adjusted, penalization")

#+ fig.width=11, fig.height=18
fn_300_fig_regression(pmul_)

#' # Results by covariate
#' 

all_ = dplyr::bind_rows(uni_,mul_,mul2_,bmul_,pmul_)

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





