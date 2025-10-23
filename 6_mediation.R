#' ---
#' title: Determinants of PFAS exposure 
#' subtitle: Mediation analysis
#' author: jriou
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
#' bibliography: misc/bib.bib  
#' ---

#+ results="hide", warnings="false", echo="false"
analysis_date = "2025-05-23"
controls = readRDS(paste0("savepoints/savepoint_",analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds"))
pei = shesp_3 %>% 
  dplyr::select(did,
                pfas_cluster,
                all_of(paste0(controls$pfas_list_retained,"_nonzero")),
                all_of(paste0(controls$pfas_list_retained,"_quant")),
                all_of(controls$covariates_list))  %>% 
  dplyr::mutate(y=as.integer(pfas_cluster %in% c("2_Intermediate","3_High")),
                .after = pfas_cluster)

controls$lab_pfas_cluster <- c( "#1 (low)", "#2 (intermediate)","#3 (high)")
controls$col_pfas_cluster <- c( "#46B8DA", "#5CB85C","#EEA236")

# apply multiple imputation
predmat = mice::quickpred(pei,
                          include = c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label"), # multiple imputation based on variables listed in `impute`
                          mincor = 1, # exclude all other covariates
                          exclude = "y")
dat_imp = mice::mice(pei, 
                     m = 5, 
                     predictorMatrix = predmat,
                     print = FALSE)
dat_imp1 = mice::complete(dat_imp)

# formula
fn_208_binom = function(dat_,covs,...) {
  prior_ = prior("normal(0,1)", class="b")
  rhs_ = covs
  rhs_ = paste(rhs_[!duplicated(rhs_)], collapse=" + ")
  form_ = bf(paste0("y ~ ", rhs_))
  print(form_)
  mod_ = brms::brm(form_,
                   family = bernoulli(link = "logit"),
                   data = dat_imp, 
                   prior = prior_,
                   chains = 4,
                   cores = 4,
                   iter = 2000,
                   warmup = 1000,
                   control = list(adapt_delta = 0.95))
  return(mod_)
}

# R squared 
r2_out = NULL

## Age and sex
cov1.1 = c("age","gender")
mod1.1 = fn_208_binom(dat_imp1,cov1.1)
r2_out = bayes_R2(mod1.1) %>% as_tibble() %>% mutate(model="1.1",.before="Estimate") %>% bind_rows(r2_out)

## Sociodemo without center
cov1.2 = c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes")
mod1.2 = fn_208_binom(dat_imp1,cov1.2)
r2_out = bayes_R2(mod1.2) %>% as_tibble() %>% mutate(model="1.2",.before="Estimate") %>% bind_rows(r2_out)

## Sociodemo with center
cov1.3 = c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label")
mod1.3 = fn_208_binom(dat_imp1,cov1.3)
r2_out = bayes_R2(mod1.3) %>% as_tibble() %>% mutate(model="1.3",.before="Estimate") %>% bind_rows(r2_out)

## Top 10 
cov1.4 = c("age","gender","cluster","center_label",            
           "a1_foreign","ex1_handcream_freq","ex1_fish_freq","ex2_skiwax" ,             
           "ex2_prof_solvent_vapours","ex1_bodylotion_freq" )
mod1.4 = fn_208_binom(dat_imp1,cov1.4)
r2_out = bayes_R2(mod1.4) %>% as_tibble() %>% mutate(model="1.4",.before="Estimate") %>% bind_rows(r2_out)

## Sociodemo plus top 10
cov1.5 = c("age","gender","a1_foreign","cluster","a1_education_3classes","a1_household_monthly_income_3classes","center_label",
           "ex1_handcream_freq","ex1_fish_freq","ex2_skiwax" ,             
           "ex2_prof_solvent_vapours","ex1_bodylotion_freq")
mod1.5 = fn_208_binom(dat_imp1,cov1.5)
r2_out = bayes_R2(mod1.5) %>% as_tibble() %>% mutate(model="1.5",.before="Estimate") %>% bind_rows(r2_out)

## Age and sex
cov1.6 = c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label","cluster")
mod1.6 = fn_208_binom(dat_imp1,cov1.6)
r2_out = bayes_R2(mod1.6) %>% as_tibble() %>% mutate(model="1.6",.before="Estimate") %>% bind_rows(r2_out)


save(r2_out,mod1.1,mod1.2,mod1.3,mod1.4,mod1.5,mod1.6,file=file.path(controls$savepoint,"6_mod_rsquared.Rdata"))

# Mediation

library(mediation)
dat_imp2 = dat_imp1 %>% 
  as_tibble() %>% 
  dplyr::select(-did,-pfas_cluster,-starts_with("pfas")) %>% 
  # collapse categorical variables into binary
  dplyr::mutate(cluster=if_else(cluster %in% c("2_dairy_focused","4_plant_based"),1,0)) %>% 
  dplyr::mutate(ex1_alcool_weekly=if_else(cluster %in% c("At least weekly","Daily and more"),1,0)) %>% 
  dplyr::mutate(ex1_cig=if_else(cluster %in% c("Past smoker","Current smoker"),1,0)) %>% 
  dplyr::mutate(a1_education_3classes=if_else(a1_education_3classes %in% c("University","Apprenticeship / professional baccalaureate"),1,0)) %>% 
  dplyr::mutate(a1_household_monthly_income_3classes=if_else(a1_household_monthly_income_3classes %in% c("CHF 4,500 to 9,000",">CHF 9,000"),1,0)) %>% 
  # transform all into 0/1
  dplyr::mutate(across(where(is.character),~as.factor(.x))) %>% 
  dplyr::mutate(across(where(is.factor),~as.numeric(.x)-1)) %>%
  dplyr::select(-y,-age_group) 
dat_imp2b = bind_cols(y=dat_imp1$y,dat_imp2)

varlist = controls$covariates_list[8:37]

## center_label

### BART selection
ovar_ = "center_label"
sel_ = fn_205_bart_selection_supp(dat_in=dat_imp2, dependent=ovar_, covariates = varlist)
med_ = fn_208_mediation(dat_in=dat_imp2b, ovar_, sel_) 
med_
sum(med_$proportion_mediated)
out_ = med_

## gender

### BART selection
ovar_ = "gender"
sel_ = fn_205_bart_selection_supp(dat_in=dat_imp2, dependent=ovar_, covariates = varlist)
med_ = fn_208_mediation(dat_in=dat_imp2b, ovar_, sel_) 
med_
sum(med_$proportion_mediated)
out_ = bind_rows(out_,med_)

## age

### BART selection
ovar_ = "age"
sel_ = fn_205_bart_selection_supp(dat_in=dat_imp2, dependent=ovar_, covariates = varlist)
med_ = fn_208_mediation(dat_in=dat_imp2b, ovar_, sel_) 
med_
sum(med_$proportion_mediated)
out_ = bind_rows(out_,med_)

## foreign

### BART selection
ovar_ = "a1_foreign"
sel_ = fn_205_bart_selection_supp(dat_in=dat_imp2, dependent=ovar_, covariates = varlist)
med_ = fn_208_mediation(dat_in=dat_imp2b, ovar_, sel_) 
med_
sum(med_$proportion_mediated)
out_ = bind_rows(out_,med_)

## income

### BART selection
ovar_ = "a1_household_monthly_income_3classes"
sel_ = fn_205_bart_selection_supp(dat_in=dat_imp2, dependent=ovar_, covariates = varlist)
med_ = fn_208_mediation(dat_in=dat_imp2b, ovar_, sel_) 
med_
sum(med_$proportion_mediated)
out_ = bind_rows(out_,med_)


## education

### BART selection
ovar_ = "a1_education_3classes"
sel_ = fn_205_bart_selection_supp(dat_in=dat_imp2, dependent=ovar_, covariates = varlist)
med_ = fn_208_mediation(dat_in=dat_imp2b, ovar_, sel_) 
med_
sum(med_$proportion_mediated)
out_ = bind_rows(out_,med_)


saveRDS(out_,file=file.path(controls$savepoint,"6_mediation.rds"))
