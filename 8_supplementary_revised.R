#' ---
#' title: "Diet, lifestyle, and sociodemographic influences on serum concentration of per- and polyfluoroalkyl substances (PFASs): insights from human biomonitoring in Switzerland"
#' subtitle: "Supplementary Information — Appendix Text S1"
#' author: "Riou et al."
#' date: "9 October 2025"
#' output:
#'   pdf_document: 
#'     keep_tex: yes
#' ---

#+ results="hide", warnings="false", echo="false"
analysis_date = "2026-02-23"
controls = readRDS(paste0("savepoints/savepoint_",analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds")) %>% 
  mutate(cluster = case_when(cluster=="1_balanced" ~ "Traditional",
                             cluster=="2_dairy_focused" ~ "Dairy-focused",
                             cluster=="3_high_fiber" ~ "Mediterranean",
                             cluster=="4_plant_based" ~ "Plant-based",
                             cluster=="5_meat_centered" ~ "Western"))

controls$lab_pfas_cluster <- c( "#1 (low, 30%)", "#2 (intermediate, 37%)","#3 (high, 33%)")
controls$col_pfas_cluster <- c( "#46B8DA", "#5CB85C","#EEA236")


controls$covariates_list_women <- controls$covariates_list[-3]
controls$covariates_list_women <- c(controls$covariates_list_women,"hm_kids","hm_breastfeeding","hm_menopause")
# keep only relevant variables and remove men
shesp_3_women = shesp_3  %>% 
  dplyr::filter(gender=="Females")


#' 
#' # Additional methods 
#'  
#' ## Diet cluster
#' 
#' The diet clusters were generated based on 20 variables from the food frequency questionnaire (available in the Supplementary Text S2) using the Partitioning Around Medoids (PAM) method.  Specifically, a single‑imputation run of the Multivariate Imputation by Chained Equations (MICE) algorithm was performed to replace 1 missing answer on whole grain consumption. We then computed the pairwise distances between participants using the Manhattan distance, which is appropriate for ordinal frequency scales and less sensitive to outliers than the Euclidean distance. We determined  the optimal number of clusters ($k$) using the gap statistic. In the final analysis the number of clusters was fixed to five (k = 5).

# knitr::include_graphics(file.path(controls$savepoint,"best_k_plot_gap_distman.pdf"))

#' **Figure.** Optimal number of diet clusters determined by the gap statistic using Manhattan distance. 
#' 
#' ## PFAS cluster
#' 
#' A similar approach was applied to PFAS serum concentrations. The concentrations were first log-transformed, after replacing undetected values by 0.01, then standardized. We then applied the PAM algorithm with a Manhattan distance, and determined the number of clusters based on the gap statistic. In that case, we fixed the number of clusters to three (k=3).

# knitr::include_graphics(file.path(controls$savepoint,"pfas_best_k_plot_gap_distman.pdf"))

#' **Figure.** Optimal number of PFAS clusters determined by the gap statistic using Manhattan distance. 
#' 
#' 
#' # Additional description
#' 
#' ## Socio-demographic characteristics
#' 

shesp_3 %>%
  dplyr::select(age,
                age_group,
                gender,
                a1_foreign,
                a1_education_3classes,
                a1_household_monthly_income_3classes,
                center_label
  ) %>%
  frq_table(cap="Socio-demographic characteristics of participants.",
            missing="ifany",
            labels = TRUE)

#' 
#' ## Food frequency questionnaire
#' 
#' The food frequency questionnaire consists of estimates of frequency for 20 foodstuffs.
#' Each foodstuff is assigned a number from 1 to 8 :
#' 
#' -  1 - Rarely / Never;
#' -  2 - Once a month;
#' -  3 - Every 2 weeks; 
#' -  4 - 1-2 times a week; 
#' -  5 - 3-6 times a week;
#' -  6 - Once a day;
#' -  7 - 2-3 times a day; 
#' -  8 - 4 times a day and more.

shesp_3  %>%  
  dplyr::select(any_of(controls$covariates_ffq)) %>% 
  frq_table(cap="Potential exposures from the food frequency questionnaire.",
            missing="ifany",
            labels=TRUE) 

#' ## Diet cluster
#' 

shesp_3  %>%  
  dplyr::select(age,gender,a1_foreign,
                a1_education_3classes,
                a1_household_monthly_income_3classes,
                center_label,cluster) %>% 
  frq_table(cap="Socio-demographic characteristics of dietary clusters.",
            missing="ifany",
            by="cluster",
            labels=TRUE) 

#' 
#' ## Other explanatory variables
#' 

shesp_3  %>%  
  dplyr::select(any_of(controls$covariates_list[-c(1:7)])) %>% 
  frq_table(cap="Other characteristics of participants.",
            missing="ifany",
            labels = TRUE) 

#' ## PFAS
#' 
#' We consider 32 different measurements, some of which are combinations of multiple substances. For each, we
#' describe whether the substance is detected and, if it is detected, the quantification (mean, IQR and 95th percentile). 

fn_100_summary_pfas(shesp_3,cap="PFAS measurements.")

#' In `r length(controls$pfas_allzero)` cases, the substance was detected in fewer than 5 participants. We focus the 
#' analysis on the other `r length(controls$pfas_list_retained)` substances, each detected in at least 5 participants.

fn_100_summary_pfas(shesp_3,
                    cap="PFAS measurements (retained).",
                    exclusions=controls$pfas_allzero)

# #' Accounting for the ponderation on age, sex and center, we find the following distributions.
# 
# fn_100_summary_pfas(shesp_3,
#                     cap="PFAS measurements, after ponderation on age and sex (retained).",
#                     exclusions=controls$pfas_allzero,
#                     ponderation=TRUE) 

#' By PFAS cluster 

fn_100_summary_pfas(shesp_3,
                    cap="PFAS measurements, by PFAS cluster.",
                    exclusions=controls$pfas_allzero,
                    groupby="pfas_cluster")

#' We also look at the correlation between different PFASs.

fn_103_correlogram_pfas(shesp_3, 
                        exclusions=controls$pfas_allzero)


#'  
#' # Variables associated with PFAS cluster
#' 

uni_ = readRDS(file.path(controls$savepoint,"4_univariable_pfas_cluster.rds")) %>% filter(!grepl("age_group",parameter))
mul_ = readRDS(file.path(controls$savepoint,"4_multivariable_pfas_cluster.rds"))  %>% filter(!grepl("age_group",parameter))
mul2_ = readRDS(file.path(controls$savepoint,"4_multivariable2_pfas_cluster.rds")) %>% filter(!grepl("age_group",parameter))
bmul_ = readRDS(file.path(controls$savepoint,"4_bart_multivariable_pfas_cluster.rds")) %>% 
  dplyr::mutate(adjusted="Adjusted, BART selection") %>% 
  dplyr::filter(parameter!="Intercept")

#' ## Univariable
#' 

fn_302_table_regression_by_covariate_cluster(cov_="ALL",coefs_=uni_)

#' ## Adjusted on age and gender
#' 

fn_302_table_regression_by_covariate_cluster(cov_="ALL",coefs_=mul_)

#' ## Adjusted on all socio-demographic variables
#' 

fn_302_table_regression_by_covariate_cluster(cov_="ALL",coefs_=mul2_)

#' ## Adjusted on top 10 variables selected by BART
#' 

fn_302_table_regression_by_covariate_cluster(cov_="ALL",coefs_=bmul_)



#'  
#' # Variables associated with each PFAS
#' 

uni_ = readRDS(file.path(controls$savepoint,"2_univariable_pfas.rds")) %>% filter(!grepl("age_group",parameter))
mul_ = readRDS(file.path(controls$savepoint,"2_multivariable_pfas.rds"))  %>% filter(!grepl("age_group",parameter)) %>% mutate(adjusted="Adjusted, simple")
mul2_ = readRDS(file.path(controls$savepoint,"2_multivariable2_pfas.rds")) %>% filter(!grepl("age_group",parameter)) %>% mutate(adjusted="Adjusted, extended")
mul2_ = readRDS(file.path(controls$savepoint,"2_multivariable2_pfas.rds")) %>% filter(!grepl("age_group",parameter)) %>% mutate(adjusted="Adjusted, extended")
bmul_ = readRDS(file.path(controls$savepoint,"2_bart_multivariable_pfas.rds")) %>% 
  dplyr::mutate(adjusted="Adjusted, BART selection")
all_ = dplyr::bind_rows(uni_,mul_,mul2_,bmul_)


#' ## Socio-demographic
#' 
#' ### Age 

asso = fn_302_table_regression_by_covariate3(cov_="age",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Sex

asso = fn_302_table_regression_by_covariate3(cov_="gender",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Foreign nationality

asso = fn_302_table_regression_by_covariate3(cov_="a1_foreign",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Education

asso = fn_302_table_regression_by_covariate3(cov_="a1_education_3classes",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Household monthly income

asso = fn_302_table_regression_by_covariate3(cov_="a1_household_monthly_income_3classes",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Study site

asso = fn_302_table_regression_by_covariate3(cov_="center_label",coefs_=all_)
asso[[1]]
asso[[2]]

#' ## Food
#' 

#' ### Diet cluster
#' 
asso = fn_302_table_regression_by_covariate3(cov_="cluster",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Overall fish consumption
#' 
asso = fn_302_table_regression_by_covariate3(cov_="ex1_fish_freq",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Tuna (frequent consumption of)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_tuna",coefs_=all_)
asso[[1]]
asso[[2]]

#' ###  Salmon (frequent consumption of)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_salmon",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Pangasius (frequent consumption of)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_pangasius",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Trout (frequent consumption of)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_trout",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Shrimp (frequent consumption of)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_shrimp",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Local fish from Swiss lakes (frequent consumption of)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_swiss_fish",coefs_=all_)
asso[[1]]
asso[[2]]

#' ## Diet and lifestyle
#' 
#' ### Vegetarian or vegan

asso = fn_302_table_regression_by_covariate3(cov_="ex1_veggie_vegan",coefs_=all_)
asso[[1]]
asso[[2]]

#' ###  Tap water

asso = fn_302_table_regression_by_covariate3(cov_="ex1_tapwater",coefs_=all_)
asso[[1]]
asso[[2]]

#' ###  Alcohol

asso = fn_302_table_regression_by_covariate3(cov_="ex1_alcool_weekly",coefs_=all_)
asso[[1]]
asso[[2]]

#' ###  Tobacco

asso = fn_302_table_regression_by_covariate3(cov_="ex1_cig",coefs_=all_)
asso[[1]]
asso[[2]]

#' ## Environmental
#' 
#' ### Use of ski fart

asso = fn_302_table_regression_by_covariate3(cov_="ex2_skiwax",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Use of water-repellent spray

asso = fn_302_table_regression_by_covariate3(cov_="ex2_spray_imp",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Hot meals from disposable package (frequency)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_hot_meals_package_freq",coefs_=all_)
asso[[1]]
asso[[2]]


#' ### Residence in urban or industrial area compared to rural area

asso = fn_302_table_regression_by_covariate3(cov_="ex1_residence_urban",coefs_=all_,cov_dataname="ex1_residence_urban")
asso[[1]]
asso[[2]]

#' ### Synthetic or natural carpet at home 

asso = fn_302_table_regression_by_covariate3(cov_="ex1_carpet",coefs_=all_)
asso[[1]]
asso[[2]]



#' ## Cosmetics 
#' 
#' ### Use of body lotion (frequency)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_bodylotion_freq",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Use of deodorant roll (frequency)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_deo_roll_freq",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Use of deodorant spray (frequency)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_deo_spray_freq",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Use of hair spray (frequency)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_hairspray_freq",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Use of perfume (frequency)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_perfume_freq",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Use of hand cream (frequency)

asso = fn_302_table_regression_by_covariate3(cov_="ex1_handcream_freq",coefs_=all_)
asso[[1]]
asso[[2]]



#' ## Professional expositions
#' 
#' ### Smoke

asso = fn_302_table_regression_by_covariate3(cov_="ex2_prof_smoke",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Exhaust gas

asso = fn_302_table_regression_by_covariate3(cov_="ex2_prof_exhaust_gas",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Solvent vapours

asso = fn_302_table_regression_by_covariate3(cov_="ex2_prof_solvent_vapours",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Cleaning agents

asso = fn_302_table_regression_by_covariate3(cov_="ex2_prof_cleaning_agents",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Dust

asso = fn_302_table_regression_by_covariate3(cov_="ex2_prof_dust",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Ski wax

asso = fn_302_table_regression_by_covariate3(cov_="ex2_prof_skiwax",coefs_=all_)
asso[[1]]
asso[[2]]

#' ### Impregnation sprays

asso = fn_302_table_regression_by_covariate3(cov_="ex2_prof_spray_imp",coefs_=all_)
asso[[1]]
asso[[2]]


#'  
#' # Sensitivity analysis focused on women
#' 
#' ## Description
#' 
#' ### Socio-demographic characteristics
#' 

shesp_3_women %>%
  dplyr::select(age,
                age_group,
                a1_foreign,
                a1_education_3classes,
                a1_household_monthly_income_3classes,
                center_label
  ) %>%
  frq_table(cap="Socio-demographic characteristics of female participants.",
            missing="ifany",
            labels = TRUE)

#' ### Diet cluster
#' 

shesp_3_women  %>%  
  dplyr::select(age,a1_foreign,
                a1_education_3classes,
                a1_household_monthly_income_3classes,
                center_label,cluster) %>% 
  frq_table(cap="Socio-demographic characteristics of dietary clusters.",
            missing="ifany",
            by="cluster",
            labels=TRUE) 

#' 
#' ### Other explanatory variables
#' 

shesp_3_women  %>%  
  dplyr::select(any_of(controls$covariates_list_women[-c(1:7)])) %>% 
  frq_table(cap="Other characteristics of participants.",
            missing="ifany",
            labels = TRUE) 

#' ### PFAS
#' 

fn_100_summary_pfas(shesp_3_women,
                    cap="PFAS measurements (retained).",
                    exclusions=controls$pfas_allzero)

#'  
#' ## Variables associated with PFAS cluster
#' 

uni_ = readRDS(file.path(controls$savepoint,"4_univariable_pfas_cluster_women.rds")) %>% filter(!grepl("age_group",parameter))
mul_ = readRDS(file.path(controls$savepoint,"4_multivariable_pfas_cluster_women.rds"))  %>% filter(!grepl("age_group",parameter))
mul2_ = readRDS(file.path(controls$savepoint,"4_multivariable2_pfas_cluster_women.rds")) %>% filter(!grepl("age_group",parameter))

controls$lab_pfas_cluster <- c( "#1 (low, 43%)", "#2 (intermediate, 40%)","#3 (high, 17%)")

#' ## Univariable
#' 

fn_302_table_regression_by_covariate_cluster(cov_="ALL",coefs_=uni_)

#' ## Adjusted on age and gender
#' 

fn_302_table_regression_by_covariate_cluster(cov_="ALL",coefs_=mul_)

#' ## Adjusted on all socio-demographic variables
#' 

fn_302_table_regression_by_covariate_cluster(cov_="ALL",coefs_=mul2_)
