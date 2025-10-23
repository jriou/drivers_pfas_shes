#' ---
#' title: Determinants of PFAS exposure 
#' subtitle: Data description
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
#' params:
#'   analysis_date:
#' bibliography: misc/bib.bib  
#' ---

#+ results="hide", warnings="false", echo="false"
if(is.null(controls)) controls = readRDS(file.path("savepoints/savepoint_",params$analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_1 = readRDS(file.path(controls$savepoint,"shesp_1.rds"))
shesp_2_stages = readRDS(file.path(controls$savepoint,"shesp_2.rds"))
shesp_2 = shesp_2_stages[[length(shesp_2_stages)]]
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds"))

#' # Exclusions
#' 
#' ## Stage 1
#' 
#' Data include four population types, but only participants from the random sample (`R`) are 
#' representative of the general population. There is also the convenience sample (`C`), the 
#' self-selected sample (`S`) and the participants recruited before the first COVID wave (`W`). 
#' 
shesp_2_stages$stage_0 %>%  
  dplyr::select(population,
                age,
                gender,
                center_label,
                date,
  ) %>% 
  frq_table(by="population",
            cap="Stage 1. Only random sample (`R`).")


#' ## Stage 2
#' 
#' Here we remove participants without PFAS measurements (`pfas_not_done=1`).
  
shesp_2_stages$stage_1 %>%  
  dplyr::select(pfas_not_done,
                age,
                gender,
                center_label,
                date,
  ) %>% 
  frq_table(by="pfas_not_done",
            cap="Stage 2. Remove participants without PFAS measurements (`pfas_not_done=1`).")


#' # Missing values
#' 
#' ## Missing values in socio-demographic characteristics
#' 

shesp_2 %>%
  dplyr::select(age,
                gender,
                a1_foreign,
                a1_education,
                a1_household_monthly_income,
                center_label
  ) %>%
  frq_table(cap="Socio-demographic characteristics of participants.")


#' ## Missing values in past SARS-CoV-2 immunity-conferring event
#' 

shesp_2  %>%  
  dplyr::select(corona_past_disease_confirmed,
                corona_past_disease_suspected,
                corona_past_disease_date,
                corona_past_hospit,
                corona_past_icu,
                corona_vaccine_sarscov2
  ) %>% 
  frq_table(cap="Variables related to SARS-CoV-2 immunity-conferring event.") 


#' # Description
#' 
#' We include a total of `r nrow(shesp_3)` participants. 
#' In the following, we describe all included variables and select which will be used in the analysis.
#' 
#' ## PFAS
#' 
#' We consider 32 different measurements, some of which are combinations of multiple substances. For each, we
#'  describe whether the substance is detected and, if it is detected, the quantification (mean, IQR and 95th percentile). 

fn_100_summary_pfas(shesp_3,cap="PFAS measurements.")

controls$pfas_allzero = c(2,3,10,11,18:27,29:32)
controls$pfas_list_retained = controls$pfas_list[-controls$pfas_allzero]
controls$pfas_labels_retained = controls$pfas_labels[-controls$pfas_allzero]

#' In `r length(controls$pfas_allzero)` cases, the substance was detected in fewer than 5 participants. We focus the 
#' analysis on the other `r length(controls$pfas_list_retained)` substances, each detected in at least 5 participants.

fn_100_summary_pfas(shesp_3,
                    cap="PFAS measurements (retained).",
                    exclusions=controls$pfas_allzero)

#' Accounting for the ponderation on age, sex and center, we find the following distributions.

fn_100_summary_pfas(shesp_3,
                    cap="PFAS measurements, after ponderation on age and sex (retained).",
                    exclusions=controls$pfas_allzero,
                    ponderation=TRUE)

#' We also look at the correlation between different PFASs.

fn_103_correlogram_pfas(shesp_3, 
                        exclusions=controls$pfas_allzero)


#' ## Socio-demographics 
#' 
#' The socio-demographic variables are shown in the following table. The distributions appear relatively balanced.
#' We create several new variables to reduce complexity:
#' 
#' - education is aggregated in 3 classes instead of 8 (considering `Other`as missing)
#' 
#' - income is aggregated in 3 classes instead of 6.
#' 
#' We consider all of these variables for the analysis.

shesp_3  %>%  
  dplyr::select(age,
                age_group,
                gender,
                a1_foreign,
                a1_education_3classes,
                a1_household_monthly_income_3classes,
                center_label
  ) %>% 
  frq_table(cap="Socio-demographic characteristics.",
            missing="ifany") 

controls$covariates_imputation = c("a1_education_3classes",
                                   "a1_household_monthly_income_3classes")


#' ## SARS-CoV-2 immunity
#' 
#' Anti-S IgG are associated with both natural infection and vaccination, while anti-N IgG are
#' only associated with natural infection. We will focus on anti-S IgG for the main analysis, as
#' the other measurements have a lot of missing values. More details about the laboratory
#' methods are available in [@fenwickHighthroughputCellVirusfree2021].

shesp_3  %>%  
  dplyr::select(starts_with("bu")) %>% 
  frq_table(cap="Laboratory measures of SARS-CoV-2 immunity.",
            missing="ifany") 

#' ## Potentially relevant exposures
#' 
#' ### Food frequency questionnaire
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
  dplyr::select(any_of(controls$covariates_ffq),cluster) %>% 
  frq_table(cap="Potential exposures from the food frequency questionnaire.",
            missing="ifany",
            labels=FALSE) 

#' In order to reduce the number of dimensions, and avoid colinearity issues, we used a clustering
#' algorithm to replace these 20 variables by a grouping.
#' We opted for the PAM (Partitioning Around Medoids) algorithm with the Manhattan distance.
#' We selected the number of groups based on the gap statistic, resulting in 5 clusters.
#' For clarity, we attributed a short name to each cluster based on the foodstuffs associated.

fn_104_ffq_heatmap(shesp_3,cluster_name="cluster")

shesp_3  %>%  
  dplyr::select(any_of(controls$covariates_ffq),cluster) %>% 
  frq_table(cap="Potential exposures from the food frequency questionnaire.",
            missing="ifany",
            by="cluster",
            labels=FALSE) 

shesp_3  %>%  
  dplyr::select(age,gender,a1_foreign,
                a1_education_3classes,
                a1_household_monthly_income_3classes,
                center_label,cluster) %>% 
  frq_table(cap="Socio-demographic characteristics of dietary clusters.",
            missing="ifany",
            by="cluster",
            labels=FALSE) 

#' ### Detailed fish consumption
#' 
#' To continue on foodstuffs, there are also more detailed questions about fish consumption that are
#' relevant to the issue.
#' 

controls$covariates_fish = c("ex1_fish_freq","ex1_tuna","ex1_salmon","ex1_pangasius","ex1_trout","ex1_shrimp","ex1_organic_fish","ex1_swiss_fish")

shesp_3  %>%  
  dplyr::select(any_of(controls$covariates_fish)) %>% 
  frq_table(cap="Potentially relevant exposures.",
            missing="ifany",
            labels=FALSE) 

#' ### Other potential exposures
#' 

shesp_3  %>%  
  dplyr::select(starts_with("scre"),starts_with("ex1"),starts_with("ex2"),starts_with("hm"),starts_with("ffq")) %>% 
  dplyr::select(-any_of(controls$covariates_ffq),-any_of(controls$covariates_fish)) %>% 
  frq_table(cap="Potentially relevant exposures.",
            missing="ifany",
            labels=FALSE) 


controls$covariates_list = c("age",
                             "age_group",
                             "gender",
                             "a1_foreign",
                             "a1_education_3classes",
                             "a1_household_monthly_income_3classes",
                             "center_label",
                             
                             "ex1_residence_urban",
                             "ex2_skiwax",
                             "ex2_spray_imp",
                             "ex1_carpet",
                             "ex1_hot_meals_package_freq",
                             
                             "ex1_alcool_weekly",
                             "ex1_cig",
                             
                             "cluster",
                             "ex1_veggie_vegan",
                             "ex1_fish_freq",
                             "ex1_tuna",
                             "ex1_salmon",
                             "ex1_pangasius",
                             "ex1_trout",
                             "ex1_shrimp",
                             "ex1_swiss_fish",
                             "ex1_tapwater",

                             "ex1_deo_spray_freq",
                             "ex1_deo_roll_freq",
                             "ex1_bodylotion_freq",
                             "ex1_handcream_freq",
                             "ex1_perfume_freq",
                             "ex1_hairspray_freq",

                             "ex2_prof_smoke",
                             "ex2_prof_exhaust_gas",
                             "ex2_prof_solvent_vapours",
                             "ex2_prof_cleaning_agents",
                             "ex2_prof_dust",
                             "ex2_prof_skiwax",
                             "ex2_prof_spray_imp"
)

controls$covariates_imputation = c("a1_education_3classes",
                                   "a1_household_monthly_income_3classes",
                                   "ex2_garden_pesticides",
                                   "ex2_skiwax",
                                   "ex1_alcool_weekly",
                                   "ex1_veggie_vegan",
                                   "ex1_residence_urban",
                                   "ex2_prof_smoke",
                                   "ex2_prof_exhaust_gas",
                                   "ex2_prof_solvent_vapours",
                                   "ex2_prof_cleaning_agents",
                                   "ex2_prof_dust	",
                                   "ex2_prof_skiwax",
                                   "ex2_prof_spray_imp")


#' ## Past SARS-CoV-2 infection
#' 
#' About 10% of participants had a confirmed SARS-CoV-2 infection, this rises to 14% if we include
#' suspected infections. There are very few hospitalizations and ICU stays linked to SARS-CoV-2, 
#' so we don't have a good way to assess the severity of the infection. Vaccination is about as
#' expected with regards to the general population. We retain as potential covariates past 
#' confirmed or suspected infection and vaccination (the last one will require imputation).

shesp_3  %>%  
  dplyr::select(corona_past_disease_confirmed,
                corona_past_disease_suspected,
                corona_past_hospit,
                corona_past_icu,
                corona_vaccine_sarscov2
  ) %>% 
  frq_table(cap="Past SARS-CoV-2 infection.",
            missing="ifany") 



controls$covariates_list2 = c(controls$covariates_list,
                              "corona_past_disease_suspected",
                              "corona_vaccine_sarscov2")
controls$covariates_imputation2 = c(controls$covariates_imputation,
                                    "corona_vaccine_sarscov2")

