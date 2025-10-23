#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 302 - table for regression coefficients
#' author: Julien Riou
#' date: 2023-11-23
#' ---

fn_302_table_regression_by_covariate = function(cov_,coefs_,summarise_results=FALSE, cov_dataname=NULL, aOR1=TRUE, aOR2=TRUE, ...) {
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    # manual changes
    dplyr::mutate(parameter2_lab=gsub("Age group: ","",parameter2_lab),
                  parameter2_lab=gsub("Education: ","",parameter2_lab),
    )
  corr_cova_labels$parameter2_lab = factor(corr_cova_labels$parameter2_lab,
                                           levels=corr_cova_labels$parameter2_lab)
  if(cov_=="ALL") {
    cov_ = corr_cova_labels$parameter2
  }
  cov_label = corr_cova_labels %>% 
    dplyr::filter(parameter2==cov_) %>% 
    dplyr::pull(parameter2_lab) %>% 
    as.character()
  
  # manual selection for age (as this string appears in many column names)
  if(cov_[[1]]=="age") {
    coefs_ = coefs_ %>% 
      dplyr::filter(parameter %in% c("age","hu_age")) 
  }
  
  # select coefs
  adj_levels = unique(coefs_$adjusted)
  coefs_2 = coefs_ %>% 
    dplyr::filter(grepl(cov_,parameter)) %>% 
    ## order type
    dplyr::mutate(adjusted=factor(adjusted,levels=adj_levels)) %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<1 | exp_beta_lob>1,"*","")) %>% 
    ## text version of results
    dplyr::mutate(exp_beta_lab=qsumaway(exp_beta,exp_beta_lob,exp_beta_upb,away,digits=2)) %>% 
    dplyr::mutate(exp_beta_lab=paste0(away,exp_beta_lab,away)) 
  
  # manual adaptations for age
  if(cov_[[1]]=="age") {
    coefs_2 = coefs_2 %>% 
      dplyr::mutate(parameter2_lab="Per 10-year increase") %>% 
      dplyr::mutate(exp_beta=exp_beta^10,
                    exp_beta_lob=exp_beta_lob^10,
                    exp_beta_upb=exp_beta_upb^10) %>% 
      dplyr::mutate(exp_beta_lab=qsumaway(exp_beta,exp_beta_lob,exp_beta_upb,away,digits=2)) %>% 
      dplyr::mutate(exp_beta_lab=paste0(away,exp_beta_lab,away)) 
  }
  
  # format coefs
  coefs_nonzero = coefs_2 %>% 
    dplyr::filter(type=="Odds ratio of detection") %>% 
    ## wide format
    dplyr::select(dependent_variable,dependent_variable_lab,parameter2_lab,adjusted,exp_beta_lab) %>% 
    tidyr::pivot_wider(names_from=adjusted,values_from=exp_beta_lab) 
  names(coefs_nonzero)[names(coefs_nonzero)=="Unadjusted"] = "OR"
  names(coefs_nonzero)[names(coefs_nonzero)=="Adjusted, simple"] = "aOR1"
  names(coefs_nonzero)[names(coefs_nonzero)=="Adjusted, extended"] = "aOR2"
  names(coefs_nonzero)[names(coefs_nonzero)=="Adjusted, BART selection"] = "aOR3"
  names(coefs_nonzero)[names(coefs_nonzero)=="Adjusted, penalization"] = "aOR4"
  if(!any(grepl("aOR3",names(coefs_nonzero)))) coefs_nonzero = coefs_nonzero %>% dplyr::mutate(aOR3="--",.after = aOR2)
  coefs_quant = coefs_2 %>% 
    dplyr::filter(type=="Fold-change") %>% 
    ## wide format
    dplyr::select(dependent_variable,dependent_variable_lab,parameter2_lab,adjusted,exp_beta_lab) %>% 
    tidyr::pivot_wider(names_from=adjusted,values_from=exp_beta_lab)
  names(coefs_quant)[names(coefs_quant)=="Unadjusted"] = "FC"
  names(coefs_quant)[names(coefs_quant)=="Adjusted, simple"] = "aFC1"
  names(coefs_quant)[names(coefs_quant)=="Adjusted, extended"] = "aFC2"
  names(coefs_quant)[names(coefs_quant)=="Adjusted, BART selection"] = "aFC3"
  names(coefs_quant)[names(coefs_quant)=="Adjusted, penalization"] = "aFC4"
  if(!any(grepl("aFC3",names(coefs_quant)))) coefs_quant = coefs_quant %>% dplyr::mutate(aFC3="--",.after = aFC2)
  # remove aORs and aFCs corresponding to adjustment variables (age and sex for simple adjustment, age, sex and others for extended adjustment)
  if(!aOR1) {
    coefs_nonzero = coefs_nonzero %>%
      dplyr::mutate(aOR1="--") %>%
      dplyr::relocate(aOR1,.after = OR)
    coefs_quant = coefs_quant %>%
      dplyr::mutate(aFC1="--") %>%
      dplyr::relocate(aFC1,.after = FC)
  }
  if(!aOR2) {
    coefs_nonzero = coefs_nonzero %>%
      dplyr::mutate(aOR2="--") %>%
      dplyr::relocate(aOR2,.after = aOR1)
    coefs_quant = coefs_quant %>%
      dplyr::mutate(aFC2="--") %>%
      dplyr::relocate(aFC2,.after = aFC1)
  }
  # replace NA by -- for variables not selected by BART 
  coefs_nonzero =  coefs_nonzero %>%
    dplyr::mutate(aOR3=ifelse(is.na(aOR3),"--",aOR3))
  coefs_quant = coefs_quant %>%
    dplyr::mutate(aFC3=ifelse(is.na(aFC3),"--",aFC3))
  
  # format desc
  tab_1 = coefs_nonzero %>% 
    dplyr::select(-dependent_variable) 
  ## manual relabel FFQ
  if(grepl("_freq",cov_)) {
    tab_1$parameter2_lab = "Per doubling frequency"
  }
  
  # keep only BART-selected OR penalized 95%CrI away from 1
  if(summarise_results) {
      tab_1 = tab_1  %>% 
        dplyr::mutate(OR_bart=aOR3!="--") %>% 
        dplyr::mutate(OR_pen=grepl("[**]",aOR4)) %>% 
        dplyr::group_by(dependent_variable_lab) %>% 
        dplyr::mutate(keep=max(OR_bart,OR_pen)) %>% 
        dplyr::filter(keep==1) %>% 
        dplyr::select(-OR_bart,-OR_pen,-keep) %>% 
        dplyr::ungroup()
  }
  
  tab_1 = tab_1 %>% 
    dplyr::arrange(dependent_variable_lab,parameter2_lab) %>% 
    dplyr::mutate(dependent_variable_lab=if_else(dependent_variable_lab!=lag(dependent_variable_lab,1),
                                                 dependent_variable_lab ,
                                                 "",
                                                 dependent_variable_lab)) %>% 
    # tidyr::replace_na(list(OR="--",aOR="--",apOR="--")) %>% 
    gt::gt() %>% 
    gt::tab_style(style = cell_text(weight = "bold"),locations = cells_column_labels()) %>% 
    gt::fmt_markdown(column=c(OR,aOR1,aOR2,aOR3,aOR4)) %>% 
    gt::cols_label(dependent_variable_lab="PFAS",
                   parameter2_lab=cov_label) %>% 
    gt::tab_caption(md(paste0("**Table.** Association between PFAS detection and ", tolower(cov_label),"."))) %>% 
    gt::cols_align("center") %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate (95% CrI)"),
                     locations = gt::cells_column_labels(columns = OR)) %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate adjusted on age and gender (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aOR1)) %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate adjusted on age, gender, nationality, education, income and center (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aOR2)) %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate adjusted on the top 10 covariates selected by BART (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aOR3)) %>% 
    gt::tab_footnote(footnote = md("Penalized odds ratio estimate adjusted on all other variables (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aOR4))
  
  tab_2 = coefs_quant %>% 
    dplyr::select(-dependent_variable)
  ## manual relabel FFQ
  if(grepl("_freq",cov_)) {
    tab_2$parameter2_lab = "Per doubling frequency"
  }
  
  # keep only lines with paFC CrI not including 1 or BART selected
  if(summarise_results) {
      tab_2 = tab_2  %>% 
        dplyr::mutate(FC_bart=aFC3!="--") %>% 
        dplyr::mutate(FC_pen=grepl("[**]",aFC4)) %>% 
        dplyr::group_by(dependent_variable_lab) %>% 
        dplyr::mutate(keep=max(FC_bart,FC_pen)) %>% 
        dplyr::filter(keep==1) %>% 
        dplyr::select(-FC_bart,-FC_pen,-keep) %>% 
        dplyr::ungroup()
  }
  tab_2 = tab_2  %>% 
    dplyr::mutate(dependent_variable_lab=if_else(dependent_variable_lab!=lag(dependent_variable_lab,1),
                                                 dependent_variable_lab ,
                                                 "",
                                                 dependent_variable_lab)) %>% 
    gt::gt() %>% 
    gt::tab_style(style = cell_text(weight = "bold"),locations = cells_column_labels()) %>% 
    gt::fmt_markdown(column=c(FC,aFC1,aFC2,aFC3,aFC4)) %>% 
    gt::cols_label(dependent_variable_lab="PFAS",
                   parameter2_lab=cov_label) %>% 
    gt::tab_caption(md(paste0("**Table.** Association between PFAS quantity (among detected) and ", tolower(cov_label),"."))) %>% 
    gt::cols_align("center") %>% 
    gt::tab_footnote(footnote = md("Fold-change estimate (95% CrI)"),
                     locations = gt::cells_column_labels(columns = FC)) %>% 
    gt::tab_footnote(footnote = md("Fold-change estimate adjusted on age and gender (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aFC1)) %>% 
    gt::tab_footnote(footnote = md("Fold-change estimate adjusted on age, gender, nationality, education, income and center (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aFC2)) %>% 
    gt::tab_footnote(footnote = md("Fold-change estimate adjusted on the top 10 covariates selected by BART (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aFC3)) %>% 
    gt::tab_footnote(footnote = md("Penalized fold-change estimate adjusted on all other variables (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aFC4))
  
  return(list(tab_1,tab_2))
}


fn_302_table_regression_by_covariate_cluster = function(cov_,coefs_,age_per10=TRUE, ...) {
  # labels
  labs2 <- list(age="Age (per 10-year)",
                genderMales='Sex: male',
                a1_foreignForeignnationality='Nationality: non-Swiss',
                a1_education_3classesApprenticeshipDprofessionalbaccalaureate = "Education: apprenticeship",
                a1_education_3classesUniversity ="Education: university",
                `a1_household_monthly_income_3classesCHF4500to9000`="Income: intermediate", 
                `a1_household_monthly_income_3classes>CHF9000`="Income: high", 
                center_labelLausanne="Center: Lausanne",
                
                cluster2_dairy_focused="Cluster: dairy-focused",
                cluster3_high_fiber="Cluster: mediterranean",
                cluster4_plant_based="Cluster: plant-based",
                cluster5_meat_centered="Cluster: western",
                
                ex1_fish_freq="Overall fish consumption",
                ex1_tunaFrequentconsumptionoftuna='Tuna',
                ex1_salmonFrequentconsumptionofsalmon='Salmon',
                ex1_pangasiusFrequentconsumptionofpangasius='Pangasius',
                ex1_troutFrequentconsumptionoftrout='Trout',
                ex1_shrimpFrequentconsumptionofshrimp='Shrimp',
                ex1_swiss_fishConsumptionoffishfromSwisslakes="Local fish from swiss lakes",
                
                ex1_veggie_veganVegetarianDvegandiet='Vegetarian or vegan',
                ex1_tapwaterMorethan2Loftapwaterperday='Tap water (>2L per day)',
                ex1_alcool_weeklyAtleastweekly='Alcohol: weekly',
                ex1_alcool_weeklyDailyandmore='Alcohol: daily',
                ex1_cigPastsmoker='Tobacco: past user',
                ex1_cigCurrentsmoker='Tobacco: current user',
                
                ex2_skiwaxExpositiontoskiwax='Ski wax',
                ex2_spray_impExpositiontoimpregnationsprays='Impregnation sprays',
                ex1_hot_meals_package_freqAtleastweekly="Disposable food packaging",
                ex1_residence_urbanResidenceinurbanDindustrialarea='Residence in urban/industrial area ',
                ex1_carpetSyntheticornaturalcarpet="Synthetic or natural carpet at home",
                
                ex1_deo_spray_freqAtleastweekly="Deodorant spray",
                ex1_deo_roll_freqAtleastweekly="Deodorant roll",
                ex1_bodylotion_freqAtleastweekly="Body lotion",
                ex1_handcream_freqAtleastweekly="Hand cream",
                ex1_perfume_freqAtleastweekly="Perfume",
                ex1_hairspray_freqAtleastweekly="Hair spray",
                
                ex2_prof_smokeProfessionalexpositiontosmoke='Smoke at work',
                ex2_prof_exhaust_gasProfessionalexpositiontoexhaustgas='Exhaust gas at work',
                ex2_prof_solvent_vapoursProfessionalexpositiontosolventvapours='Solvent vapours at work',
                ex2_prof_cleaning_agentsProfessionalexpositiontocleaningagents='Cleaning agents at work',
                ex2_prof_dustProfessionalexpositiontodust='Dust at work',
                ex2_prof_skiwaxProfessionalexpositiontoskiwax='Ski wax at work',
                ex2_prof_spray_impProfessionalexpositiontoimpregnationsprays='Impregnation sprays at work'
  )
  corr_pfas_labels = tibble(dependent_variable=c("mu2","mu3"),
                            dependent_variable_lab=controls$lab_pfas_cluster[2:3])
  corr_cova_labels = tibble(parameter=names(labs2),
                            parameter_lab=unname(unlist(labs2))) 
  
  cov_label = corr_cova_labels %>% 
    dplyr::filter(parameter==cov_) %>% 
    dplyr::pull(parameter_lab) %>% 
    as.character()
  
  if(age_per10) {
    coefs_ = coefs_ %>% 
      dplyr::mutate(exp_beta=ifelse(parameter=="age",exp(10*Estimate),exp_beta),
                    exp_beta_lob =ifelse(parameter=="age",exp(10*`l-95% CI`),exp_beta_lob),
                    exp_beta_upb =ifelse(parameter=="age",exp(10*`u-95% CI`),exp_beta_upb),
                    exp_beta_upb =ifelse(parameter=="age",exp(10*`u-95% CI`),exp_beta_upb))
  }
  
  # select coefs
  adj_levels = unique(coefs_$adjusted)
  coefs_2 = coefs_ %>% 
    ## order type
    dplyr::mutate(adjusted=factor(adjusted,levels=adj_levels)) %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<1 | exp_beta_lob>1,"*","")) %>% 
    ## text version of results
    dplyr::mutate(exp_beta_lab=qsumaway(exp_beta,exp_beta_lob,exp_beta_upb,away,digits=2)) %>% 
    dplyr::mutate(exp_beta_lab=paste0(away,exp_beta_lab,away)) 
  
 
  # format coefs
  coefs_nonzero = coefs_2 %>% 
    ## wide format
    dplyr::select(dependent_variable_lab,parameter_lab,exp_beta_lab) %>% 
    tidyr::pivot_wider(names_from=dependent_variable_lab,values_from=exp_beta_lab) 
  names(coefs_nonzero)[1] = "Variable"
  names(coefs_nonzero)[2:3] = paste("OR for",names(coefs_nonzero)[2:3])
  
  # format desc
  tab_1 = coefs_nonzero %>%
    # gtsummary::modify_caption(") %>% 
    gt::gt() %>% 
    gt::tab_style(style = cell_text(weight = "bold"),locations = cells_column_labels()) %>% 
    gt::fmt_markdown() %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate (95% CrI)"),
                     locations = gt::cells_column_labels(columns = 2:3))  %>% 
    gt::tab_caption(gt::md("**Table**. Results from multinomial regression models showing the association between individual variables and three PFAS clusters, expressed as odds ratios (ORs) relative to cluster #1 (low exposure)."))
  
  
  return(tab_1)
}


fn_302_table_regression_by_covariate2 = function(cov_,coefs_,summarise_results=FALSE, cov_dataname=NULL, aOR1=TRUE, aOR2=TRUE, ...) {
  # labels
  labs2 <- list(age="Age",
                genderMales='Gender: male',
                a1_foreignForeignnationality='Nationality: non-Swiss',
                a1_education_3classesApprenticeshipDprofessionalbaccalaureate = "Education: apprenticeship",
                a1_education_3classesUniversity ="Education: university",
                `a1_household_monthly_income_3classesCHF4500to9000`="Income: intermediate", 
                `a1_household_monthly_income_3classes>CHF9000`="Income: high", 
                center_labelLausanne="Center: Lausanne",
                
                cluster2_dairy_focused="Cluster: dairy focused",
                cluster3_high_fiber="Cluster: high-fiber",
                cluster4_plant_based="Cluster: plant-based",
                cluster5_meat_centered="Cluster: meat-centered",
                
                ex1_fish_freq="Overall fish consumption",
                ex1_tunaFrequentconsumptionoftuna='Tuna',
                ex1_salmonFrequentconsumptionofsalmon='Salmon',
                ex1_pangasiusFrequentconsumptionofpangasius='Pangasius',
                ex1_troutFrequentconsumptionoftrout='Trout',
                ex1_shrimpFrequentconsumptionofshrimp='Shrimp',
                ex1_swiss_fishConsumptionoffishfromSwisslakes="Local fish from swiss lakes",
                
                ex1_veggie_veganVegetarianDvegandiet='Vegetarian or vegan',
                ex1_tapwaterMorethan2Loftapwaterperday='Tap water (>2L per day)',
                ex1_alcool_weeklyAtleastweekly='Alcohol: weekly',
                ex1_alcool_weeklyDailyandmore='Alcohol: daily',
                ex1_cigPastsmoker='Tobacco: past user',
                ex1_cigCurrentsmoker='Tobacco: current user',
                
                ex2_skiwaxExpositiontoskiwax='Ski wax',
                ex2_spray_impExpositiontoimpregnationsprays='Impregnation sprays',
                ex1_hot_meals_package_freqAtleastweekly="Disposable food packaging",
                ex1_residence_urbanResidenceinurbanDindustrialarea='Residence in urban/industrial area ',
                ex1_carpetSyntheticornaturalcarpet="Synthetic or natural carpet at home",
                
                ex1_deo_spray_freqAtleastweekly="Deodorant spray",
                ex1_deo_roll_freqAtleastweekly="Deodorant roll",
                ex1_bodylotion_freqAtleastweekly="Body lotion",
                ex1_handcream_freqAtleastweekly="Hand cream",
                ex1_perfume_freqAtleastweekly="Perfume",
                ex1_hairspray_freqAtleastweekly="Hair spray",
                
                ex2_prof_smokeProfessionalexpositiontosmoke='Smoke at work',
                ex2_prof_exhaust_gasProfessionalexpositiontoexhaustgas='Exhaust gas at work',
                ex2_prof_solvent_vapoursProfessionalexpositiontosolventvapours='Solvent vapours at work',
                ex2_prof_cleaning_agentsProfessionalexpositiontocleaningagents='Cleaning agents at work',
                ex2_prof_dustProfessionalexpositiontodust='Dust at work',
                ex2_prof_skiwaxProfessionalexpositiontoskiwax='Ski wax at work',
                ex2_prof_spray_impProfessionalexpositiontoimpregnationsprays='Impregnation sprays at work'
  )
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(labs2),
                            parameter2_lab=unname(unlist(labs2))) 
  
  
  if(cov_=="ALL") {
    cov_ = corr_cova_labels$parameter
  }
  cov_label = corr_cova_labels %>% 
    dplyr::filter(parameter2==cov_) %>% 
    dplyr::pull(parameter2_lab) %>% 
    as.character()
  
  # manual selection for age (as this string appears in many column names)
  if(cov_[[1]]=="age") {
    coefs_ = coefs_ %>% 
      dplyr::filter(parameter %in% c("age","hu_age")) 
  }
  
  # select coefs
  adj_levels = unique(coefs_$adjusted)
  coefs_2 = coefs_ %>% 
    dplyr::filter(grepl(cov_,parameter)) %>% 
    ## order type
    dplyr::mutate(adjusted=factor(adjusted,levels=adj_levels)) %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<1 | exp_beta_lob>1,"*","")) %>% 
    ## text version of results
    dplyr::mutate(exp_beta_lab=qsumaway(exp_beta,exp_beta_lob,exp_beta_upb,away,digits=2)) %>% 
    dplyr::mutate(exp_beta_lab=paste0(away,exp_beta_lab,away)) 
  
  # manual adaptations for age
  if(cov_[[1]]=="age") {
    coefs_2 = coefs_2 %>% 
      dplyr::mutate(parameter2_lab="Per 10-year increase") %>% 
      dplyr::mutate(exp_beta=exp_beta^10,
                    exp_beta_lob=exp_beta_lob^10,
                    exp_beta_upb=exp_beta_upb^10) %>% 
      dplyr::mutate(exp_beta_lab=qsumaway(exp_beta,exp_beta_lob,exp_beta_upb,away,digits=2)) %>% 
      dplyr::mutate(exp_beta_lab=paste0(away,exp_beta_lab,away)) 
  }
  
  # format coefs
  coefs_nonzero = coefs_2 %>% 
    dplyr::filter(type=="Odds ratio of detection") %>% 
    ## wide format
    dplyr::select(dependent_variable,dependent_variable_lab,parameter2_lab,adjusted,exp_beta_lab) %>% 
    tidyr::pivot_wider(names_from=adjusted,values_from=exp_beta_lab) 
  names(coefs_nonzero)[names(coefs_nonzero)=="Unadjusted"] = "OR"
  names(coefs_nonzero)[names(coefs_nonzero)=="Adjusted, simple"] = "aOR1"
  names(coefs_nonzero)[names(coefs_nonzero)=="Adjusted, extended"] = "aOR2"

  coefs_quant = coefs_2 %>% 
    dplyr::filter(type=="Fold-change") %>% 
    ## wide format
    dplyr::select(dependent_variable,dependent_variable_lab,parameter2_lab,adjusted,exp_beta_lab) %>% 
    tidyr::pivot_wider(names_from=adjusted,values_from=exp_beta_lab)
  names(coefs_quant)[names(coefs_quant)=="Unadjusted"] = "FC"
  names(coefs_quant)[names(coefs_quant)=="Adjusted, simple"] = "aFC1"
  names(coefs_quant)[names(coefs_quant)=="Adjusted, extended"] = "aFC2"

  # format desc
  tab_1 = coefs_nonzero %>% 
    dplyr::select(-dependent_variable)  %>% 
    dplyr::arrange(dependent_variable_lab,parameter2_lab) %>% 
    dplyr::mutate(dependent_variable_lab=if_else(dependent_variable_lab!=lag(dependent_variable_lab,1),
                                                 dependent_variable_lab ,
                                                 "",
                                                 dependent_variable_lab)) %>% 
    # tidyr::replace_na(list(OR="--",aOR="--",apOR="--")) %>% 
    gt::gt() %>% 
    gt::tab_style(style = cell_text(weight = "bold"),locations = cells_column_labels()) %>% 
    gt::fmt_markdown(column=c(OR,aOR1,aOR2)) %>% 
    gt::cols_label(dependent_variable_lab="PFAS",
                   parameter2_lab="Variable") %>% 
    gt::tab_caption(md(paste0("**Table.** Association with PFAS detection."))) %>% 
    gt::cols_align("center") %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate (95% CrI)"),
                     locations = gt::cells_column_labels(columns = OR)) %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate adjusted on age and gender (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aOR1)) %>% 
    gt::tab_footnote(footnote = md("Odds ratio estimate adjusted on age, gender, nationality, education, income and center (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aOR2)) 
  
  tab_2 = coefs_quant %>% 
    dplyr::select(-dependent_variable) %>% 
    dplyr::mutate(dependent_variable_lab=if_else(dependent_variable_lab!=lag(dependent_variable_lab,1),
                                                 dependent_variable_lab ,
                                                 "",
                                                 dependent_variable_lab)) %>% 
    gt::gt() %>% 
    gt::tab_style(style = cell_text(weight = "bold"),locations = cells_column_labels()) %>% 
    gt::fmt_markdown(column=c(FC,aFC1,aFC2)) %>% 
    gt::cols_label(dependent_variable_lab="PFAS",
                   parameter2_lab="Variable") %>% 
    gt::tab_caption(md(paste0("**Table.** Association with PFAS quantity upon detection."))) %>% 
    gt::cols_align("center") %>% 
    gt::tab_footnote(footnote = md("Fold-change estimate (95% CrI)"),
                     locations = gt::cells_column_labels(columns = FC)) %>% 
    gt::tab_footnote(footnote = md("Fold-change estimate adjusted on age and gender (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aFC1)) %>% 
    gt::tab_footnote(footnote = md("Fold-change estimate adjusted on age, gender, nationality, education, income and center (95% CrI)"),
                     locations = gt::cells_column_labels(columns = aFC2)) 
  
  return(list(tab_1,tab_2))
}
