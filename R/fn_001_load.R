#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 001 - Formatting and selection of relevant variables
#' author: jriou
#' date: 2024-06-12
#' ---

fn_001_load = function(...) {
  
  # Demographics ------------------------------------------------------------
  
  demo = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.1.1_Demographics_v2.0_20230406_cleaned.csv"),
                                      col_select = -1)) %>% 
    dplyr::mutate(gender=if_else(gender==2,"Females","Males"),
                  gender=factor(gender,levels=c("Males","Females")),
                  center_label=if_else(center==0,"Bern","Lausanne"),
                  center_label=factor(center_label,levels=c("Bern","Lausanne"))
    )
  
  # Screening ---------------------------------------------------------------
  
  scre = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.1.2_Screening_v3.0_20240312_cleaned.csv"),
                                      col_select = -1)) %>% 
    # relevant for PFAS/immunity: smoking, food
    dplyr::mutate(did=did,
                  date=lubridate::as_date(lubridate::dmy_hms(start_visit)),
                  # scre_smoker=if_else(pef_sq_scr3_c_tobacco==3,1,0),
                  # scre_smoked_24h=if_else(pef_sq_scr3_d_24hsmoke==2,1,0),
                  # scre_food_nuts=pef_sq_opea_q1_a,
                  # scre_food_fresh_fruit=pef_sq_opea_q1_b,
                  # scre_food_fresh_vegetables=pef_sq_opea_q1_c,
                  # scre_food_fresh_spices=pef_sq_opea_q1_d,
                  # scre_food_fresh_crude_oil=pef_sq_opea_q1_e,
                  # scre_food_fresh_processed=pef_sq_opea_q1_f,
                  # scre_food_fresh_deep_fried=pef_sq_opea_q1_g,
                  # scre_food_fresh_coffee_tea=pef_sq_opea_q1_h,
                  # scre_food_fresh_red_wine=pef_sq_opea_q1_i,
                  .keep="none"
    )
  
  # A1 ----------------------------------------------------------------------
  
  a1 = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.2.1_A1_v2.0_20230406_cleaned.csv"),
                                      col_select = -1)) %>% 
    # relevant for PFAS/immunity: foreign nationality, university education, extreme incomes
    dplyr::mutate(did=did,
                  a1_foreign=if_else(a1_q3_nat==3,"Foreign nationality","Swiss nationality"),
                  a1_foreign=factor(a1_foreign,levels=c("Swiss nationality","Foreign nationality")),
                  a1_education=case_when(a1_q7_training==1 ~ "Primary school",
                                         a1_q7_training==2 ~ "Secondary school",
                                         a1_q7_training==3 ~ "High school",
                                         a1_q7_training==4 ~ "Apprenticeship / professional baccalaureate",
                                         a1_q7_training==5 ~ "University: Bachelors",
                                         a1_q7_training==6 ~ "University: Masters/Licence",
                                         a1_q7_training==7 ~ "University: Doctorate/PhD",
                                         a1_q7_training==8 ~ "Other",
                                         a1_q7_training==9 ~ NA_character_,
                                         TRUE ~ NA_character_,
                                         .ptype = factor(levels = c("Primary school","Secondary school","High school","Apprenticeship / professional baccalaureate","University: Bachelors","University: Masters/Licence","University: Doctorate/PhD","Other"),ordered=FALSE)),
                  a1_household_monthly_income=case_when(a1_q8_salary == 1 ~ "<CHF 3,000",
                                                        a1_q8_salary == 2 ~ "CHF 3,000 to 4,500",
                                                        a1_q8_salary == 3 ~ "CHF 4,500 to 6,000",
                                                        a1_q8_salary == 4 ~ "CHF 6,000 to 9,000",
                                                        a1_q8_salary == 5 ~ "CHF 9,000 to 11,000",
                                                        a1_q8_salary == 6 ~ ">CHF 11,000",
                                                        a1_q8_salary == 7 ~ NA_character_,
                                                        TRUE ~ NA_character_,
                                                        .ptype = factor(levels = c("<CHF 3,000","CHF 3,000 to 4,500","CHF 4,500 to 6,000","CHF 6,000 to 9,000","CHF 9,000 to 11,000",">CHF 11,000"),ordered=FALSE)),
                  .keep="none"
    ) 
  
  # EX1 ---------------------------------------------------------------------
  
  ex1 = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.2.2_EX1_v2.0_20230406_cleaned.csv"),
                                      col_select = -1)) %>% 
    # relevant for PFAS/immunity: household size, hot meals in package, veggie diet, specific foods (1 never to 8 four times a day)
    dplyr::mutate(did=did,
                  ex1_hot_meals_package=ex1_q40_mcdo,
                  ex1_q41_diet=ex1_q41_diet,
                  across(starts_with("ex1_q34_q1_q"),as.numeric),
                  across(starts_with("ex1_q34_q2_q"),as.numeric),
                  across(starts_with("ex1_q28_q"),as.factor),
                  across(starts_with("ex1_q29_q"),as.factor),
                  ex1_tuna=ex1_q35_fish___1,
                  ex1_salmon=ex1_q35_fish___2,
                  ex1_pangasius=ex1_q35_fish___3,
                  # ex1_tilapia=ex1_q35_fish___4,
                  ex1_trout=ex1_q35_fish___5,
                  ex1_shrimp=ex1_q35_fish___6,
                  # ex1_fishfingers=ex1_q35_fish___7,   
                  across(starts_with("ex1_q25_q"),as.numeric),
                  ex1_residence_city=ex1_q3_main_home___1,
                  ex1_residence_suburb=ex1_q3_main_home___2,
                  ex1_residence_industrial=ex1_q3_main_home___3,
                  ex1_residence_countryside=ex1_q3_main_home___4,
                  ex1_residence_nearfield=ex1_q3_main_home___5,
                  ex1_tapwater=ex1_q44_q1_water,
                  ex1_q10_q1_lounge=ex1_q10_q1_lounge,
                  ex1_q10_q2_bedroom=ex1_q10_q2_bedroom,
                  ex1_q10_q3_rest=ex1_q10_q3_rest,
                  ex1_distance_field=ex1_q3_q1_main_home,
                  ex1_alcool=ex1_q44_q9_alcool,
                  .keep="none") %>% 
    dplyr::rename(ex1_q34_q2_q11_cornflakes=ex1_q34_q2_q11_bio,
                  ex1_q34_q2_q13_patat=ex1_q34_q2_q13_bio) %>% 
    dplyr::rename_with(~ gsub("q34_q1_q[[:digit:]]*_?","",.x), starts_with("ex1_q34_q1_q")) %>% 
    dplyr::rename_with(~ gsub("q34_q2_q[[:digit:]]*_?","organic_",.x), starts_with("ex1_q34_q2_q")) %>% 
    dplyr::rename_with(~ gsub("q25_q[[:digit:]]*_?","",.x), starts_with("ex1_q25_q")) %>% 
    dplyr::rename(ex1_handcream=ex1_hand,
                  ex1_bodylotion=ex1_body,
                  ex1_hairspray=ex1_hair,
                  ex1_toothpaste=ex1_tooth)
  
  
  # CORONA ------------------------------------------------------------------
  
  corona = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.2.5_Corona_v2.0_20230406_cleaned.csv"),
                                      col_select = -1)) %>% 
    # relevant for PFAS/immunity: past disease, severity, vaccines
    dplyr::mutate(did=did,
                  corona_past_disease_confirmed=if_else(cov_q1_disease==1,1,0),
                  corona_past_disease_suspected=if_else(cov_q1_disease %in% 1:2,1,0),
                  corona_past_disease_date=lubridate::dmy(cov_q1_disease_when),
                  corona_past_positive_test=as.numeric(if_any(starts_with("cov_q3_result_test"), ~ .==1)),
                  corona_past_hospit=cov_q5_treatment___4,
                  corona_past_icu=if_else(cov_q5_treatement_hosp_b==0 | cov_q5_treatement_hosp_c==0,1,0,0), #TODO check that 0 actually corresponds to "yes"
                  corona_vaccine_sarscov2=if_else(cov_q7_vaccinations_q==0,1,0), #TODO check that 0 actually corresponds to "yes"
                  .keep="none") %>% 
    # check consistency and clear out NAs
    dplyr::mutate(corona_past_disease_confirmed=if_else(corona_past_disease_confirmed==1 | corona_past_positive_test==1 | corona_past_hospit==1 | corona_past_icu==1,1,0,0),
                  corona_past_disease_suspected=if_else(corona_past_disease_suspected==1 | corona_past_disease_confirmed==1 | corona_past_positive_test==1 | corona_past_hospit==1 | corona_past_icu==1,1,0,0),
    )
  
  
  # EX2 ---------------------------------------------------------------------
  
  ex2 = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.2.7_EX2_v2.0_20230406_cleaned.csv"),
                                      col_select = -1)) %>% 
    # relevant for PFAS/immunity: exposure to herbicides, ski fart, or impregnation sprays, profession (https://www.kubb-tool.bfs.admin.ch/fr), professional exposures
    dplyr::mutate(did=did,
                  ex2_garden_pesticides=case_when(ex2_q9_garden==1 ~ 1, ex2_q9_garden==2 ~ 0, ex2_q9_garden==3 ~ NA_integer_),
                  ex2_skiwax=case_when(ex2_q10_ski==1 ~ 1, ex2_q10_ski==2 ~ 0, ex2_q10_ski==3 ~ NA_integer_),
                  ex2_spray_imp=case_when(ex2_q11_spray_imp==1 ~ 1, ex2_q11_spray_imp==2 ~ 0, ex2_q11_spray_imp==3 ~ NA_integer_),
                  # ex2_q2_q1_job=ex2_q2_q1_job,
                  # ex2_q2_q4_job=ex2_q2_q4_job,
                  # ex2_q2_q4_job_sector1=substr(ex2_q2_q4_job,1,1),
                  # ex2_q2_q4_job_sector2=substr(ex2_q2_q4_job,1,2),
                  ex2_prof_smoke=case_when(ex2_q6_q1_smoke==1 ~ 1, ex2_q6_q1_smoke==2 ~ 0, ex2_q6_q1_smoke==3 ~ NA_integer_),
                  ex2_prof_exhaust_gas=case_when(ex2_q6_q2_gas==1 ~ 1, ex2_q6_q2_gas==2 ~ 0, ex2_q6_q2_gas==3 ~ NA_integer_),
                  ex2_prof_solvent_vapours=case_when(ex2_q6_q3_solvent==1 ~ 1, ex2_q6_q3_solvent==2 ~ 0, ex2_q6_q3_solvent==3 ~ NA_integer_),
                  ex2_prof_cleaning_agents=case_when(ex2_q6_q4_clean==1 ~ 1, ex2_q6_q4_clean==2 ~ 0, ex2_q6_q4_clean==3 ~ NA_integer_),
                  ex2_prof_dust=case_when(ex2_q6_q5_dust==1 ~ 1, ex2_q6_q5_dust==2 ~ 0, ex2_q6_q5_dust==3 ~ NA_integer_),
                  ex2_prof_polycarbonate=case_when(ex2_q7_q1_polycarb==1 ~ 1, ex2_q7_q1_polycarb==2 ~ 0, ex2_q7_q1_polycarb==3 ~ NA_integer_),
                  ex2_prof_epoxy_resin=case_when(ex2_q7_q2_resin==1 ~ 1, ex2_q7_q2_resin==2 ~ 0, ex2_q7_q2_resin==3 ~ NA_integer_),
                  ex2_prof_skiwax=case_when(ex2_q7_q4_wax==1 ~ 1, ex2_q7_q4_wax==2 ~ 0, ex2_q7_q4_wax==3 ~ NA_integer_),
                  ex2_prof_spray_imp=case_when(ex2_q7_q5_spray==1 ~ 1, ex2_q7_q5_spray==2 ~ 0, ex2_q7_q5_spray==3 ~ NA_integer_),
                  .keep="none")
  
  
  
  # BLOODURINE --------------------------------------------------------------
  
  bu = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.4.1_BloodUrine_param_v2.1_20240320.csv"),
                                      col_select = -1)) %>% 
    # relevant for PFAS/immunity: SARS-CoV-2-related measurements
    dplyr::mutate(did=did,
                  bu_antiS_igg_nonzero=if_else(labo_cov2g=="#p",1,0),
                  bu_antiS_igg_quant=if_else(labo_cov2r=="#non",NA_real_,suppressWarnings(as.numeric(labo_cov2r)),NA_real_), # TODO: not sure about the unit
                  bu_antiN_igg_nonzero=if_else(labo_cov2ng=="#p",1,0),
                  bu_antiN_igg_quant=as.numeric(labo_cov2nr),
                  bu_ratio_iga_nonzero=if_else(labo_cov2a=="#p",1,0),
                  bu_ratio_iga_quant=as.numeric(labo_cov2ar),
                  bu_ige_pos=as.numeric(labo_cov2e), 
                  .keep="none") 
  
  
  # PFAS --------------------------------------------------------------------
  
  pfas_raw = suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.4.5_PFAS_v1.0_20230901.csv"),
                                               col_select = -1))
  pfas = pfas_raw %>% 
    dplyr::select(did,pfas_not_done)
  for(i in 1:length(controls$pfas_list)) {
    tmp = pull(pfas_raw[,controls$pfas_list[i]])
    out = suppressWarnings(tibble(nonzero = if_else(tmp=="<0.1",0,1),
                                  quant = case_when(tmp=="<0.1" ~ NA_real_,
                                                    tmp==">30" ~ 30,
                                                    TRUE ~ as.numeric(tmp))))
    names(out) = paste(controls$pfas_list[i],names(out),sep="_")
    pfas = dplyr::bind_cols(pfas,out)
    rm(tmp,out)
  }
  
  # HM ---------------------------------------------------------------------
  
  hm = 
    suppressMessages(readr::read_csv2(file.path(controls$data_path,"1.2.3_HM_v2.0_20230406_cleaned.csv"),
                                      col_select = -1)) %>% 
    dplyr::mutate(did=did,
                  hm_kids=hm_q121_q1_kid,
                  hm_kids=tidyr::replace_na(hm_kids,0),
                  hm_breastfeeding=if_else(hm_q122_q1_q3_kid1f==1 | hm_q122_q1_q3_kid2f==1 | hm_q122_q1_q3_kid3==1 | hm_q122_q1_q3_kid4==1,1,0,0),                  
                  hm_breastfeeding=factor(hm_breastfeeding,levels=0:1,labels=c("No breastfeeding","Breastfeeding")),
                  hm_q122_q1_q4_kid1f=tidyr::replace_na(hm_q122_q1_q4_kid1f,0),
                  hm_q122_q1_q4_kid2f=tidyr::replace_na(hm_q122_q1_q4_kid2f,0),
                  hm_q122_q1_q4_kid3=tidyr::replace_na(hm_q122_q1_q4_kid3,0),
                  hm_q122_q1_q4_kid4=tidyr::replace_na(hm_q122_q1_q4_kid4,0),
                  hm_breastfeeding_months=hm_q122_q1_q4_kid1f+hm_q122_q1_q4_kid2f+hm_q122_q1_q4_kid3+hm_q122_q1_q4_kid4,
                  hm_breastfeeding_months=tidyr::replace_na(hm_breastfeeding_months,0),
                  hm_breastfeeding_months=as.numeric(hm_breastfeeding_months),
                  .keep="none") %>% 
    dplyr::select(-hm_q122_q1_q4_kid1f,-hm_q122_q1_q4_kid2f,-hm_q122_q1_q4_kid3,-hm_q122_q1_q4_kid4)
  
  # ANT ---------------------------------------------------------------------
  
  ant = 
    suppressMessages(readr::read_delim(file.path(controls$data_path,"1.3.4_Anthropo_v2.0_20230411_cleaned.csv"),
                                       col_select=-1,delim=";")) %>% 
    dplyr::mutate(did=did,
                  ant_height=pef_anthro_4_height,
                  ant_weight= pef_anthro_4_weight,
                  ant_waist_circ=pef_anthro_4_waist_circ,
                  ant_bmi=pef_anthro_4_bmi,
                  .keep="none") 
  
  # Other data prepared by BAG (céline Fragnière) ---------------------------------------------------------------------
  
  od = 
    suppressMessages(readxl::read_xlsx(file.path(controls$data_path,"other_data_prepared_by_bag","SHeS pp - Publi PFAS I - data for stat analysis V2.xlsx"),
                                       sheet=1,
                                       skip=1,
                                       na = "-")) %>% 
    dplyr::select(did,ex1_swiss_fish=`35.1. Poissons des lacs suisses?`)
  
  # Join all ----------------------------------------------------------------
  
  dat_out = 
    demo %>% 
    dplyr::left_join(scre,by = join_by(did)) %>% 
    dplyr::left_join(a1,by = join_by(did)) %>% 
    dplyr::left_join(ex1,by = join_by(did)) %>% 
    dplyr::left_join(corona,by = join_by(did)) %>% 
    dplyr::left_join(ex2,by = join_by(did)) %>% 
    dplyr::left_join(bu,by = join_by(did)) %>% 
    dplyr::left_join(pfas,by = join_by(did)) %>% 
    dplyr::left_join(hm,by = join_by(did))  %>% 
    dplyr::left_join(ant,by = join_by(did)) %>% 
    dplyr::left_join(od,by = join_by(did))
  
  
  # Save --------------------------------------------------------------------
  
  write_rds(dat_out,file=file.path(controls$savepoint,"shesp_1.rds"))
  
  return(dat_out)
}