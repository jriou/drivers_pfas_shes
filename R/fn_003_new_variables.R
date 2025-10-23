#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 003 - New variables
#' author: jriou
#' date: 2024-06-12
#' ---

fn_003_new_variables = function(dat_in) {
  
  # Select dataset
  dat_0 = dat_in
  
  # Demographics: female as 0/1 and age groups
  
  dat_1 = dat_0 %>% 
    dplyr::mutate(female=if_else(gender=="Females","Female","Male"),
                  female=factor(female,levels=c("Male","Female")),
                  age_group=cut(age,c(20,30,40,50,60,70),include.lowest=TRUE,right=FALSE)) %>% 
    dplyr::relocate(female,age_group,.after=gender)
  
  # Screening: no new variable
  dat_2 = dat_1
  
  
  # A1: aggregate education and income in larger classes
  dat_3 = dat_2 %>% 
    dplyr::mutate(a1_education_3classes=case_when(a1_education %in% c("Primary school","Secondary school","High school","Other") ~ "Up to high school",
                                                  a1_education == "Apprenticeship / professional baccalaureate" ~ "Apprenticeship / professional baccalaureate",
                                                  a1_education %in% c("University: Bachelors","University: Masters/Licence","University: Doctorate/PhD") ~ "University",
                                                  a1_education == "Other" ~ NA_character_,
                                                  TRUE ~ NA_character_,
                                                  .ptype = factor(levels = c("Up to high school","Apprenticeship / professional baccalaureate","University"),ordered=FALSE)),
                  a1_household_monthly_income_3classes=case_when(a1_household_monthly_income %in% c("<CHF 3,000","CHF 3,000 to 4,500") ~ "<CHF 4,500",
                                                                 a1_household_monthly_income %in% c("CHF 4,500 to 6,000","CHF 6,000 to 9,000") ~ "CHF 4,500 to 9,000",
                                                                 a1_household_monthly_income %in% c("CHF 9,000 to 11,000",">CHF 11,000") ~ ">CHF 9,000",
                                                                 TRUE ~ NA_character_,
                                                                 .ptype = factor(levels = c("<CHF 4,500","CHF 4,500 to 9,000",">CHF 9,000"),ordered=FALSE))
    ) %>% 
    dplyr::relocate(a1_education_3classes,a1_household_monthly_income_3classes,.after=a1_household_monthly_income)
  
  
  
  # EX1: categorize consumption of organic products
  # composite groups: 
  # 1 "High consumption of non-organic foods"
  # 2 "Intermediate consumption of non-organic foods"
  # 3 "Low consumption of non-organic foods"
  dat_4 = dat_3
  ex1_organic_cat = c("ex1_organic_milk","ex1_organic_yog","ex1_organic_cheese","ex1_organic_egg",
                      "ex1_organic_poultry","ex1_organic_meat","ex1_organic_pork","ex1_organic_tofu",
                      "ex1_organic_fish","ex1_organic_starch","ex1_organic_cornflakes","ex1_organic_whole_grain",
                      "ex1_organic_patat","ex1_organic_lentil","ex1_organic_fruit","ex1_organic_vegetable",  
                      "ex1_organic_salad","ex1_organic_walnut","ex1_organic_sugar","ex1_organic_chips" )
  if(FALSE) { # REMOVE VARIABLES ABOUT ORGANIC PRODUCTS FOR NOW
    tmp = NULL
    for(tag in ex1_organic_cat) {
      dd = dplyr::bind_cols(dat_4[,tag],dat_4[,paste0(gsub("organic_","",tag))]) %>% 
        dplyr::rename(part_organic=1,consumption=2) %>% 
        dplyr::mutate(out=dplyr::case_when(consumption>=6 & part_organic %in% 1:2 ~ "Group 3",
                                           consumption>=6 & part_organic==3 ~ "Group 2",
                                           consumption>=6 & part_organic==4 ~ "Group 1",
                                           consumption%in% 4:5 & part_organic %in% 1:3 ~ "Group 2",
                                           consumption%in% 4:5 & part_organic==4 ~ "Group 1",
                                           consumption<=3 & part_organic %in% 1:2 ~ "Group 2",
                                           consumption<=3 & part_organic %in% 3:4 ~ "Group 1",
                                           TRUE ~ NA_character_,
                                           .ptype = factor(levels = c("Group 1","Group 2","Group 3"),ordered=FALSE))) %>% 
        dplyr::select(out)
      names(dd) = tag
      tmp = dplyr::bind_cols(tmp,dd)
    }
    
    dat_4 = dat_4 %>% 
      dplyr::select(-all_of(ex1_organic_cat)) %>% 
      dplyr::bind_cols(tmp) 
  } else {
    dat_4 = dat_4 %>% 
      dplyr::select(-all_of(ex1_organic_cat))
  }
  
  # EX1: how to treat questions based on frequency
  # 1, Rarely / Never | 2, Once a month | 3, Every 2 weeks | 4, 1-2 times a week | 5, 3-6 times a week | 6, Once a day | 7, 2-3 times a day | 8, 4 times a day and more
  ex1_dichotomize = c("ex1_milk","ex1_yog","ex1_cheese","ex1_egg","ex1_fish","ex1_fruit",
                      "ex1_poultry","ex1_tofu","ex1_meat","ex1_pork","ex1_starch","ex1_cornflakes","ex1_whole_grain","ex1_patat",
                      "ex1_lentil","ex1_vegetable","ex1_salad","ex1_walnut","ex1_sugar","ex1_chips","ex1_pork")
  
  if(controls$ffq_grouping) {
    # option 1: categorize 
    tmp = dat_4 %>% 
      dplyr::select(all_of(ex1_dichotomize)) %>% 
      dplyr::mutate(across(everything(), ~ case_when(.%in% 1:3 ~ "Less than weekly",
                                                     . %in% 4:5 ~ "At least weekly",
                                                     . %in% 6:8 ~ "Daily and more",
                                                     TRUE ~ NA_character_,
                                                     .ptype = factor(levels = c("Less than weekly","At least weekly","Daily and more"),ordered=FALSE)))) %>% 
      dplyr::rename_with(~paste0(.x,"_weekly"))
  } else {
    # option 2: continuous
    tmp = dat_4 %>% 
      dplyr::select(all_of(ex1_dichotomize)) %>% 
      # dplyr::mutate(across(everything(), ~ factor(.,ordered=FALSE))) %>% 
      dplyr::rename_with(~paste0(.x,"_freq"))
  }
  
  dat_4 = dat_4 %>% 
    dplyr::bind_cols(tmp) %>% 
    dplyr::select(-all_of(ex1_dichotomize))
  
  # EX1 Hot meals from disposable packages
  
  # 1, More than 7 meals per week | 2, 5-7 meals per week | 3, 2-4 meals per week | 4, 1 meal per week | 5, 1-2 meals per month | 6, Less than one meal per month | 7, Never
  tmp = dat_4 %>% 
    dplyr::select(ex1_hot_meals_package) %>% 
    dplyr::mutate(across(everything(), ~ case_when(. %in% 1:4 ~ "At least weekly",
                                                   . %in% 5:7 ~ "Less than weekly",
                                                   TRUE ~ NA_character_,
                                                   .ptype = factor(levels = c("Less than weekly","At least weekly"),ordered=FALSE)))) %>% 
    dplyr::rename_with(~paste0(.x,"_freq"))
  dat_4 = dat_4 %>% 
    dplyr::select(-ex1_hot_meals_package) %>% 
    dplyr::bind_cols(tmp) 
  
  # EX1: cosmetics
  
  # 1, Never | 2, Less than once per month | 3, 1-3 times per month | 4, 1-3 times per week | 5, 4-7 times per week | 6, More than once per day
  ex1_dichotomize2 = c("ex1_deo_spray","ex1_deo_roll","ex1_bodylotion","ex1_handcream","ex1_perfume","ex1_hairspray","ex1_toothpaste")
  
  tmp = dat_4 %>% 
    dplyr::select(all_of(ex1_dichotomize2)) %>% 
    dplyr::mutate(across(everything(), ~ case_when(. %in% 1:3 ~ "Less than weekly",
                                                   . %in% 4:6 ~ "At least weekly",
                                                   TRUE ~ NA_character_,
                                                   .ptype = factor(levels = c("Less than weekly","At least weekly"),ordered=FALSE)))) %>% 
    dplyr::rename_with(~paste0(.x,"_freq"))
  
  dat_4 = dat_4 %>% 
    dplyr::select(-all_of(ex1_dichotomize2)) %>% 
    dplyr::bind_cols(tmp) 
  
  # EX1: fish types
  dat_4 = dat_4 %>% 
    dplyr::mutate(ex1_tuna=factor(ex1_tuna,levels=0:1,labels=c("No/rarely tuna","Frequent consumption of tuna")),
                  ex1_salmon=factor(ex1_salmon,levels=0:1,labels=c("No/rarely salmon","Frequent consumption of salmon")),
                  ex1_pangasius=factor(ex1_pangasius,levels=0:1,labels=c("No/rarely pangasius","Frequent consumption of pangasius")),
                  ex1_trout=factor(ex1_trout,levels=0:1,labels=c("No/rarely trout","Frequent consumption of trout")),
                  ex1_shrimp=factor(ex1_shrimp,levels=0:1,labels=c("No/rarely shrimp","Frequent consumption of shrimp"))
    )
  
  # EX1: categorize residence 
  dat_4 = dat_4 %>% 
    dplyr::mutate(ex1_residence_urban=case_when(ex1_residence_city==1 | ex1_residence_suburb==1 | ex1_residence_industrial==1 ~ "Residence in urban/industrial area",
                                                ex1_residence_countryside==1 | ex1_residence_nearfield==1 ~ "Residence in rural area"),
                  ex1_residence_urban=as.factor(ex1_residence_urban)) %>% 
    dplyr::select(-ex1_residence_city,-ex1_residence_suburb,-ex1_residence_industrial,-ex1_residence_countryside,-ex1_residence_nearfield)
  
  # Ex1: diet 
  if(TRUE) {
  dat_4 = dat_4 %>% 
    dplyr::mutate(ex1_veggie_vegan=case_when(ex1_q41_diet==1 ~ "Omnivorous diet",
                                             ex1_q41_diet %in% 2:4 ~ "Vegetarian/vegan diet"),
                  ex1_veggie_vegan=as.factor(ex1_veggie_vegan)) %>% 
    dplyr::select(-ex1_q41_diet)
  }
  
  # EX1: tap water
  dat_4 = dat_4 %>% 
    dplyr::mutate(ex1_tapwater=case_when(ex1_tapwater %in% 1:2 ~ "More than 2L of tap water per day",
                                         ex1_tapwater %in% 3:8 ~ "Up to 2L of tap water per day")) %>% 
    dplyr::mutate(ex1_tapwater=factor(ex1_tapwater,levels=c("Up to 2L of tap water per day","More than 2L of tap water per day")))
  
  # EX1: alcool
  dat_4 = dat_4 %>% 
    dplyr::mutate(ex1_alcool_weekly=case_when(ex1_alcool %in% 1:6 ~ "Daily and more",
                                              ex1_alcool %in% 7 ~ "At least weekly",
                                              ex1_alcool %in% 8 ~ "Less than weekly",
                                              TRUE ~ NA_character_,
                                              .ptype = factor(levels = c("Less than weekly","At least weekly","Daily and more"),ordered=FALSE))) %>% 
    dplyr::select(-ex1_alcool)
  
  # EX1: carpet 
  dat_4 = dat_4 %>%
    dplyr::mutate(ex1_carpet=if_else(ex1_q10_q1_lounge %in% 7:8 | ex1_q10_q2_bedroom %in% 7:8 | ex1_q10_q3_rest %in% 7:8 ,"Synthetic or natural carpet","No carpet"),
                  ex1_carpet=as.factor(ex1_carpet)) %>% 
    dplyr::select(-ex1_q10_q1_lounge,-ex1_q10_q2_bedroom,-ex1_q10_q3_rest)
  
  # EX1: swiss fish (comes from manual recoding of free text entries by Céline Fragnière)
  dat_4 = dat_4 %>% 
    dplyr::mutate(ex1_swiss_fish=if_else(ex1_swiss_fish=="oui","Consumption of fish from Swiss lakes","No consumption of fish from Swiss lakes"),
                  ex1_swiss_fish=factor(ex1_swiss_fish,levels=c("No consumption of fish from Swiss lakes","Consumption of fish from Swiss lakes")))
  
  # EX1: categorize distance to nearest field/vineyard
  if(FALSE) {
    dat_4 = dat_4 %>% 
      dplyr::mutate(ex1_distance_field=if_else(ex1_distance_field<200,1,0,0)) %>% 
      dplyr::mutate(ex1_distance_field=factor(ex1_distance_field,levels=c(0,1),labels=c("More than 200m to nearest field/vineyard","Less than 200m to nearest field/vineyard")))
  }
  
  # EX1: tobacco
  tmp = dat_4 %>% 
    dplyr::select(age,starts_with("ex1_q28_q"),starts_with("ex1_q29_q")) %>% 
    dplyr::mutate(across(ends_with("_start"),~as.numeric(as.character(.)))) %>% 
    dplyr::mutate(across(ends_with("_stop"),~as.numeric(as.character(.)))) %>% 
    ## one person starting at 99
    dplyr::mutate(ex1_q28_q1_cig_start=ifelse(ex1_q28_q1_cig_start>age,NA,ex1_q28_q1_cig_start), 
                  ex1_q28_q2_cigarillos_start=ifelse(ex1_q28_q2_cigarillos_start>age,NA,ex1_q28_q2_cigarillos_start)) %>% 
    ## cigs per day
    dplyr::mutate(number_day = as.numeric(as.character(dat_in$ex1_q29_q1_q1_cigarettes)),
                  number_day = ifelse(number_day==11315,NA,number_day), # remove mistake
                  number_week = as.numeric(as.character(dat_in$ex1_q29_q1_q2_cigarettes)),
                  number_day = ifelse(is.na(number_day) & number_week>7,number_week/7,number_day)) %>% 
    ## cigarettes, cigars and pipe
    dplyr::mutate(ex1_cig=case_when(ex1_q28_q1_cig==2 | ex1_q28_q2_cigarillos==2 | ex1_q28_q3_pipe==2 ~ "Current smoker",
                                    ex1_q28_q1_cig==3 | ex1_q28_q2_cigarillos==3 | ex1_q28_q3_pipe==3 ~ "Past smoker",
                                    TRUE ~ "Never smoked"),
                  ex1_cig=factor(ex1_cig,levels=c("Never smoked","Past smoker","Current smoker"),ordered=FALSE),
                  cig_start=pmin(ex1_q28_q1_cig_start,ex1_q28_q2_cigarillos_start,ex1_q28_q3_pipe_start,na.rm=TRUE),
                  cig_stop=pmax(age,ex1_q28_q1_cig_stop,ex1_q28_q2_cigarillos_stop,ex1_q28_q3_pipe_stop,na.rm=TRUE),
                  ex1_cig_years=cig_stop-cig_start,
                  ex1_cig_packyears=number_day/20*ex1_cig_years) %>% 
    ## other tobacco products
    dplyr::mutate(ex1_other_tobacco=case_when(ex1_q28_q4_ecig==2 | ex1_q28_q5_heat==2 | ex1_q28_q6_shisha==2 | ex1_q28_q7_other2==2 ~ "Current user of other tobacco products",
                                              ex1_q28_q4_ecig==3 | ex1_q28_q5_heat==3 | ex1_q28_q6_shisha==3 | ex1_q28_q7_other2==3  ~ "Past user of other tobacco products",
                                              TRUE ~ "Never used other tobacco products"),
                  ex1_other_tobacco=factor(ex1_other_tobacco,levels=c("Never used other tobacco products","Past user of other tobacco products","Current user of other tobacco products"),ordered=FALSE),
                  cig_start=pmin(ex1_q28_q4_ecig_start,ex1_q28_q5_heat_start,ex1_q28_q6_shisha_start,ex1_q28_q7_other_start,na.rm=TRUE),
                  cig_start=ifelse(cig_start>age,age,cig_start),
                  cig_stop=pmax(age,ex1_q28_q4_ecig_stop,ex1_q28_q6_shisha_stop,ex1_q28_q7_other_stop,na.rm=TRUE),
                  ex1_other_tobacco_years=cig_stop-cig_start) %>% 
    ## replace NA by 0 for duration of smoking of non smokers
    dplyr::mutate(ex1_cig_years=ifelse(ex1_cig_years==0 & ex1_cig!="Never smoked",NA,ex1_cig_years),
                  ex1_cig_packyears=ifelse(ex1_cig_packyears==0 & ex1_cig!="Never smoked",NA,ex1_cig_packyears),
                  ex1_cig_years=ifelse(ex1_cig=="Never smoked",0,ex1_cig_years),
                  ex1_cig_packyears=ifelse(ex1_cig=="Never smoked",0,ex1_cig_packyears)) %>% 
    dplyr::select(ex1_cig,ex1_cig_years,ex1_cig_packyears,ex1_other_tobacco,ex1_other_tobacco_years)
  
  dat_4 = dat_4 %>% 
    dplyr::bind_cols(tmp) %>% 
    dplyr::select(-starts_with("ex1_q28_q"),-starts_with("ex1_q29_q"))
  
  # CORONA: time since past disease
  dat_5 = dat_4 %>% 
    dplyr::mutate(corona_time_since_past_disease = as.numeric(date - corona_past_disease_date)) %>% 
    dplyr::relocate(corona_time_since_past_disease,.after=corona_past_disease_date) 
  
  # EX2: rename
  dat_6 = dat_5 %>% 
    dplyr::mutate(ex2_garden_pesticides=factor(ex2_garden_pesticides,levels=0:1,labels=c("No exposition to pesticides","Exposition to pesticides")),
                  ex2_skiwax=factor(ex2_skiwax,levels=0:1,labels=c("No exposition to ski wax","Exposition to ski wax")),
                  ex2_spray_imp=factor(ex2_spray_imp,levels=0:1,labels=c("No exposition to impregnation sprays","Exposition to impregnation sprays")),
                  ex2_prof_smoke=factor(ex2_prof_smoke,levels=0:1,labels=c("No professional exposition to smoke","Professional exposition to smoke")),
                  ex2_prof_exhaust_gas=factor(ex2_prof_exhaust_gas,levels=0:1,labels=c("No professional exposition to exhaust gas","Professional exposition to exhaust gas")),
                  ex2_prof_solvent_vapours=factor(ex2_prof_solvent_vapours,levels=0:1,labels=c("No professional exposition to solvent vapours","Professional exposition to solvent vapours")),
                  ex2_prof_cleaning_agents=factor(ex2_prof_cleaning_agents,levels=0:1,labels=c("No professional exposition to cleaning agents","Professional exposition to cleaning agents")),
                  ex2_prof_dust=factor(ex2_prof_dust,levels=0:1,labels=c("No professional exposition to dust","Professional exposition to dust")),
                  ex2_prof_polycarbonate	=factor(ex2_prof_polycarbonate	,levels=0:1,labels=c("No professional exposition to polycarbonate","Professional exposition to polycarbonate")),
                  ex2_prof_epoxy_resin=factor(ex2_prof_epoxy_resin,levels=0:1,labels=c("No professional exposition to epoxy resin","Professional exposition to epoxy resin")),
                  ex2_prof_skiwax=factor(ex2_prof_skiwax,levels=0:1,labels=c("No professional exposition to ski wax","Professional exposition to ski wax")),
                  ex2_prof_spray_imp=factor(ex2_prof_spray_imp,levels=0:1,labels=c("No professional exposition to impregnation sprays","Professional exposition to impregnation sprays")),
    )
  
  # BLOODURINE: no new variable
  dat_7 = dat_6
  
  # PFAS: no new variable
  dat_8 = dat_7
  
  # Save only last stage
  dat_out = dat_8
  write_rds(dat_out,file=file.path(controls$savepoint,"shesp_3.rds"))
  
  return(dat_out)
}
