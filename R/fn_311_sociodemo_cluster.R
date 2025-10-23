#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 311 - figure for cluster by sociodemographic variable
#' author: Julien Riou
#' date: 2025-05-01
#' ---

fn_311_sociodemo_cluster = function(dat_in, ...) {
  # labels
  labs <- list(`age_group[20,30)`="Age group: 20-29",
               `age_group[30,40)`="Age group: 30-39",
               `age_group[40,50)`="Age group: 40-49",
               `age_group[50,60)`="Age group: 50-59",
               `age_group[60,70]`="Age group: 60-69",
               genderFemales='Sex: female',
               genderMales='Sex: male',
               `a1_foreignSwiss nationality`='Nationality: Swiss',
               `a1_foreignForeign nationality`='Nationality: non-Swiss',
               `a1_education_3classesUp to high school` = "Education: up to high school",
               `a1_education_3classesApprenticeship / professional baccalaureate` = "Education: apprenticeship",
               a1_education_3classesUniversity ="Education: university",
               `a1_household_monthly_income_3classes<CHF 4,500`="Income: low", 
               `a1_household_monthly_income_3classesCHF 4,500 to 9,000`="Income: intermediate", 
               `a1_household_monthly_income_3classes>CHF 9,000`="Income: high",
               center_labelBern="Center: Bern",
               center_labelLausanne="Center: Lausanne")
  corr_cova_labels = tibble(parameter=names(labs),
                            parameter_lab=unname(unlist(labs))) %>% 
    dplyr::mutate(parameter_lab=factor(parameter_lab,levels=unname(unlist(labs)))) %>% 
    dplyr::mutate(parameter_lab_cat = c(rep("Age group",5),
                                        rep("Sex",2),
                                        rep("Nationality",2),
                                        rep("Education",3),
                                        rep("Income",3),
                                        rep("Center",2)    )) %>%
    dplyr::mutate(parameter_lab_cat = fct_inorder(parameter_lab_cat))
  
  dd_ = dat_in %>% 
    dplyr::select(age_group,
                  gender,
                  a1_foreign,
                  a1_education_3classes,
                  a1_household_monthly_income_3classes,
                  center_label,
                  pfas_cluster)  %>% 
    pivot_longer(1:6) %>% 
    group_by(name,value,pfas_cluster) %>% 
    summarise(n=n(), .groups="drop_last") %>%
    filter(!is.na(value)) %>% 
    mutate(tot=sum(n),
           p=n/tot,
           pfas_cluster=factor(pfas_cluster,levels=c("3_High","2_Intermediate","1_Low"),labels=rev(controls$lab_pfas_cluster)),
           parameter=paste0(name,value)) %>% 
    left_join(corr_cova_labels,by = join_by(parameter))
  
  g = dd_ %>% 
    ggplot() +
    geom_col(aes(x=p,y=parameter_lab,fill=pfas_cluster),width=1,
             colour="black") +
    scale_fill_manual(values=rev(controls$col_pfas_cluster)) +
    scale_x_continuous( labels = scales::label_percent(accuracy = 1), expand=expansion(0,0))+ 
    scale_y_discrete(limits=rev, expand=expansion(0,0))+ 
    facet_grid(parameter_lab_cat~.,scales = "free_y",space = "free_y")  +
    labs(x="Proportion",y=NULL,fill="PFAS cluster:") +
    theme(legend.position="bottom",
          legend.direction="vertical",
          strip.text = element_blank(),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), "cm")) +
    guides(fill = guide_legend(reverse = TRUE)) 
  
  return(g)
}

fn_311_center_differences = function(dat_in, vars, ...) {
  # labels
  labs <- list(cluster1_balanced="Cluster: balanced",
               cluster2_dairy_focused="Cluster: dairy-focused",
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
               ex1_hot_meals_package_freq="Disposable food packaging",
               ex1_residence_urbanResidenceinurbanDindustrialarea='Residence in urban/industrial area ',
               ex1_carpetSyntheticornaturalcarpet="Synthetic or natural carpet at home",
               
               ex1_deo_spray_freq="Deodorant spray",
               ex1_deo_roll_freq="Deodorant roll",
               ex1_bodylotion_freq="Body lotion",
               ex1_handcream_freq="Hand cream",
               ex1_perfume_freq="Perfume",
               ex1_hairspray_freq="Hair spray",
               
               ex2_prof_smokeProfessionalexpositiontosmoke='Smoke at work',
               ex2_prof_exhaust_gasProfessionalexpositiontoexhaustgas='Exhaust gas at work',
               ex2_prof_solvent_vapoursProfessionalexpositiontosolventvapours='Solvent vapours at work',
               ex2_prof_cleaning_agentsProfessionalexpositiontocleaningagents='Cleaning agents at work',
               ex2_prof_dustProfessionalexpositiontodust='Dust at work',
               ex2_prof_skiwaxProfessionalexpositiontoskiwax='Ski wax at work',
               ex2_prof_spray_impProfessionalexpositiontoimpregnationsprays='Impregnation sprays at work'
  )
  corr_pfas_labels = tibble(dependent_variable=c("mu2","mu3"),
                            dependent_variable_lab=paste("Cluster",controls$lab_pfas_cluster[2:3]))
  corr_cova_labels = tibble(parameter=names(labs2),
                            parameter_lab=unname(unlist(labs2))) %>% 
    dplyr::mutate(parameter_lab=factor(parameter_lab,levels=unname(unlist(labs2)))) %>% 
    dplyr::mutate(parameter_lab_cat = c(rep("Socio-demographics",8),
                                        rep("Diet cluster",4),
                                        rep("Seafood",7),
                                        rep("Diet and lifestyle",6),
                                        rep("Environment",5),
                                        rep("Cosmetics",6),
                                        rep("Occupational",7)
    )) %>%
    dplyr::mutate(parameter_lab_cat = fct_inorder(parameter_lab_cat))
  
  dd_ = dat_in %>% 
    dplyr::select(center_label,
                  all_of(names(vars)))  %>% 
    dplyr::mutate(across(everything(), function(x) as.factor(x))) %>% 
    pivot_longer(2:(length(vars)+1)) %>% 
    group_by(name,value,center_label) %>% 
    summarise(n=n(), .groups="drop_last") %>%
    filter(!is.na(value)) %>% 
    mutate(tot=sum(n),
           p=n/tot,
           # pfas_cluster=factor(pfas_cluster,levels=c("3_High","2_Intermediate","1_Low"),labels=rev(controls$lab_pfas_cluster)),
           parameter=paste0(name,value)) %>% 
    ungroup() %>% 
    left_join(corr_cova_labels,by = join_by(parameter))
  
  g = dd_ %>% 
    ggplot() +
    geom_col(aes(x=p,y=parameter_lab,fill=center_label),width=1,
             colour="black") +
    scale_fill_manual(values=rev(controls$col_pfas_cluster)) +
    scale_x_continuous( labels = scales::label_percent(accuracy = 1), expand=expansion(0,0))+ 
    scale_y_discrete(limits=rev, expand=expansion(0,0))+ 
    facet_grid(parameter_lab_cat~.,scales = "free_y",space = "free_y")  +
    labs(x="Proportion",y=NULL,fill="PFAS cluster:") +
    theme(legend.position="bottom",
          legend.direction="vertical",
          strip.text = element_blank(),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), "cm"))
  
  return(g)
}
