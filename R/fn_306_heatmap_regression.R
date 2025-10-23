#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 306 - heatmap for regression coefficients
#' author: Julien Riou
#' date: 2023-11-01
#' ---

fn_306_heatmap_regression = function(coefs, type="FC", uplim=NULL, age_per10=TRUE, bart=NULL, horseshoe=NULL, ...) {
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=unname(unlist(controls$labs))))
  
  if(age_per10) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(parameter=="age",exp(10*Estimate),exp_beta),
                    exp_beta_lob =ifelse(parameter=="age",exp(10*`l-95% CI`),exp_beta_lob),
                    exp_beta_upb =ifelse(parameter=="age",exp(10*`u-95% CI`),exp_beta_upb))
  }
  if(!is.null(uplim)) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(exp_beta>=uplim,uplim-.01,exp_beta),
                    exp_beta=ifelse(exp_beta<1/uplim,1/uplim+.01,exp_beta))
  }
  
  # selected from bart
  if(!is.null(bart)) {
    sel_bart = bart %>% 
      dplyr::select(parameter,dependent_variable) %>% 
      dplyr::mutate(sel_bart=TRUE) %>% 
      dplyr::filter(parameter!="Intercept")
    coefs = coefs %>% 
      dplyr::left_join(sel_bart,by = join_by(parameter, dependent_variable))
  }
  # selected from horseshoe
  if(!is.null(horseshoe)) {
    sel_horseshoe = horseshoe %>% 
      dplyr::mutate(sel_horseshoe=exp_beta_upb<1 | exp_beta_lob>1) %>% 
      dplyr::filter(sel_horseshoe) %>% 
      dplyr::select(parameter,dependent_variable,sel_horseshoe) %>% 
      dplyr::filter(parameter!="Intercept")
    coefs = coefs %>% 
      dplyr::left_join(sel_horseshoe,by = join_by(parameter, dependent_variable))
  }
  
  
  # format results
  dd_ = coefs %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<=1 | exp_beta_lob>=1,"95% CrI away from 1","95% CrI includes 1")) 
  # dplyr::mutate(away=ifelse(exp_beta<=.98 | exp_beta>=1.02,"95% CrI away from 1","95% CrI includes 1"))
  
  if(type=="OR") {
    dd_ =  dd_ %>% 
      filter(type %in% c("Odds ratio of detection","Odds ratio"))
    if(is.null(uplim)) {
      uplim = round(max(2,dd_$exp_beta),2)
      botlim = round(min(.5,dd_$exp_beta),2)  
    } else {
      botlim = round(1/uplim,2)
    }
    
    g1 = dd_ %>% 
      dplyr::mutate(exp_beta=ifelse(exp_beta<botlim,botlim,exp_beta),
                    exp_beta=ifelse(exp_beta>uplim,uplim,exp_beta)) %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
      # geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,colour=factor(away)),fill="transparent") +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
      scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=controls$pfas_order) +
      scale_colour_manual(values=c("black","transparent"),guide="none") +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom",
            legend.key.width = unit(2, "cm"),
            background.grid=element_blank()) +
      labs(x=NULL,y="PFAS",fill="Odds ratio") 
    
    g2 = dd_ %>% 
      filter(away=="95% CrI away from 1") %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
      geom_text(aes(x=parameter2_lab,y=dependent_variable_lab,label=sprintf("%.2f",exp_beta))) +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +        scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=pfas_order) +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom",
            legend.key.width = unit(2, "cm"),
            background.grid=element_blank()) +
      labs(x=NULL,y="PFAS",fill="Odds ratio") 
    
    if(!is.null(bart)) {
      g3 = dd_ %>% 
        dplyr::filter(sel_bart|sel_horseshoe) %>% 
        filter(away=="95% CrI away from 1") %>% 
        ggplot() +
        geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
        geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,colour=factor(away)),fill="transparent") +
        # geom_text(aes(x=parameter2_lab,y=dependent_variable_lab,label=sprintf("%.2f",exp_beta))) +
        # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
        scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                             mid="white",midpoint=0,low="dodgerblue",high="firebrick") +        scale_x_discrete(expand=expansion(0,0)) +
        scale_y_discrete(expand=expansion(0,0),limits=pfas_order) +
        scale_colour_manual(values=c("black","transparent"),guide="none") +
        theme(axis.text.x=element_text(angle=45,hjust=1),
              legend.position="bottom",
              # legend.key.width = unit(2, "cm"),
              background.grid=element_blank()) +
        labs(x=NULL,y="PFAS",fill="Odds ratio") 
    } else {
      g3 = NULL
    }
  }
  
  if(type=="FC") {
    dd_ =  dd_ %>% 
      filter(type=="Fold-change")
    if(is.null(uplim)) {
      uplim = round(max(2,dd_$exp_beta),2)
      botlim = round(min(.5,dd_$exp_beta),2)  
    } else {
      botlim = round(1/uplim,2)
    }
    
    g1 = dd_ %>% 
      dplyr::mutate(exp_beta=ifelse(exp_beta<botlim,botlim,exp_beta),
                    exp_beta=ifelse(exp_beta>uplim,uplim,exp_beta)) %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
      # geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,colour=factor(away)),fill="transparent") +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
      scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=pfas_order) +
      scale_colour_manual(values=c("black","transparent"),guide="none") +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom",
            legend.key.width = unit(2, "cm"),
            panel.grid = element_blank()) +
      labs(x=NULL,y="PFAS",fill="Fold-change") 
    
    g2 = dd_ %>% 
      filter(away=="95% CrI away from 1") %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
      # geom_text(aes(x=parameter2_lab,y=dependent_variable_lab,label=sprintf("%.2f",exp_beta))) +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
      scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=pfas_order) +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom",
            legend.key.width = unit(2, "cm"),
            panel.grid = element_blank()) +
      labs(x=NULL,y="PFAS",fill="Fold-change") 
    
    if(!is.null(bart)) {
      g3 = dd_ %>% 
        dplyr::filter(sel_bart|sel_horseshoe) %>% 
        filter(away=="95% CrI away from 1") %>% 
        ggplot() +
        geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
        geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,colour=factor(away)),fill="transparent") +
        # geom_text(aes(x=parameter2_lab,y=dependent_variable_lab,label=sprintf("%.2f",exp_beta))) +
        # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
        scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                             mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
        scale_x_discrete(expand=expansion(0,0)) +
        scale_y_discrete(expand=expansion(0,0),limits=pfas_order) +
        scale_colour_manual(values=c("black","transparent"),guide="none") +
        theme(axis.text.x=element_text(angle=45,hjust=1),
              legend.position="bottom",
              # legend.key.width = unit(2, "cm"),
              panel.grid = element_blank()) +
        labs(x=NULL,y="PFAS",fill="Fold-change") 
    } else {
      g3 = NULL
    }
    
  }
  return(list(g1,g2,g3))
}


fn_306_heatmap_regression2 = function(coefs, type="FC", uplim=NULL, age_per10=TRUE, bart=NULL, horseshoe=NULL, ...) {
  # labels
  labs2 <- list(age="Age",
                genderMales='Sex: male',
                a1_foreignForeignnationality='Nationality: non-swiss',
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
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(labs2),
                            parameter2_lab=unname(labs2)) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=unname(unlist(labs2))))
  
  if(age_per10) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(parameter=="age",exp(10*Estimate),exp_beta),
                    exp_beta_lob =ifelse(parameter=="age",exp(10*`l-95% CI`),exp_beta_lob),
                    exp_beta_upb =ifelse(parameter=="age",exp(10*`u-95% CI`),exp_beta_upb))
  }
  if(!is.null(uplim)) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(exp_beta>=uplim,uplim-.01,exp_beta),
                    exp_beta=ifelse(exp_beta<1/uplim,1/uplim+.01,exp_beta))
  }
  
  # selected from bart
  if(!is.null(bart)) {
    sel_bart = bart %>% 
      dplyr::select(parameter,dependent_variable) %>% 
      dplyr::mutate(sel_bart=TRUE) %>% 
      dplyr::filter(parameter!="Intercept")
    coefs = coefs %>% 
      dplyr::left_join(sel_bart,by = join_by(parameter, dependent_variable))
  }
  
  # format results
  dd_ = coefs %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<=1 | exp_beta_lob>=1,"95% CrI away from 1","95% CrI includes 1")) 
  # dplyr::mutate(away=ifelse(exp_beta<=.98 | exp_beta>=1.02,"95% CrI away from 1","95% CrI includes 1"))
  
  if(type=="OR") {
    dd_ =  dd_ %>% 
      filter(type %in% c("Odds ratio of detection","Odds ratio")) 
    if(is.null(uplim)) {
      uplim = round(max(2,dd_$exp_beta),2)
      botlim = round(min(.5,dd_$exp_beta),2)  
    } else {
      botlim = round(1/uplim,2)
    }
    
    g3 = dd_ %>% 
      dplyr::filter(sel_bart) %>% 
      filter(away=="95% CrI away from 1") %>% 
      ggplot() +
      geom_tile(aes(y=parameter2_lab,x=dependent_variable_lab,fill=exp_beta),colour="black") +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +       
      scale_x_discrete(expand=expansion(0,0),limits=pfas_order) +
      scale_y_discrete(expand=expansion(0,0),limits=rev) +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom") +
      labs(x=NULL,y="PFAS",fill="Odds ratio") 
  }
  
  if(type=="FC") {
    dd_ =  dd_ %>% 
      filter(type=="Fold-change")
    if(is.null(uplim)) {
      uplim = round(max(2,dd_$exp_beta),2)
      botlim = round(min(.5,dd_$exp_beta),2)  
    } else {
      botlim = round(1/uplim,2)
    }
    
    g3 = dd_ %>% 
      dplyr::filter(sel_bart) %>% 
      filter(away=="95% CrI away from 1") %>% 
      ggplot() +
      geom_tile(aes(y=parameter2_lab,x=dependent_variable_lab,fill=exp_beta),colour="black") +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
      scale_x_discrete(expand=expansion(0,0),limits=pfas_order) +
      scale_y_discrete(expand=expansion(0,0),limits=rev) +
      scale_colour_manual(values=c("black","transparent"),guide="none") +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom") +
      labs(x=NULL,y="PFAS",fill="Fold-change") 
    
  }
  return(g3)
}


fn_306_heatmap_regression3 = function(coefs, uplim=NULL, age_per10=TRUE, bart=NULL, horseshoe=NULL, keepall=FALSE, ...) {
  # labels
  labs2 <- list(age="Age",
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
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels2) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels2))
  corr_cova_labels = tibble(parameter2=names(labs2),
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
  
  if(age_per10) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(parameter=="age",exp(10*Estimate),exp_beta),
                    exp_beta_lob =ifelse(parameter=="age",exp(10*`l-95% CI`),exp_beta_lob),
                    exp_beta_upb =ifelse(parameter=="age",exp(10*`u-95% CI`),exp_beta_upb))
  }
  if(!is.null(uplim)) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(exp_beta>=uplim,uplim-.01,exp_beta),
                    exp_beta=ifelse(exp_beta<1/uplim,1/uplim+.01,exp_beta))
  }
  
  # selected from bart
  if(!is.null(bart)) {
    sel_bart = bart %>% 
      dplyr::select(parameter,dependent_variable) %>% 
      dplyr::mutate(sel_bart=TRUE) %>% 
      dplyr::filter(parameter!="Intercept")
    coefs = coefs %>% 
      dplyr::left_join(sel_bart,by = join_by(parameter, dependent_variable))
  }
  

  # format results
  dd_ = coefs %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<=1 | exp_beta_lob>=1,"95% CrI away from 1","95% CrI includes 1")) %>% 
    dplyr::mutate(away=factor(away, levels=c("95% CrI includes 1","95% CrI away from 1")),
                  type=factor(type, levels=c("Odds ratio of detection","Fold-change"),
                              labels=c("Odds-ratio of detection",
                                       "Fold-change upon detection")),
                  dependent_variable_lab=factor(dependent_variable_lab,levels=pfas_order)
    )
  if(!keepall) {
    dd_ = dd_ %>% 
      dplyr::filter(away=="95% CrI away from 1") 
  }
  
  uplim = round(max(2,dd_$exp_beta),2)
  botlim = round(min(.5,dd_$exp_beta),2)  
  
  # Plot
  g = ggplot(dd_) +
    geom_tile(aes(y=parameter_lab,x=dependent_variable_lab,fill=exp_beta)) +
    geom_tile(data=filter(dd_, sel_bart),
              aes(y=parameter_lab,x=dependent_variable_lab,fill=exp_beta),colour="black",linewidth=.5) +
    scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                         mid="white",midpoint=0,low="dodgerblue",high="firebrick") +       
    scale_x_discrete(expand=expansion(0,add=0.7),limits=rev) +
    scale_y_discrete(expand=expansion(0,add=.7),limits=rev) +
    facet_grid(parameter_lab_cat~type,scales = "free",space = "free") +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.position="bottom",
          legend.key.size = unit(1,"cm")) +
    labs(x="PFAS",y="",fill="Estimate (OR or FC):", colour="")   +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        label.position = "bottom",
        ticks = FALSE,
        # Make the key *wide*; adjust if needed for your device size
        barwidth = unit(5, "cm"),
        barheight = unit(.3, "cm")
      )
    ) 
  
  
return(g)
}

