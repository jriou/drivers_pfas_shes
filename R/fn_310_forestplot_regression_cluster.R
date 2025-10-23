#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 310 - forestplot for regression coefficients
#' author: Julien Riou
#' date: 2025-05-01
#' ---

fn_310_forestplot_regression_cluster = function(coefs, uplim=11, age_per10=TRUE, bart=NULL, horseshoe=NULL, bart_stars=NULL, ...) {
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
  corr_pfas_labels = tibble(dependent_variable=c("mu2","mu3"),
                            dependent_variable_lab=controls$lab_pfas_cluster[2:3])
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
  
  if(!is.null(bart_stars)) {
    corr_cova_labels = corr_cova_labels %>% 
      dplyr::mutate(bart_star = ifelse(parameter %in% bart_stars,"Yes","No"))
  }
  
  if(age_per10) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(parameter=="age",exp(10*Estimate),exp_beta),
                    exp_beta_lob =ifelse(parameter=="age",exp(10*`l-95% CI`),exp_beta_lob),
                    exp_beta_upb =ifelse(parameter=="age",exp(10*`u-95% CI`),exp_beta_upb))
  }
  if(!is.null(uplim)) {
    coefs = coefs %>% 
      dplyr::mutate(exp_beta=ifelse(exp_beta>=uplim,NA,exp_beta),
                    exp_beta_lob=ifelse(exp_beta_lob>=uplim,uplim-.01,exp_beta_lob),
                    exp_beta_upb=ifelse(exp_beta_upb>=uplim,uplim-.01,exp_beta_upb),
                    exp_beta=ifelse(exp_beta<1/uplim,1/uplim+.01,exp_beta),
                    exp_beta_lob=ifelse(exp_beta_lob<1/uplim,1/uplim+.01,exp_beta_lob),
                    exp_beta_upb=ifelse(exp_beta_upb<1/uplim,1/uplim+.01,exp_beta_upb))
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
    dplyr::left_join(corr_cova_labels,by = join_by(parameter),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<=1 | exp_beta_lob>=1,"95% CrI away from 1","95% CrI includes 1")) 
  # dplyr::mutate(away=ifelse(exp_beta<=.98 | exp_beta>=1.02,"95% CrI away from 1","95% CrI includes 1"))
  
  # lines placement
  cat_placement = corr_cova_labels %>% 
    group_by(parameter_lab_cat) %>% 
    summarise(n=n()) %>% 
    mutate(x_end=cumsum(n)+.5,
           x_start=x_end-n,
           midpoint=(x_end-x_start)/2)
  
    g = dd_ %>%
      ggplot() +
      geom_hline(yintercept=1,linetype=2) +
      geom_pointrange(aes(y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb,x=parameter_lab,
                          colour=dependent_variable_lab,shape=bart_star),size=.5,
                      position=position_dodge(.5)) +
      # geom_vline(xintercept=cat_placement$x_end) +
      facet_grid(parameter_lab_cat~.,scales = "free_y",space = "free_y") +
      scale_y_log10(expand=expansion(0,0),breaks=c(.1,.2,.5,1,2,5,10),minor_breaks =NULL) +
      scale_x_discrete(limits=rev,expand=expansion(0,add=0.5)) +
      scale_colour_manual(values=controls$col_pfas_cluster[2:3]) +
      theme(legend.position="bottom",
            legend.direction="vertical", legend.margin=margin()) +
      labs(x=NULL,
           y=paste("aOR"),
           colour="PFAS cluster:",
           shape="BART top 10:") +
      coord_flip()

  return(g)
}
