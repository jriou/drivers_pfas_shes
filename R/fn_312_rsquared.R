#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 312 - figure for r squared and mediation
#' author: Julien Riou
#' date: 2025-05-01
#' ---

fn_312_rsquared = function(rr_, ...) {
  
  r2_out_labs = tibble(model=c("1.1","1.2","1.3","1.4","1.5","1.6"),
                       adj = c("Age + sex",
                               "Socio-demographic\n(age + sex + education +\n income + nationality)",
                               "Socio-demographic + center",
                               "BART top 10*",
                               "BART top 10 + socio-demographic",
                               "Socio-demographic + center + diet cluster")) %>% 
    dplyr::mutate(adj2=factor(adj, levels=c("Age + sex",
                                            "Socio-demographic\n(age + sex + education +\n income + nationality)",
                                            "Socio-demographic + center",
                                            "Socio-demographic + center + diet cluster",
                                            "BART top 10*",
                                            "BART top 10 + socio-demographic")))
  
  g1 = rr_ %>% 
    dplyr::left_join(r2_out_labs, by="model") %>% 
    dplyr::filter(!(model %in% c("1.5"))) %>% 
    ggplot() +
    geom_col(aes(x=adj2,y=Estimate,fill=as.numeric(adj2)),colour="black") +
    geom_errorbar(aes(x=adj2,ymin=Q2.5,ymax=Q97.5), width=0) +
    scale_fill_viridis_c(option="C",guide="none") +      
    scale_x_discrete(limits=rev) +
    scale_y_continuous(expand=expansion(0,add=c(0,.05))) +
    labs(x="Explanatory variables",       
         y = expression(paste("Pseudo-", R^2)), 
         fill=NULL)  +
    coord_flip()
  return(g1)
}



fn_312_mediation = function(mm_, print.table=FALSE, ...) {
  # labels
  labs2 <- list(age="Age",
                gender='Sex',
                a1_foreign='Nationality',
                a1_education_3classes = "Education",
                a1_household_monthly_income_3classes="Income", 
                center_label="Center",
                
                cluster="Diet cluster",
                
                ex1_fish_freq="Overall fish consumption",
                ex1_tuna='Tuna',
                ex1_salmon='Salmon',
                ex1_pangasius='Pangasius',
                ex1_trout='Trout',
                ex1_shrimp='Shrimp',
                ex1_swiss_fish="Local fish from swiss lakes",
                
                ex1_veggie_vegan='Vegetarian or vegan',
                ex1_tapwater='Tap water (>2L per day)',
                ex1_alcool_weekly='Alcohol: weekly or daily',
                ex1_cig='Tobacco: past or current user',

                ex2_skiwax='Ski wax',
                ex2_spray_imp='Impregnation sprays',
                ex1_hot_meals_package_freq="Disposable food packaging",
                ex1_residence_urban='Residence in urban/industrial area ',
                ex1_carpet="Synthetic or natural carpet at home",
                
                ex1_deo_spray_freq="Deodorant spray",
                ex1_deo_roll_freq="Deodorant roll",
                ex1_bodylotion_freq="Body lotion",
                ex1_handcream_freq="Hand cream",
                ex1_perfume_freq="Perfume",
                ex1_hairspray_freq="Hair spray",
                
                ex2_prof_smoke='Smoke at work',
                ex2_prof_exhaust_gas='Exhaust gas at work',
                ex2_prof_solvent_vapours='Solvent vapours at work',
                ex2_prof_cleaning_agents='Cleaning agents at work',
                ex2_prof_dust='Dust at work',
                ex2_prof_skiwax='Ski wax at work',
                ex2_prof_spray_imp='Impregnation sprays at work'
  )
  corr_cova_labels = tibble(mediator=names(labs2),
                            parameter_lab=unname(unlist(labs2))) %>% 
    dplyr::mutate(parameter_lab=factor(parameter_lab,levels=unname(unlist(labs2))))
  
  treatment_labels = tibble(
    treatment=names(labs2)[1:6],
    treatment_lab=unname(unlist(labs2[1:6]))
  ) %>% 
    mutate(treatment_lab=factor(treatment_lab,levels=c("Age","Sex","Nationality","Education","Income","Center")))
  
  mm_2 = mm_ %>% 
    dplyr::left_join(corr_cova_labels,by = join_by(mediator)) %>% 
    dplyr::left_join(treatment_labels,by = join_by(treatment)) %>%
    dplyr::mutate(away=if_else(proportion_mediated_lob>0,1,0),
                  parameter_lab2=if_else(proportion_mediated_lob>0 ,parameter_lab,"Other")) %>%
    dplyr::mutate(rank=min_rank(-proportion_mediated),
                  parameter_lab3=if_else(rank<=8,parameter_lab,"Other"))
  
  order_med = mm_2 %>%
    arrange(-proportion_mediated) %>% 
    dplyr::pull(parameter_lab2) %>% 
    unique()
  order_med = c(setdiff(order_med, "Other"), "Other")
  
  default_colors = scales::hue_pal()(length(order_med) - 1)
  levels_except_other = setdiff(order_med, "Other")
  named_colors = setNames(default_colors, levels_except_other)
  named_colors["Other"] = "grey80"
  mm_2$parameter_lab2 = factor(mm_2$parameter_lab2,levels=order_med)
  
  g = ggplot(mm_2) +
    geom_col(aes(x=treatment_lab,y=proportion_mediated,fill=parameter_lab2),colour="black") +
    # scale_fill_discrete(limits=order_med)  +
    # geom_hline(yintercept=0,linetype=2, linewidth=1) +
    scale_y_continuous(labels=scales::percent, breaks=c(-.25,0,.25,.5,.75)) +
    scale_fill_manual(values = named_colors, limits=order_med) +
    labs(x="Explanatory variable",y="Proportion mediated",fill="Mediators:") +
    theme(axis.text.x = element_text(angle=45,hjust=1))
    
  if(print.table) { 
  tab1 = mm_2 %>% 
    dplyr::group_by(treatment_lab,proportion_mediated) %>% 
    arrange(treatment_lab,-proportion_mediated) %>% 
    dplyr::mutate(prop_test = qsum(proportion_mediated,proportion_mediated_lob,proportion_mediated_upb,2)) %>% 
    ungroup() %>% 
    dplyr::select(treatment_lab,parameter_lab,prop_test)
  
  print(tab1,n=60)
  }
  return(g)
}