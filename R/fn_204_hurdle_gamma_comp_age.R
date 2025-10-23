#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 204 - Comparing age implementation
#' author: Julien Riou
#' date: 2024-10-03
#' ---

fn_204_hurdle_gamma_comp_age = function(mod_age,mod_age_group, ...) {
  uni_age = bind_rows(fn_202_hurdle_gamma_extract(mod_age,"age"),fn_202_hurdle_gamma_extract(mod_age_group,"age_group"))
  loo1 = brms::loo(mod_age)
  loo2 = brms::loo(mod_age_group)
  loo_comp = brms::loo_compare(loo1, loo2)
  newdat_age = tidyr::expand_grid(age=20:70,pfas=controls$pfas_list_retained) %>% 
    dplyr::mutate(age_group=cut(age,c(20,30,40,50,60,70),include.lowest=TRUE,right=FALSE)) %>% 
    dplyr::mutate(pfas=gsub("pfas_","",pfas),
                  pfas=gsub("_quant","",pfas),
                  pfas=gsub("pfos_total","pfostotal",pfas)) 
  newdat_age$pred_age = brms::posterior_predict(mod_age,newdata = newdat_age) %>% 
    apply(2,mean)
  newdat_age$pred_age_group = brms::posterior_predict(mod_age_group,newdata = newdat_age) %>% 
    apply(2,mean)
  rr = list(uni_age=uni_age,loo1=loo1,loo2=loo2,loo_comp=loo_comp,pred=newdat_age)
  return(rr)
}

fn_204_plot_comp_age = function(x, ...) {
   gg = ggplot(x$pred) +
      geom_line(aes(x=age,y=pred_age,colour="Continuous")) +
      geom_line(aes(x=age,y=pred_age_group,colour="Discrete")) +
      # geom_jitter(data=filter(pei_long,!grepl("nonzero",pfas)),aes(x=age,y=y,colour=pfas), alpha=.3) +
      scale_y_continuous(trans="pseudo_log") +
      labs(x="Age",y="Mean posterior serum PFAS concentration",colour="Age implementation") +
      facet_wrap(~pfas, scales = "free")
  
  return(gg)
}

fn_204_plot_cont_age = function(x, ...) {
  dd = x$uni_age %>% 
    dplyr::filter(parameter=="age") %>% 
    fn_202_hurdle_gamma_format_labels() %>% 
    mutate(labels=ifelse(pfas=="average_pfas","[Average]",labels))
  pfas_order = dd %>% 
    dplyr::arrange(-Estimate) %>% 
    dplyr::pull(labels)
  # pfas_order = pfas_order[!pfas_order=="[Average]"]
  
  # Show coefficients
  g1 = dd %>% 
    ggplot() +
    geom_hline(yintercept=1,linetype=2) +
    geom_pointrange(aes(x=labels,y=exp(Estimate),ymin=exp(`l-95% CI`),ymax=exp(`u-95% CI`),colour=labels)) +
    geom_pointrange(data=filter(dd,pfas=="average_pfas"),aes(x=labels,y=exp(Estimate),ymin=exp(`l-95% CI`),
                                                             ymax=exp(`u-95% CI`)),colour="grey30") +
    scale_colour_discrete(limits=rev(pfas_order), guide="none") +
    scale_x_discrete(limits=rev(pfas_order)) +
    labs(x="PFAS",y="Fold change per year of age", size=NULL,colour="PFAS") +
    coord_flip()
    
  # Plot showing the increase over age
  g2 = x$uni_age %>% 
    expand_grid(age=20:70) %>% 
    filter(parameter=="age") %>% 
    separate(pfas,",|_",into=c("pfas","1","2")) %>% 
    fn_202_hurdle_gamma_format_labels() %>% 
    dplyr::mutate(pred=exp(Estimate*(age-20)),
                  pred_lb=exp(`l-95% CI`*(age-20)),
                  pred_ub=exp(`u-95% CI`*(age-20))) %>% 
    mutate(mean=ifelse(pfas=="average","Mean","")) %>% 
    ggplot() +
    geom_hline(yintercept=1,linetype=2) +
    geom_line(aes(x=age,y=pred,colour=labels,size=mean,linetype=mean)) +
    geom_ribbon(aes(x=age,ymin=pred_lb,ymax=pred_ub,alpha=mean)) + 
    scale_size_manual(values=c(1,1.6),guide="none") +
    scale_linetype_manual(values=c(1,2),guide="none") +
    scale_alpha_manual(values=c(0,.3),guide="none") +
    scale_colour_discrete(limits=pfas_order) +
    labs(x="Age",y="Predicted serum concentration relative to age 20", size=NULL,colour="      PFAS")
    
  gg = cowplot::plot_grid(g1,g2, rel_widths = c(1,1.5), labels=LETTERS)
  
  return(gg)
}
