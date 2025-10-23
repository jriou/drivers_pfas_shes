#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 301 - figure for regression coefficients
#' author: Julien Riou
#' date: 2023-11-21
#' ---

fn_301_fig_regression_by_covariate = function(cov_,coefs_,dat_,fade_overlap=TRUE, ...) {
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=unname(unlist(controls$labs))))
  
  # format results
  dd_ = coefs_ %>% 
    dplyr::filter(grepl(cov_,parameter)) %>% 
    ## order type
    dplyr::mutate(adjusted=factor(adjusted,levels=c("Unadjusted","Adjusted","Adjusted penalized"))) %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2)) %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<1 | exp_beta_lob>1,"95% CrI away from 1","95% CrI includes 1")) 
  # 
  # # table
  # dd_1 = dd_ %>% 
  #   dplyr::filter(type=="Odds ratio of detection")
  # nonzero = dat_ %>% 
  #   dplyr::select(contains(cov_),all_of(paste0(pfas_list,"_nonzero")))
  # nonzero %>% 
  #   gtsummary::tbl_summary(
  #     type = list(everything()~"dichotomous"),
  #     missing="no",
  #     label=pfas_labels,
  #     by=1) 
  # 
  # # nonzero %>% 
  #   tbl_strata("age_group") %>% 
  #   gtsummary::tbl_summary(
  #     type = list(everything()~"dichotomous"),
  #     missing="no",
  #     label=pfas_labels) 
  # 
  # 
  # baseline plot
  if(fade_overlap) {
    g = ggplot(dd_) +
      geom_pointrange(aes(x=dependent_variable_lab,y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb, 
                          colour=adjusted, alpha=away),
                      position=position_dodge(.7)) +
      scale_alpha_manual(values=c(1,0.2),guide="none") 
  } else {
    g = ggplot(dd_) +
      geom_pointrange(aes(x=parameter2_lab,y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb,colour=type),
                      position=position_dodge(.5))
  }
  
  # options
  g = g +
    facet_grid(parameter2_lab~type) +
    scale_y_log10(breaks=c(.05,.1,.2,.5,1,2,5,10,20),minor_breaks=NULL) +
    scale_x_discrete(limits=rev) +
    coord_flip(ylim=c(1/30,30)) +
    geom_hline(yintercept = 1, linetype=2) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    labs(x=NULL,colour="Model:",alpha="95% credible interval:",y=expression(exp(beta))) +
    theme(legend.position="bottom")
    
  return(g)
}
