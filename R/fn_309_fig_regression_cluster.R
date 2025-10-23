#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 309 - figure for regression coefficients on pfas clusters
#' author: Julien Riou
#' date: 2025-05-12
#' ---

fn_309_fig_regression_cluster = function(coefs, fade_overlap=TRUE, ...) {
  # labels
  corr_pfas_labels = tibble(dependent_variable=c("mu2","mu3"),
                            dependent_variable_lab=c("Intermediate PFAS levels","High PFAS levels")) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=c("Intermediate PFAS levels","High PFAS levels")))
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=unname(unlist(controls$labs))))
  
  # format results
  dd_ = coefs %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::mutate(parameter2=gsub("mu2Intermediate_","",parameter2),
                  parameter2=gsub("mu3High_","",parameter2)) %>% 
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<1 | exp_beta_lob>1,"95% CrI away from 1","95% CrI includes 1")) %>% 
    ## remove intercepts
    dplyr::filter(parameter2!="Intercept")
  
  # baseline plot
  if(fade_overlap) {
    g = ggplot(dd_) +
      geom_pointrange(aes(x=parameter2_lab,y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb, colour=type, alpha=away),
                      position=position_dodge(.5)) +
      scale_alpha_manual(values=c(1,0.2)) 
  } else {
    g = ggplot(dd_) +
      geom_pointrange(aes(x=parameter2_lab,y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb,colour=type),
                      position=position_dodge(.5))
  }
  
  # options
  g = g +
    facet_wrap(~dependent_variable_lab,nrow=1) +
    scale_y_log10(breaks=c(.05,.1,.2,.5,1,2,5,10,20),minor_breaks=NULL)+
    scale_x_discrete(limits=rev) +
    coord_flip(ylim=c(1/30,30)) +
    geom_hline(yintercept = 1, linetype=2) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    labs(x=NULL,colour="Type of association:",alpha="95% credible interval:",y=expression(exp(beta))) +
    theme(legend.position="bottom",
          legend.background=element_blank())
    
  return(g)
}
