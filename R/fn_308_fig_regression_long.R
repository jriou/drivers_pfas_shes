#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 308 - figure for regression coefficients from long regression
#' author: Julien Riou
#' date: 2024-12-05
#' ---

fn_308_fig_regression_long = function(coefs, fade_overlap=TRUE, ...) {
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels)) %>% 
    dplyr::mutate(dependent_variable =gsub("pfas_","",dependent_variable ),
                  dependent_variable =gsub("_quant","",dependent_variable ),
                  dependent_variable =gsub("pfos_total","pfostotal",dependent_variable )) 
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=unname(unlist(controls$labs))))
  
  # format results
  dd_ = coefs %>% 
    ## format parameter
    dplyr::mutate(parameter2=gsub("hu_","",parameter)) %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=ifelse(exp_beta_upb<1 | exp_beta_lob>1,"95% CrI away from 1","95% CrI includes 1")) %>% 
    ## filters
    dplyr::filter(grepl("average",pfas)) 
  
  
  # ggplot(dd_) +
  #   geom_pointrange(aes(x=parameter2_lab,y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb, colour=type, alpha=away),
  #                   position=position_dodge(.5)) +
  #   # facet_wrap(~type) %>%
  #   scale_alpha_manual(values=c(1,0.2)) +
  #   scale_y_log10(breaks=c(.05,.1,.2,.5,1,2,5,10,20),minor_breaks=NULL)+
  #   scale_x_discrete(limits=rev) +
  #   coord_flip(ylim=c(1/6,6)) +
  #   geom_hline(yintercept = 1, linetype=2) +
  #   theme(axis.text.x=element_text(angle=45,hjust=1)) +
  #   # labs(x=NULL,colour="Type of association:",alpha="95% credible interval:",y=expression(exp(beta))) +
  #   theme(legend.position="bottom",
  #         legend.background=element_blank())
  
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
    facet_wrap(~type,nrow=1) +
    scale_y_log10(breaks=c(.05,.1,.2,.5,1,2,5,10,20),minor_breaks=NULL)+
    scale_x_discrete(limits=rev) +
    coord_flip(ylim=c(1/3,3)) +
    geom_hline(yintercept = 1, linetype=2) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    labs(x=NULL,colour="Type of association:",alpha="95% credible interval:",y=expression(exp(beta))) +
    theme(legend.position="bottom",
          legend.background=element_blank())
    
  return(g)
}
