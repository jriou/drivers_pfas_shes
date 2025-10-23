#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 300 - figure for regression coefficients adapted for anti-S IgG
#' author: Julien Riou
#' date: 2024-01-22
#' ---

fn_304_fig_regression_antiS = function(coefs, fade_overlap=TRUE, ...) {
  # labels
  corr_pfas_labels = tibble(parameter2=controls$pfas_list,
                            parameter2_lab=controls$pfas_labels) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=unname(unlist(controls$labs))))
  corr_cova_labels = dplyr::bind_rows(
    dplyr::mutate(corr_pfas_labels,
                  parameter2=paste0(parameter2,"_nonzero"),
                  parameter2_lab=paste(parameter2_lab,"(detection)")),
    dplyr::mutate(corr_pfas_labels,
                  parameter2=paste0(parameter2,"_quant"),
                  parameter2_lab=paste(parameter2_lab,"(quantity)"))
  )
  
  # format results
  dd_ = coefs %>% 
    dplyr::filter(grepl("pfas",parameter)) %>% 
    ## labels for dep variable
    dplyr::mutate(dependent_variable_lab="Anti-S IgG") %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    ## check overlap with 1
    dplyr::mutate(away=factor(ifelse(exp_beta_upb<1 | exp_beta_lob>1,"95% CrI away from 1","95% CrI includes 1"),
                              levels=c("95% CrI away from 1","95% CrI includes 1")),
                  parameter2_lab=factor(parameter2_lab,levels=corr_cova_labels$parameter2_lab)) %>% 
    dplyr::filter(!is.na(parameter2_lab))
  
  # baseline plot
  if(fade_overlap) {
    g = ggplot(dd_) +
      geom_pointrange(aes(x=parameter2_lab,y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb, colour=type, alpha=away),
                      position=position_dodge(.5)) +
      scale_alpha_manual(values=c(.2,1)) 
  } else {
    g = ggplot(dd_) +
      geom_pointrange(aes(x=parameter2_lab,y=exp_beta,ymin=exp_beta_lob,ymax=exp_beta_upb,colour=type),
                      position=position_dodge(.5))
  }
  
  # options
  g = g +
    # facet_wrap(~dependent_variable_lab,ncol=5) +
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
