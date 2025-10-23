#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 306 - heatmap for regression coefficients
#' author: Julien Riou
#' date: 2023-11-01
#' ---

fn_307_heatmap_regression_selection = function(coefs, pen, selection_limit, type="FC", uplim=NULL, ...) {
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    dplyr::mutate(parameter2_lab=factor(parameter2_lab,levels=unname(unlist(controls$labs))))
  
  # format results
  dd_ = coefs %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many")
  
  # penalized
  pp_ = pen %>% 
    ## labels for PFAS
    dplyr::left_join(corr_pfas_labels,by = join_by(dependent_variable)) %>% 
    ## labels for covariates
    dplyr::left_join(corr_cova_labels,by = join_by(parameter2),relationship = "many-to-many") %>% 
    dplyr::select(parameter, dependent_variable , pen_beta = Estimate)
  
  # select 
 dd_ = dplyr::left_join(dd_,pp_,by = join_by(parameter, dependent_variable)) %>% 
   dplyr::mutate(away=abs(pen_beta)>selection_limit)
  
  if(type=="OR") {
    dd_ =  dd_ %>% 
      filter(type=="Odds ratio of detection")
    if(is.null(uplim)) {
      uplim = round(max(2,dd_$exp_beta),2)
      botlim = round(min(.5,dd_$exp_beta),2)  
    } else {
      botlim = round(1/uplim,2)
    }
    
    g1 = dd_ %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta,colour=factor(away))) +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,colour=factor(away)),fill="transparent") +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
      scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=rev) +
      scale_colour_manual(values=c("black","transparent"),guide="none") +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom",
            background.grid=element_blank()) +
      labs(x=NULL,y="PFAS",fill="Odds ratio") 
    
    g2 = dd_ %>% 
      filter(away==TRUE) %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
      geom_text(aes(x=parameter2_lab,y=dependent_variable_lab,label=sprintf("%.2f",exp_beta))) +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +        scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=rev) +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom",
            panel.grid=element_blank()) +
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
    
    g1 = dd_ %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta,colour=factor(away))) +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,colour=factor(away)),fill="transparent") +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
      scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=rev) +
      scale_colour_manual(values=c("black","transparent"),guide="none") +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            legend.position="bottom") +
      labs(x=NULL,y="PFAS",fill="Fold-change") 
    
    g2 = dd_ %>% 
      filter(away==TRUE) %>% 
      ggplot() +
      geom_tile(aes(x=parameter2_lab,y=dependent_variable_lab,fill=exp_beta)) +
      geom_text(aes(x=parameter2_lab,y=dependent_variable_lab,label=sprintf("%.2f",exp_beta))) +
      # scale_fill_viridis_c(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim)) +
      scale_fill_gradient2(breaks=c(0.05,.5,.1,.2,.5,1,2,5,10,20,30),trans="log",limits=c(botlim,uplim),
                           mid="white",midpoint=0,low="dodgerblue",high="firebrick") +  
      scale_x_discrete(expand=expansion(0,0)) +
      scale_y_discrete(expand=expansion(0,0),limits=rev) +
      theme(axis.text.x=element_text(angle=45,hjust=1)) +
      labs(x=NULL,y="PFAS",fill="Fold-change") 
  }
  return(list(g1,g2))
}
