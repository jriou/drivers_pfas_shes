#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 102 - Histogram for 1 PFAS
#' author: Julien Riou
#' date: 2024-08-13
#' ---


fn_102_hist_pfas = function(dat_in, var, cap="", ...) {
  
  # select and format data
  nonzero = dat_in %>% 
    dplyr::select(all_of(paste0(var,"_nonzero")))
  names(nonzero) = var
  nonzero = nonzero %>% 
    tidyr::pivot_longer(everything())
  
  id_nonzero = nonzero %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(n=sum(value,na.rm=TRUE)) %>%
    dplyr::arrange(-n)
  
  quant = dat_in %>% 
    dplyr::select(all_of(paste0(var,"_quant")))
  names(quant) = var
  
  long_quant = quant %>% 
    dplyr::select(all_of(pull(nonzero,name))) %>% 
    tidyr::pivot_longer(everything()) %>% 
    dplyr::mutate(is_nonzero=pull(nonzero,value)) %>%
    dplyr::mutate(is_nonzero=factor(is_nonzero,labels=c("Non-detectable","Detectable"))) %>% 
    dplyr::filter(!is.na(is_nonzero)) %>% 
    tidyr::replace_na(list(value=0)) %>% 
    dplyr::filter(name %in% id_nonzero$name) 
  
  out = long_quant %>% 
    ggplot() +
    geom_histogram(aes(x=value,fill=is_nonzero), bins=10, colour="black") +
    scale_fill_manual(values=c(col1,col2)) +
    labs(x="Measure (ng/mL)",title=paste0("Figure. ",cap),fill=NULL) +
    theme(legend.position.inside=c(.8,.8))
  
  return(out)
}


fn_102_ridge_pfas = function(dat_in, vars, ...) {
  
  corr_pfas_labels = tibble(name=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels2) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels2))
  
  # select and format data
  nonzero = dat_in %>% 
    dplyr::select(all_of(paste0(vars,"_nonzero")))
  names(nonzero) = vars
  nonzero = nonzero %>% 
    tidyr::pivot_longer(everything())
  
  id_nonzero = nonzero %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(n=sum(value,na.rm=TRUE)) %>%
    dplyr::arrange(-n)
  
  quant = dat_in %>% 
    dplyr::select(all_of(paste0(vars,"_quant")))
  names(quant) = vars
  
  long_quant = quant %>% 
    dplyr::select(all_of(pull(nonzero,name))) %>% 
    tidyr::pivot_longer(everything()) %>% 
    dplyr::mutate(is_nonzero=pull(nonzero,value)) %>%
    dplyr::mutate(is_nonzero=factor(is_nonzero,labels=c("Non-detectable","Detectable"))) %>% 
    dplyr::filter(!is.na(is_nonzero)) %>% 
    tidyr::replace_na(list(value=0)) %>% 
    dplyr::filter(name %in% id_nonzero$name) %>% 
    dplyr::left_join(corr_pfas_labels,by = join_by(name))
  
  pfas_order = long_quant %>% 
    dplyr::group_by(dependent_variable_lab) %>% 
    dplyr::summarise(mean=mean(value)) %>% 
    dplyr::arrange(mean) %>% 
    dplyr::pull(dependent_variable_lab)
  
  controls$pfas_order <<- pfas_order
  
  pfas_detection = long_quant %>% 
    dplyr::group_by(dependent_variable_lab) %>% 
    dplyr::summarise(Detectable=paste0(round(mean(is_nonzero=="Detectable")*100,1),"%"))
  
  out = long_quant %>% 
    filter(is_nonzero=="Detectable") %>% 
    ggplot() +
    ggridges::geom_density_ridges(aes(x=value,y=dependent_variable_lab,fill=dependent_variable_lab),alpha = 0.7, from=0) +
    geom_text(data=pfas_detection,aes(y=dependent_variable_lab,label=Detectable),x=-.2,nudge_x=30,nudge_y=.3,size=1.9) +
    scale_x_continuous(trans="pseudo_log",expand=expansion(c(.1,0)), 
                       breaks=c(0,10,20,40,60), minor_breaks = NULL) +
    scale_y_discrete(limits=pfas_order) +
    scale_fill_viridis_d(limits=pfas_order,guide="none") +
    labs(x="Serum concentration (ng/mL)",y="PFAS",fill=NULL) +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=7))
  
  return(out)
}


fn_102_radar_pfas_cluster = function(dat_in, vars, ...) {
  
  corr_pfas_labels = tibble(name=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels2) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels2))
  
  # select and format data
  quant = dat_in %>% 
    dplyr::select(pfas_cluster,all_of(paste0(vars,"_quant")))
  names(quant) = c("pfas_cluster",vars)
  
  long_quant = quant %>% 
    tidyr::pivot_longer(-1) %>% 
    tidyr::replace_na(list(value=0)) %>%
    # dplyr::filter(!is.na(value)) %>% 
    dplyr::left_join(corr_pfas_labels,by = join_by(name))
  
  pfas_order = long_quant %>% 
    dplyr::group_by(dependent_variable_lab) %>% 
    dplyr::summarise(mean=mean(value)) %>% 
    dplyr::arrange(mean) %>% 
    dplyr::pull(dependent_variable_lab)
  
  pfas_range = long_quant %>% 
    dplyr::group_by(dependent_variable_lab) %>% 
    dplyr::mutate(value=(value-mean(value))/sd(value)) %>%
    dplyr::group_by(pfas_cluster,dependent_variable_lab) %>% 
    dplyr::summarise(median=median(value),bot=quantile(value,0.25),upp=quantile(value,0.75),.groups="drop")
  
  out = pfas_range %>% 
    ggplot(aes(x=dependent_variable_lab,colour=pfas_cluster,group=pfas_cluster)) +
    geom_line(aes(y=median),linewidth=.8, alpha=.8)+ 
    geom_ribbon(aes(ymin=bot,ymax=upp,fill=pfas_cluster),alpha=.2,colour=NA) +
    scale_y_continuous(label=NULL) +
    scale_x_discrete(limits=rev(pfas_order)) +
    scale_colour_manual(values=controls$col_pfas_cluster,labels=controls$lab_pfas_cluster) +
    scale_fill_manual(values=controls$col_pfas_cluster,labels=controls$lab_pfas_cluster) +
    labs(y=NULL,x=NULL,colour="PFAS cluster:",fill="PFAS cluster:") +
    coord_polar(start=-3/14, clip="off") +
    theme_minimal() +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=7),
          legend.position="bottom", 
          legend.title.position="top",
          legend.text = element_text(size=6),
          legend.title = element_text(size=7)) +
    annotate('text', x = 0.5, y = c(-1, 0, 1, 2), label = c('-1','0', '1', '2'), size=2.5, colour="grey") 
  
  return(out)
}
