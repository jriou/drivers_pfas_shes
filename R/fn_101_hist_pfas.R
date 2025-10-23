#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 101 - Histograms for PFAS
#' author: Julien Riou
#' date: 2023-09-26
#' ---

fn_101_hist_pfas = function(dat_in, cap="", lower_limit=0, ...) {
  
  # labels
  # pfas_list = controls$pfas_list
  pfas_labels = controls$pfas_labels
  # pfas_list = c("pfas_pfunda","pfas_pfoa")
  pfas_labels = c("PFUnDA","PFOA")
  names(pfas_labels) = pfas_list

  # select and format data
  nonzero = dat_in %>% 
    dplyr::select(all_of(paste0(pfas_list,"_nonzero")))
  names(nonzero) = pfas_labels
  nonzero = nonzero %>% 
    tidyr::pivot_longer(everything())
  
  id_nonzero = nonzero %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(n=sum(value,na.rm=TRUE)) %>%
    dplyr::arrange(-n) %>% 
    dplyr::filter(n>lower_limit)
  
  quant = dat_in %>% 
    dplyr::select(all_of(paste0(pfas_list,"_quant")))
  names(quant) = controls$pfas_labels
  
  long_quant = quant %>% 
    dplyr::select(all_of(pull(nonzero,name))) %>% 
    tidyr::pivot_longer(everything()) %>% 
    dplyr::mutate(is_nonzero=pull(nonzero,value),
                  is_nonzero=factor(is_nonzero,labels=c("Non-detectable","Detectable"))) %>% 
    dplyr::filter(!is.na(is_nonzero)) %>% 
    tidyr::replace_na(list(value=0)) %>% 
    dplyr::filter(name %in% id_nonzero$name) %>% 
    dplyr::mutate(name = factor(name,levels=pfas_labels))
  
  out = long_quant %>% 
    ggplot() +
    geom_histogram(aes(x=value,fill=is_nonzero), bins=30, colour="black") +
    facet_wrap(~name,scales="free_x") +
    scale_fill_manual(values=c(col1,col2)) +
    labs(x="Measure (ng/mL)",title=paste0("Figure. ",cap),fill=NULL) +
    theme(legend.position=c(.9,.1))
  
  return(out)
}
