#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 100 - Summary table for PFAS
#' author: jriou
#' date: 2024-06-12
#' ---

fn_100_summary_pfas = function(dat_in, exclusions = NULL, cap=NULL, ponderation=FALSE, groupby=NULL, ...) {
  
  # labels
  pfas_list = controls$pfas_list
  pfas_labels = as.list(controls$pfas_labels)
  if(!is.null(exclusions)) {
    pfas_list = controls$pfas_list[-exclusions]
    pfas_labels = as.list(controls$pfas_labels[-exclusions])
  }
  names(pfas_labels) = pfas_list
  
  # select and format data
  nonzero = dat_in %>% 
    dplyr::select(all_of(paste0(pfas_list,"_nonzero")))
  names(nonzero) = pfas_list
  
  quant = dat_in %>% 
    dplyr::select(all_of(paste0(pfas_list,"_quant")))
  names(quant) = pfas_list
  
  # order
  pfas_order2 =  quant %>% 
    tidyr::pivot_longer(everything()) %>% 
    tidyr::replace_na(list(value=0)) %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(mean=mean(value)) %>% 
    dplyr::arrange(mean) %>% 
    dplyr::pull(name)
  
  controls$pfas_order2 <<- pfas_order2
  
  # create tables
  n1 = nonzero %>% 
    mutate(across(everything(), ~ factor(.x,levels=c("0","1")))) %>%
    gtsummary::tbl_summary(type = list(everything()~"dichotomous"),
                           statistic = list(all_dichotomous() ~ "{n} ({p}%)"),
                           label=pfas_labels,
                           include=rev(pfas_order2),
                           value = list(everything() ~ "1"),
                           missing="no") %>% 
    modify_header(list(label ~ "**Measure**",
                                  all_stat_cols() ~"**Detected**")) 
  
  n2 = quant %>% 
    gtsummary::tbl_summary(
      missing="no",
      label=pfas_labels,
      include=rev(pfas_order2),
      type = list(everything()~"continuous"),
      statistic = list(everything()~"{median} ({p25} to {p75}); P95={p95}")) %>% 
    gtsummary::modify_header(list(label ~ "**Measure**",
                                  all_stat_cols() ~"**Quantity (among detected)**")) %>% 
    gtsummary::modify_footnote(all_stat_cols() ~ "Median (IQR); 95th percentile")
  
  # ponderation
  if(ponderation) {
    nonzero$weight = dat_in$weight
    n1 = nonzero %>%
      mutate(across(starts_with("pfas"), ~ factor(.x,levels=c("0","1")))) %>%
      survey::svydesign(ids = ~1, weights = ~weight, data = .) %>% 
      gtsummary::tbl_svysummary(type = list(starts_with("pfas")~"dichotomous"),
                                statistic = list(all_dichotomous() ~ "{p}%"),
                                label=pfas_labels,
                                include=rev(pfas_order2),
                                value = list(starts_with("pfas") ~ "1"),
                                missing="no") %>% 
      gtsummary::modify_header(list(label ~ "**Measure**",
                                    all_stat_cols() ~"**Detected**"))  %>% 
      gtsummary::modify_table_body(filter, !(variable == "weight"))
    
    quant$weight = dat_in$weight
    n2 = survey::svydesign(ids = ~1, weights = ~weight, data = quant) %>% 
      gtsummary::tbl_svysummary(
        missing="no",
        label=pfas_labels,
        include=rev(pfas_order2),
        type = list(everything()~"continuous"),
        statistic = list(everything()~"{median} ({p25} to {p75}); P95={p95}")) %>% 
      gtsummary::modify_header(list(label ~ "**Measure**",
                                    all_stat_cols() ~"**Quantity (among detected)**")) %>% 
      gtsummary::modify_footnote(all_stat_cols() ~ "Median (IQR); 95th percentile") %>% 
      gtsummary::modify_table_body(filter, !(variable == "weight"))
    
  }
  
  # ponderation
  if(!is.null(groupby)) {
    nonzero$group = dat_in[[groupby]]
    n1 = nonzero %>%
      mutate(across(starts_with("pfas"), ~ factor(.x,levels=c("0","1")))) %>%
      gtsummary::tbl_summary(type = list(starts_with("pfas")~"dichotomous"),
                                statistic = list(all_dichotomous() ~ "{p}%"),
                                label=pfas_labels,
                                by = group,
                                include=rev(pfas_order2),
                                value = list(starts_with("pfas") ~ "1"),
                                missing="no") %>% 
      gtsummary::modify_header(list(label ~ "**Measure**"))
    
    quant$group = dat_in[[groupby]]
    n2 = quant %>% 
      gtsummary::tbl_summary(
        missing="no",
        label=pfas_labels,
        by = group,
        include=rev(pfas_order2),
        type = list(everything()~"continuous"),
        statistic = list(everything()~"{median} ({p25} to {p75})")) %>% 
      gtsummary::modify_header(list(label ~ "**Measure**")) %>% 
      gtsummary::modify_footnote(all_stat_cols() ~ "Median (IQR)")
    
  }
  
  # merge tables
  out = gtsummary::tbl_merge(list(n1,n2),
                             tab_spanner = c(NA_character_,NA_character_)) 
  
  # caption
  if(!is.null(cap)) {
    cap = paste0("**Table.** ",cap)
    out = out %>% 
      gtsummary::modify_caption(cap)
  }
  
  return(out)
}
