#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 305 - table with summary of important regression coefficients
#' author: Julien Riou
#' date: 2024-01-23
#' ---

fn_305_table_summary_coefficients = function(coefs_,dat_in, ...) {
  # labels
  corr_pfas_labels = tibble(dependent_variable=controls$pfas_list,
                            dependent_variable_lab=controls$pfas_labels) %>% 
    dplyr::mutate(dependent_variable_lab=factor(dependent_variable_lab,levels=controls$pfas_labels))
  corr_cova_labels = tibble(parameter2=names(controls$labs),
                            parameter2_lab=unname(unlist(controls$labs))) %>% 
    # manual changes
    dplyr::mutate(parameter2_lab=gsub("Age group: ","",parameter2_lab),
                  parameter2_lab=gsub("Education: ","",parameter2_lab),
    )
  corr_cova_labels$parameter2_lab = factor(corr_cova_labels$parameter2_lab,
                                           levels=corr_cova_labels$parameter2_lab)
  if(cov_=="ALL") {
    cov_ = corr_cova_labels$parameter2
  }
  
  
  
  return(list(tab_1,tab_2))
}
