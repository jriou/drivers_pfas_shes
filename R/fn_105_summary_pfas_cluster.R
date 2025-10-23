#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 105 - Summary table for PFAS clusters
#' author: jriou
#' date: 2025-05-01
#' ---

fn_105_summary_pfas_cluster = function(dat_in, exclusions = NULL, cap=NULL, ponderation=FALSE, ...) {
  
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
    dplyr::select(pfas_cluster,all_of(paste0(pfas_list,"_nonzero")))
  names(nonzero) = c("pfas_cluster",pfas_list)
  
  quant = dat_in %>% 
    dplyr::select(pfas_cluster,all_of(paste0(pfas_list,"_quant")))
  names(quant) = c("pfas_cluster",pfas_list)
  
  # create tables
  n1 = quant %>% 
    gtsummary::tbl_summary(
      missing="no",
      label=pfas_labels,
      type = list(everything()~"continuous"),
      statistic = list(
        everything() ~ "{p_nonmiss}%; {median} ({p25}–{p75})"
      ),
      by=pfas_cluster) %>% 
    gtsummary::modify_header(list(label ~ "",
                                  stat_1 ~ "**Low**",
                                  stat_2 ~ "**Intermediate**",
                                  stat_3 ~ "**High**")) %>% 
    gtsummary::modify_footnote(all_stat_cols() ~ "Proportion detectable; median (IQR) among detectable")
  
  # merge tables
  out = n1
  
  # caption
  if(!is.null(cap)) {
    cap = paste0("**Table.** ",cap)
    out = out %>% 
      gtsummary::modify_caption(cap)
  }
  return(out)
}
