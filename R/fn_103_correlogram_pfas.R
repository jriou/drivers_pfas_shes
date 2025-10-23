#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 103 - correlogram
#' author: Julien Riou
#' date: 2024-11-15
#' ---

fn_103_correlogram_pfas = function(dat_in, exclusions=NULL, lower_limit=5, ...) {

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
  cor_nonzero = suppressWarnings(cor(nonzero,use="pairwise.complete.obs"))
  
  quant = dat_in %>% 
    dplyr::select(all_of(paste0(pfas_list,"_quant")))
  names(quant) = pfas_list
  cor_quant = suppressWarnings(cor(quant,use="pairwise.complete.obs"))
  
  # count non-missing observations
  count_non_missing_pairs <- function(x, y) sum(!is.na(x) & !is.na(y))
  nonmissing_quant = outer(pfas_list, pfas_list, Vectorize(function(x, y) count_non_missing_pairs(quant[[x]], quant[[y]])))
  nonmissing_lim = nonmissing_quant>lower_limit
  nonmissing_lim[nonmissing_lim==0] = NA
  cor_quant = cor_quant * nonmissing_lim
  
  # create the heatmap
  g1 = cor_nonzero %>% 
    reshape2::melt() %>% 
    ggplot(aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Correlation") +
    geom_text(aes(label=round(value,2)),size=3) +
    scale_x_discrete(labels=pfas_labels) +
    scale_y_discrete(labels=pfas_labels,limits=rev) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "", y = "", title = "Figure. Correlation matrix of PFAS detection.")
  
  g2 = cor_quant %>% 
    reshape2::melt() %>% 
    filter(!is.na(value)) %>% 
    ggplot(aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Correlation") +
    geom_text(aes(label=round(value,2)),size=3) +
    scale_x_discrete(labels=pfas_labels) +
    scale_y_discrete(labels=pfas_labels,limits=rev) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "", y = "", title = "Figure. Correlation matrix of PFAS quantities (minimum 5 observations by pair).")
  
  return(g2)
}