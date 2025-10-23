#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 303 - correlogram
#' author: Julien Riou
#' date: 2024-01-12
#' ---

fn_303_correlogram = function(dat_in, pfas=NULL,  exclusions=NULL, ...) {
  require(GGally)
  dat_in %>% 
    dplyr::select(-did) %>% 
    ggpairs()
  
}