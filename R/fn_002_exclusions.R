#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 002 - Exclusions
#' author: jriou
#' date: 2024-06-12
#' ---

fn_002_exclusions = function(dat_in) {
  
  # Select dataset
  stage_0 = dat_in
  
  # Stage 1: Keep only random sample 
  stage_1 = stage_0 %>% 
    dplyr::filter(population=="R")
  
  # Stage 2: Remove participants without PFAS measurement
  stage_2 = stage_1 %>% 
    dplyr::filter(pfas_not_done==0)
  
  # Stage 3: Nothing 
  stage_3 = stage_2 # %>% ...
  
  # Stage 4: Nothing 
  stage_4 = stage_3 # %>% ...
  
  # Add stages if necessary and report below...
  
  # Put all stages in a list (for creating the flowchart)
  dat_out = list(stage_0=stage_0,
                 stage_1=stage_1,
                 stage_2=stage_2,
                 stage_3=stage_3,
                 stage_4=stage_4)
  
  # Save all stages
  write_rds(dat_out,file=file.path(controls$savepoint,"shesp_2.rds"))
  
  # Return only last stage 
  
  return(dat_out[[length(dat_out)]])
}