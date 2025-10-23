#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 004 - Ponderation
#' author: jriou
#' date: 2024-11-08
#' ---

fn_004_ponderation = function(dat_in) {
  
  # Population data ----
  population_count  = readr::read_csv("data/px-x-0103010000_399_20241108-170837.csv",
                                      locale=locale(encoding="latin1")) %>% 
    dplyr::mutate(center_label=ifelse(Canton=="Vaud","Lausanne","Bern"),
                  gender=ifelse(Sexe=="Femme","Females","Males")) %>%
    dplyr::select(center_label,gender,6:15) %>% 
    tidyr::pivot_longer(3:12) %>% 
    dplyr::mutate(age_group=case_when(name %in% c("20-24 ans","25-29 ans") ~ "[20,30)",
                                      name %in% c("30-34 ans","35-39 ans") ~ "[30,40)",
                                      name %in% c("40-44 ans","45-49 ans") ~ "[40,50)",
                                      name %in% c("50-54 ans","55-59 ans") ~ "[50,60)",
                                      name %in% c("60-64 ans","65-69 ans") ~ "[60,70]")) %>% 
    dplyr::group_by(center_label,gender,age_group) %>% 
    dplyr::summarise(population_count  = sum(value), .groups="drop")
  
  # Counts in SHeS-pilot ----
  sample_counts  = dat_in %>% 
    dplyr::group_by(center_label,gender,age_group) %>% 
    dplyr::summarise(sample_count = n(), .groups="drop")
  
  # Compute weigths ---- 
  weight_data = population_count %>%
    dplyr::left_join(sample_counts, by = c("age_group", "gender", "center_label"))  %>%
    dplyr::mutate(weight = population_count / sample_count)
  
  # Merge back in data ---- 
  dat_out = dat_in %>% 
    dplyr::left_join(weight_data, by = c("age_group", "gender", "center_label"))
  
  write_rds(dat_out,file=file.path(controls$savepoint,"shesp_3.rds"))
  
  return(dat_out)
}

