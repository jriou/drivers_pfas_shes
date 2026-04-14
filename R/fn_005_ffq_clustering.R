#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 005 - FFQ clustering
#' author: jriou based on code from morgane gafner
#' date: 2024-11-15
#' ---

fn_005_ffq_clustering = function(dat_in) {
  
  # select FFQ variables
  ex1_ffq = paste0(c("ex1_milk","ex1_yog","ex1_cheese",
                     "ex1_egg","ex1_fish","ex1_fruit",
                     "ex1_poultry","ex1_tofu","ex1_meat",
                     "ex1_pork","ex1_starch","ex1_cornflakes",
                     "ex1_whole_grain","ex1_patat","ex1_lentil",
                     "ex1_vegetable","ex1_salad","ex1_walnut",
                     "ex1_sugar","ex1_chips"),
                   "_freq")
  controls$covariates_ffq <<- ex1_ffq
  dat_1 = dat_in %>% 
    dplyr::select(any_of(ex1_ffq))
  
  # run one iteration of MICE to deal with missing values
  set.seed(111)
  dat_2 = dat_1 %>% 
    mice::mice(printFlag = FALSE) %>% 
    mice::complete(1)
  
  # pick number of clusters k 
  if(FALSE){
    # euclidian distance
    best_k_plot_gap_disteuc = factoextra::fviz_nbclust(dat_2, pam, method = "gap_stat",metric = "euclidean")
    write_rds(best_k_plot_gap_disteuc,file=file.path(controls$savepoint,"best_k_plot_gap_disteuc.rds"))
    # gower distance
    best_k_plot_gap_distgower = factoextra::fviz_nbclust(dat_2, pam, method = "gap_stat",metric = "gower")
    write_rds(best_k_plot_gap_distgower,file=file.path(controls$savepoint,"best_k_plot_gap_distgower.rds"))
    # manhattan distance
    best_k_plot_gap_distman = factoextra::fviz_nbclust(dat_2, pam, method = "gap_stat",metric = "manhattan")
    best_k_plot_gap_distman + labs(title=NULL)
    ggsave(file=file.path(controls$savepoint,"best_k_plot_gap_distman.pdf"))
    write_rds(best_k_plot_gap_distman,file=file.path(controls$savepoint,"best_k_plot_gap_distman.rds"))
    
    best_k_plot_gap_disteuc=read_rds(file.path(controls$savepoint,"best_k_plot_gap_disteuc.rds"))
    best_k_plot_gap_distgower=read_rds(file.path(controls$savepoint,"best_k_plot_gap_distgower.rds"))
    best_k_plot_gap_distman=read_rds(file.path(controls$savepoint,"best_k_plot_gap_distman.rds"))
  }
  
  # PAM with Manhattan distance
  k = 5
  controls$pam_result = cluster::pam(dat_2,metric = "manhattan", k)
  dat_out = dat_in %>% 
    dplyr::mutate(cluster = controls$pam_result$clustering)
  
  # reorder clusters
  if(FALSE) {  
    dd = dat_2 %>% 
      dplyr::mutate(cluster = controls$pam_result$clustering) %>% 
      pivot_longer(1:20) %>% 
      mutate(value=as.numeric(value)) %>% 
      group_by(cluster,name) %>% 
      summarise(value=mean(value)) %>% 
      group_by(name) %>% 
      mutate(rank=dense_rank(desc(value))) %>% 
      mutate(name = recode(name, 
                           "ex1_milk_freq" = "Milk",
                           "ex1_yog_freq" = "Yogurt",
                           "ex1_cheese_freq" = "Cheese",
                           "ex1_egg_freq" = "Eggs",
                           "ex1_poultry_freq" = "Poultry",
                           "ex1_meat_freq" = "Meat (except poultry)",
                           "ex1_pork_freq" = "Transformed meat",
                           "ex1_tofu_freq" = "Tofu", 
                           "ex1_fish_freq" = "Fish", 
                           "ex1_starch_freq" = "Starches", 
                           "ex1_cornflakes_freq" = "Breakfast cereals", 
                           "ex1_whole_grain_freq" = "Whole grain cereals", 
                           "ex1_patat_freq" = "Potatoes products",
                           "ex1_lentil_freq" = "Pulses", 
                           "ex1_fruit_freq" = "Fruits",
                           "ex1_vegetable_freq" = "Cooked vegetables", 
                           "ex1_salad_freq" = "Raw vegetables", 
                           "ex1_walnut_freq" = "Nuts", 
                           "ex1_sugar_freq" = "Sweets and deserts",
                           "ex1_chips_freq" = "Salted snacks"), 
             cluster = case_when(
               cluster == 1 ~ "Dairy focused", 
               cluster == 2 ~ "Balanced", 
               cluster == 3 ~ "High-fiber", 
               cluster == 4 ~ "Plant-based", 
               cluster == 5 ~ "Meat-centered" 
             )) %>% 
      mutate(col2=ifelse(value<5,"low","high"))
      
    ggplot(dd) +
      geom_tile(aes(y=name,x=cluster,fill=value)) +
      geom_text(aes(y=name,x=cluster,label=rank,colour=col2)) +
      scale_colour_manual(values=c("grey","black"),guide=FALSE)+
      scale_fill_viridis_c(option = "G",
                          direction = -1,
                           breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
                           labels = c("Rarely/never", "Once a month", "Every week", "1-2 times a week",
                                      "3-6 times a week", "Once a day", "2-3 times a day", "4+ times a day"),
                           name = "Consumption \n frequency") +
      scale_y_discrete(limits = c("Fruits", "Cooked vegetables", "Raw vegetables",
                                  "Milk", "Yogurt", "Cheese",
                                  "Starches", "Breakfast cereals", "Whole grain cereals",
                                  "Potatoes products", "Pulses",
                                  "Nuts", "Eggs", "Tofu",
                                  "Poultry", "Meat (except poultry)", "Transformed meat", "Fish",
                                  "Sweets and deserts", "Salted snacks"),
                       expand=expansion(0,0)) +
      scale_x_discrete(limits=rev(c("Balanced","Dairy focused","Meat-centered","High-fiber","Plant-based")),
                       expand=expansion(0,0)) +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(file=file.path(controls$savepoint,"ffq_cluster_plot.pdf"), width=10,height=6)
  }
  dat_out = dat_out %>% 
    dplyr::mutate(cluster=case_when(cluster==2 ~ "1_balanced",
                                    cluster==1 ~ "2_dairy_focused",
                                    cluster==3 ~ "3_high_fiber",
                                    cluster==4 ~ "4_plant_based",
                                    cluster==5 ~ "5_meat_centered"))
  
  write_rds(dat_out,file=file.path(controls$savepoint,"shesp_3.rds"))
  
  return(dat_out)
}

