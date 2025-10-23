#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 104 - ffq heatmap
#' author: Julien Riou
#' date: 2024-11-15
#' ---

fn_104_ffq_heatmap = function(dat_in, cluster_name=NULL, ...) {
  ex1_ffq = paste0(c("ex1_milk","ex1_yog","ex1_cheese",
                     "ex1_egg","ex1_fish","ex1_fruit",
                     "ex1_poultry","ex1_tofu","ex1_meat",
                     "ex1_pork","ex1_starch","ex1_cornflakes",
                     "ex1_whole_grain","ex1_patat","ex1_lentil",
                     "ex1_vegetable","ex1_salad","ex1_walnut",
                     "ex1_sugar","ex1_chips"),
                   "_freq")
  dd = dat_in %>% 
    dplyr::select(any_of(ex1_ffq),cluster) %>% 
    pivot_longer(1:20) %>% 
    dplyr::mutate(value=as.numeric(value)) %>% 
    dplyr::group_by(cluster,name) %>% 
    dplyr::summarise(value=mean(value,na.rm=TRUE),.groups="drop") %>% 
    dplyr::group_by(name) %>% 
    dplyr::mutate(rank=dense_rank(desc(value))) %>% 
    dplyr::mutate(name = recode(name, 
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
                    cluster == "2_dairy_focused" ~ "Dairy-focused (23%)", 
                    cluster == "1_balanced" ~ "Traditional (28%)", 
                    cluster == "3_high_fiber" ~ "Mediterranean (14%)", 
                    cluster == "4_plant_based" ~ "Plant-based (9%)", 
                    cluster == "5_meat_centered" ~ "Western (26%)" 
                  )) %>% 
    dplyr::mutate(col2=ifelse(value>4,"low","high"))
  
  g1 =  ggplot(dd) +
    geom_tile(aes(y=name,x=cluster,fill=value)) +
    geom_text(aes(y=name,x=cluster,label=rank,colour=col2), size=2.5) +
    scale_colour_manual(values=c("black","grey"),guide="none")+
    scale_fill_viridis_c(option = "G",
                         direction = -1,
                         limits = c(1,7),
                         breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
                         labels = c("Rarely/never", "Once a month", "Every week", "1-2 times a week",
                                    "3-6 times a week", "Once a day", "2-3 times a day", "4+ times a day"),
                         name = "Consumption frequency:") +
    scale_y_discrete(limits = rev(c("Fruits", "Cooked vegetables", "Raw vegetables",
                                "Milk", "Yogurt", "Cheese",
                                "Starches", "Breakfast cereals", "Whole grain cereals",
                                "Potatoes products", "Pulses",
                                "Nuts", "Eggs", "Tofu",
                                "Poultry", "Meat (except poultry)", "Transformed meat", "Fish",
                                "Sweets and deserts", "Salted snacks")),
                     expand=expansion(0,0)) +
    scale_x_discrete(limits=c("Traditional (28%)","Dairy-focused (23%)","Western (26%)","Mediterranean (14%)","Plant-based (9%)"),
                     expand=expansion(0,0)) +
    # coord_flip() +
    theme(legend.position="bottom",
          axis.text.x = element_text(angle = 45, hjust = 1,size=6),
          axis.text.y = element_text(size=6),
          legend.text = element_text(angle = 45, hjust = 1,size=6),
          axis.title=element_text(size=7),
          legend.title=element_text(size=7)) +
    labs(x="Diet cluster",y="Food items")  +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        label.position = "bottom",
        ticks = FALSE,
        # Make the key *wide*; adjust if needed for your device size
        barwidth = unit(5, "cm"),
        barheight = unit(.3, "cm")
      )
    ) 
  
  
  # g1 = dat_in %>% 
  #   dplyr::select(any_of(controls$covariates_ffq),any_of(cluster_name)) %>% 
  #   pivot_longer(1:length(controls$covariates_ffq)) %>% 
  #   mutate(value=as.numeric(value)) %>% 
  #   group_by(cluster,name) %>% 
  #   summarise(value=mean(value,na.rm=TRUE),.groups="drop") %>% 
  #   group_by(name) %>% 
  #   mutate(rank=dense_rank(desc(value))) %>% 
  #   ggplot() +
  #   geom_tile(aes(y=name,x=cluster,fill=value)) +
  #   geom_text(aes(y=name,x=cluster,label=rank)) +
  #   scale_y_discrete(expand=expansion(c(0,0))) +
  #   scale_x_discrete(expand=expansion(c(0,0))) +
  #   scale_fill_viridis_c(option = "A")+
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   labs(fill="Mean frequency",x="Cluster",y="Foodstuffs",title="Figure. Clustering (colours show mean frequency, numbers show rankings by row).")
  
  
  return(g1)
}
