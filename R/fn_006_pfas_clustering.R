#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: 006 - PFAS clustering
#' author: jriou based on code from morgane gafner
#' date: 2025-04-30
#' ---

fn_006_pfas_clustering = function(dat_in) {
  
  # select PFAS variables
  sel_pfas = paste0(c("pfas_pfba","pfas_pfhpa","pfas_pfoa","pfas_pfna","pfas_pfda","pfas_pfunda","pfas_pfdoda",    
                      "pfas_pfbs","pfas_pfhxs","pfas_pfhps","pfas_lpfos","pfas_brpfos","pfas_pfos_total","pfas_mefosaa"),"_quant")
  dat_1 = dat_in %>% 
    dplyr::select(any_of(sel_pfas)) %>% 
    dplyr::mutate(across(everything(), ~replace_na(.x, 0.01))) %>% 
    dplyr::mutate(across(everything(),~ log(.))) %>% 
    dplyr::mutate(across(everything(),~(. - mean(. , na.rm = TRUE)) / sd(., na.rm = TRUE)))
  
  # pick number of clusters k 
  if(FALSE){
    # euclidian distance
    diss_matrix = cluster::daisy(dat_1,metric = "euclidean")
    best_k_plot_gap_disteuc = factoextra::fviz_nbclust(as.matrix(diss_matrix), pam, method = "gap_stat")
    write_rds(best_k_plot_gap_disteuc,file=file.path(controls$savepoint,"pfas_best_k_plot_gap_disteuc.rds"))
    # gower distance
    diss_matrix = cluster::daisy(dat_1,metric = "gower")
    best_k_plot_gap_distgower = factoextra::fviz_nbclust(as.matrix(diss_matrix), pam, method = "gap_stat")
    write_rds(best_k_plot_gap_distgower,file=file.path(controls$savepoint,"pfas_best_k_plot_gap_distgower.rds"))
    # manhattan distance
    diss_matrix = cluster::daisy(dat_1,metric = "manhattan")
    best_k_plot_gap_distman = factoextra::fviz_nbclust(as.matrix(diss_matrix), pam, method = "gap_stat")
    best_k_plot_gap_distman + labs(title=NULL)
    ggsave(file=file.path(controls$savepoint,"pfas_best_k_plot_gap_distman.pdf"))
    write_rds(best_k_plot_gap_distman,file=file.path(controls$savepoint,"pfas_best_k_plot_gap_distman.rds"))
    
    best_k_plot_gap_disteuc=read_rds(file.path(controls$savepoint,"pfas_best_k_plot_gap_disteuc.rds"))
    best_k_plot_gap_distgower=read_rds(file.path(controls$savepoint,"pfas_best_k_plot_gap_distgower.rds"))
    best_k_plot_gap_distman=read_rds(file.path(controls$savepoint,"pfas_best_k_plot_gap_distman.rds"))
  }
  
  # PAM with Manhattan distance
  k=3
  pam_res = cluster::pam(dat_1, k, metric="manhattan")
  controls$pam_result_pfas <<- pam_res
  dat_out = dat_in %>% 
    dplyr::mutate(pfas_cluster = pam_res$clustering)
  
  # reorder clusters
  if(FALSE) {  
    factoextra::fviz_cluster(pam_res, 
                             # palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
                             geom = "point",
                             ellipse.type = "convex", 
                             ggtheme = theme_bw()
    )
    dd = dat_in %>% 
      dplyr::select(any_of(sel_pfas)) %>% 
      dplyr::mutate(cluster = pam_res$clustering) %>% 
      pivot_longer(1:14) %>% 
      mutate(value=as.numeric(value),
             value=ifelse(is.na(value),0,value)) %>% 
      group_by(cluster,name) %>% 
      summarise(value=mean(value)) %>% 
      group_by(name) %>% 
      mutate(rank=dense_rank(desc(value))) %>% 
      mutate(cluster = case_when(
        cluster == 1 ~ "Intermediate",
        cluster == 2 ~ "High",
        cluster == 3 ~ "Low")) %>%
      mutate(col2=ifelse(value<5,"low","high"))
    
    ggplot(dd) +
      geom_tile(aes(y=name,x=cluster,fill=value)) +
      geom_text(aes(y=name,x=cluster,label=rank,colour=col2)) +
      scale_colour_manual(values=c("grey","black"),guide=FALSE)+
      scale_fill_viridis_c(option = "G" , direction=-1) +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(file=file.path(controls$savepoint,"pfas_cluster_plot.pdf"), width=10,height=6)
  }
  dat_out = dat_out %>% 
    dplyr::mutate(pfas_cluster=case_when(pfas_cluster==1 ~ "2_Intermediate",
                                         pfas_cluster==2 ~ "3_High",
                                         pfas_cluster==3 ~ "1_Low"))
  
  write_rds(dat_out,file=file.path(controls$savepoint,"shesp_3.rds"))
  
  return(dat_out)
}

