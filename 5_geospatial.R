#' ---
#' title: Determinants of PFAS exposure 
#' subtitle: Geospatial analysis
#' author: jriou
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     code_folding: hide
#'     theme: cosmo
#'     highlight: pygments
#'     fig_width: 10
#'     fig_height: 8
#'     fig_caption: true
#' bibliography: misc/bib.bib  
#' ---

#+ results="hide", warnings="false", echo="false"
analysis_date = "2025-05-23"
controls = readRDS(paste0("savepoints/savepoint_",analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds"))
spat = fn_011_load_geo()
pei = shesp_3 %>% 
  dplyr::select(did,
                pfas_cluster,
                all_of(paste0(controls$pfas_list_retained,"_nonzero")),
                all_of(paste0(controls$pfas_list_retained,"_quant")),
                all_of(controls$covariates_list))  %>% 
  dplyr::left_join(spat$did_zip,by=join_by("did"))

controls$lab_pfas_cluster <- c( "#1 (low)", "#2 (intermediate)","#3 (high)")
controls$col_pfas_cluster <- c( "#46B8DA", "#5CB85C","#EEA236")

# aggregate by zip
zip = pei %>% 
  group_by(ZIP4, pfas_cluster) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(proportion = n / sum(n)) %>%
  mutate(max=ifelse(proportion==max(proportion),pfas_cluster,NA)) %>% 
  filter(!is.na(max)) %>% 
  ungroup() %>% 
  select(ZIP4,max) %>% 
  mutate(max=factor(max,
                    levels=c("1_Low","2_Intermediate","3_High"),
                    labels=controls$lab_pfas_cluster))

zip_map = spat$shapefile %>% 
  left_join(zip,relationship = "many-to-many")
zip_map_box = zip_map %>% 
  filter(!is.na(max))
bbox = st_bbox(zip_map_box)

ggplot() +
  geom_sf(data=zip_map,colour="black", fill="white") +
  geom_sf(data=zip_map_box,aes(fill=max)) +
  geom_sf(data=spat$see,fill="grey70") +
  scale_fill_manual(values=controls$col_pfas_cluster) +
  coord_sf(xlim=bbox[c(1,3)],ylim=bbox[c(2,4)]) +
  labs(x="Longitude",y="Latitude",fill="Most frequent\nPFAS cluster:")

# inla
library(INLA)
library(spdep)

pei2 = pei %>% 
  dplyr::select(pfas_cluster,
                ZIP4,
                c("age","gender","a1_foreign","a1_education_3classes","a1_household_monthly_income_3classes","center_label")) %>% 
  dplyr::mutate(y=as.integer(pfas_cluster %in% c("2_Intermediate","3_High")),
                .before = ZIP4)

# build neighbourhood matrix
zip_map =  spat$shapefile %>% 
  left_join(pei2,by=join_by(ZIP4)) %>% 
  mutate(ZIP_code=as.integer(factor(ZIP4)))
nb = spdep::poly2nb(zip_map)
adj = nb2mat(nb, style = "B", zero.policy = TRUE)
inla_zip_graph = tempfile()
nb2INLA(inla_zip_graph, nb)

zip_map_df = zip_map %>% 
  as_tibble() 
m1 = inla(
  y ~ f(ZIP_code, model = "bym2", graph = inla_zip_graph),
  family = "binomial",
  data = zip_map_df,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(m1)
tmp = m1$summary.random$ZIP_code %>%
  dplyr::select(ZIP_code = ID, mean, `0.025quant`, `0.975quant`)
zip_map %>% 
  dplyr::left_join(tmp) %>% 
  ggplot() +
  geom_sf(data=zip_map,colour="black", fill="white") +
  geom_sf(aes(fill=mean)) +
  geom_sf(data=spat$see,fill="grey70") +
  scale_fill_viridis_c() +
  coord_sf(xlim=bbox[c(1,3)],ylim=bbox[c(2,4)]) +
  labs(x="Longitude",y="Latitude", fill="Spatial\nvariation:")
ggsave(file.path(controls$savepoint,"map_m1.png"),height=5, width=7, units="in", dpi=500)


m2 = inla(
  y ~ f(ZIP_code, model = "bym2", graph = inla_zip_graph) +
    age +
    gender +
    a1_foreign +
    a1_education_3classes +
    a1_household_monthly_income_3classes +
    center_label
    ,
  family = "binomial",
  data = zip_map_df,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(m2)
tmp = m2$summary.random$ZIP_code %>%
  dplyr::select(ZIP_code = ID, mean, `0.025quant`, `0.975quant`)
zip_map %>% 
  dplyr::left_join(tmp) %>% 
  ggplot() +
  geom_sf(data=zip_map,colour="black", fill="white") +
  geom_sf(aes(fill=mean)) +
  geom_sf(data=spat$see,fill="grey70") +
  scale_fill_viridis_c() +
  coord_sf(xlim=bbox[c(1,3)],ylim=bbox[c(2,4)]) +
  labs(x="Longitude",y="Latitude", fill="Spatial\nvariation:")
ggsave(file.path(controls$savepoint,"map_m2.png"),height=5, width=7, units="in", dpi=500)


m3 = inla(
  y ~ f(ZIP_code, model = "bym2", 
        graph = inla_zip_graph, 
        hyper = list(
          prec = list(prior = "pc.prec", param = c(1, 0.01)),
          phi = list(prior = "pc", param = c(0.5, 2/3))
        )) +
    age +
    gender +
    a1_foreign +
    a1_education_3classes +
    a1_household_monthly_income_3classes
  ,
  family = "binomial",
  data = zip_map_df,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(m3)
tmp = m3$summary.random$ZIP_code %>%
  dplyr::select(ZIP_code = ID, mean, `0.025quant`, `0.975quant`)
pred_m3 = zip_map %>% 
  dplyr::left_join(tmp)
save(zip,pred_m3,zip_map,spat,bbox,file=file.path(controls$savepoint,"pred_m3.Rdata"))
ggplot(pred_m3) +
  geom_sf(data=zip_map,colour="black", fill="white") +
  geom_sf(aes(fill=mean)) +
  geom_sf(data=spat$see,fill="grey70") +
  scale_fill_viridis_c() +
  coord_sf(xlim=bbox[c(1,3)],ylim=bbox[c(2,4)]) +
  labs(x="Longitude",y="Latitude", fill="Log(odds-ratio):")
ggsave(file.path(controls$savepoint,"map_m3.png"),height=5, width=7, units="in", dpi=500)

