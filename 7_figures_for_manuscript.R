#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: Figures for manuscript
#' author: Julien Riou
#' date: 2025-05-15
#' ---


# Setup ----
analysis_date = "2025-09-04"
controls = readRDS(paste0("savepoints/savepoint_",analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds"))
pei = shesp_3 %>% 
  dplyr::select(did,
                pfas_cluster,
                all_of(paste0(controls$pfas_list_retained,"_nonzero")),
                all_of(paste0(controls$pfas_list_retained,"_quant")),
                all_of(controls$covariates_list)) 

controls$lab_pfas_cluster <- c( "#1 (low, 30%)", "#2 (intermediate, 37%)","#3 (high, 33%)")
controls$col_pfas_cluster <- c( "#46B8DA", "#5CB85C","#EEA236")

theme_update(axis.text = element_text(size = 6),
             axis.title = element_text(size=7),
             legend.text = element_text(size = 6),
             legend.title = element_text(size=7),
             strip.text = element_text(size=7))


# TABLE 1: key variables ---- 

tab1 = shesp_3 %>%
  dplyr::select(age,
                gender,
                a1_foreign,
                a1_education_3classes,
                a1_household_monthly_income_3classes,
                center_label
  ) %>%
  frq_table() 
tab1

tab1 %>%
  as_kable_extra(format = "latex",
                 booktabs = TRUE,   
                 linesep = ""
  ) %>%
  cat(file = file.path(controls$savepoint,"tab1.tex"))


# FIGURE 1: PFAS ----

controls$pfas_labels2 = str_squish(str_remove(controls$pfas_labels, "\\s*\\[[^]]*\\]"))

g1 = fn_102_ridge_pfas(pei,var=controls$pfas_list_retained)
g2 = fn_102_radar_pfas_cluster(pei,vars=controls$pfas_list_retained)

cowplot::plot_grid(g1,g2, labels=c("A","B"), ncol=2, rel_widths = c(1,1),label_size=10)
ggsave(file.path(controls$savepoint,"fig1.pdf"),height=7, width=17.8, units="cm", dpi=500)

# Table 2: PFAS

tab2 = fn_100_summary_pfas(shesp_3,
                           exclusions=controls$pfas_allzero)

tab2

tab2 %>%
  as_kable_extra(format = "latex",
                 booktabs = TRUE,   
                 linesep = ""
  ) %>%
  cat(file = file.path(controls$savepoint,"tab2.tex"))


# FIGURE 2: FFQ cluster ---- 

fn_104_ffq_heatmap(shesp_3)
ggsave(file.path(controls$savepoint,"fig2.pdf"),height=14, width=8.7, units="cm", dpi=500)

# FIGURE 3: associations with PFAS clusters ---- 
controls$lab_pfas_cluster <- c( "#1 (low)", "#2 (intermediate)","#3 (high)")

g1 = fn_311_sociodemo_cluster(pei)

uni_ = readRDS(file.path(controls$savepoint,"4_univariable_pfas_cluster.rds")) %>% filter(!grepl("age_group",parameter))
mul_ = readRDS(file.path(controls$savepoint,"4_multivariable_pfas_cluster.rds"))  %>% filter(!grepl("age_group",parameter))
pmul_ = readRDS(file.path(controls$savepoint,"4_pen_multivariable_pfas_cluster.rds")) %>% filter(!grepl("age_group",parameter))
mul2_ = readRDS(file.path(controls$savepoint,"4_multivariable2_pfas_cluster.rds")) %>% filter(!grepl("age_group",parameter))
bmul_ = readRDS(file.path(controls$savepoint,"4_bart_multivariable_pfas_cluster.rds")) %>% filter(!grepl("age_group",parameter))

pfas_order = controls$pfas_order
bart_selected = bart_stars=bmul_ %>% pull(parameter) %>% unique()
g2 = fn_310_forestplot_regression_cluster(mul2_, bart_stars = bart_selected,uplim = 20)

cowplot::plot_grid(g1,g2, labels=c("A","B"), rel_widths = c(1,1.4), label_size=10)
ggsave(file.path(controls$savepoint,"fig3.pdf"),height=18, width=17.8, units="cm", dpi=500)

# FIGURE 4: associations with individual PFAS ---- 

uni_ = readRDS(file.path(controls$savepoint,"2_univariable_pfas.rds")) %>% filter(!grepl("age_group|breastfeed",parameter))
mul2_ = readRDS(file.path(controls$savepoint,"2_multivariable2_pfas.rds")) %>% filter(!grepl("age_group|breastfeed",parameter))
bmul_ = readRDS(file.path(controls$savepoint,"2_bart_multivariable_pfas.rds")) %>% filter(!grepl("age_group|breastfeed",parameter))

fn_306_heatmap_regression3(mul2_,uplim=5.1,botlim=.1,bart=bmul_,horseshoe=pmul_)

ggsave(file.path(controls$savepoint,"fig4.pdf"),height=18, width=17.8, units="cm", dpi=500)


# FIGURE 5: further explorations

## Rsquared
l = load(file.path(controls$savepoint,"6_mod_rsquared.Rdata"))
g1 = fn_312_rsquared(r2_out)

med_all = readRDS(file.path(controls$savepoint,"6_mediation.rds"))
g2 =  fn_312_mediation(med_all)
fn_312_mediation(med_all, print.table = TRUE)
cowplot::plot_grid(g1,g2, labels=c("A","B"), ncol=2, rel_widths = c(1,1.4),label_size=10)
ggsave(file.path(controls$savepoint,"fig5.pdf"),height=8, width=17.8, units="cm", dpi=500)

## MAP
load(file.path(controls$savepoint,"pred_m3.Rdata"))

centers = spat$shapefile  %>% 
  filter(ZIP4 %in% c("1010","3007")) %>% 
  sf::st_centroid()

pred_m3 %>% 
  filter(ZIP4 %in% zip$ZIP4) %>% 
  ggplot() +
  geom_sf(data=zip_map,colour="grey", fill="white") +
  geom_sf(aes(fill=exp(mean))) +
  geom_sf(data=spat$see,fill="grey") +
  geom_sf(data=centers,colour="black",size=1.7,shape=15) +
  scale_fill_viridis_c(trans="log", breaks=c(.5,1,2), limits=c(1/2.4,2.4), option="D",direction=1) +
  coord_sf(xlim=bbox[c(1,3)],ylim=bbox[c(2,4)]) +
  labs(x="Longitude",y="Latitude", fill="Odds-ratio:")

# ggsave(file.path(controls$savepoint,"fig5.png"),height=6, width=9, units="in", dpi=500)
