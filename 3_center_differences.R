#' ---
#' title: Determinants of PFAS exposure 
#' subtitle: Data analysis 2 - center difference
#' author: Julien Riou
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
#' params:
#'   analysis_date:
#' bibliography: misc/bib.bib  
#' ---

#+ results="hide", warnings="false", echo="false"
analysis_date = "2025-05-23"
controls = readRDS(paste0("savepoints/savepoint_",analysis_date,"/controls.rds"))
source("R/setup.R")
shesp_3 = readRDS(file.path(controls$savepoint,"shesp_3.rds"))

# keep only relevant variables
pei = shesp_3 %>% 
  dplyr::select(did,
                all_of(paste0(controls$pfas_list_retained,"_nonzero")),
                all_of(paste0(controls$pfas_list_retained,"_quant")),
                all_of(controls$covariates_list)) 

#' # Differences across centers
#' 
#' The objective is to understand the differences across centers regarding all exposition variables,
#' and this way try to explain the differences in PFAS levels between participants in Bern and Lausanne.
#'
#' ## Data prepatation
#' 
#' We use one iteration of MICE to create a complete database with only the center and the variables.

pei_center = pei %>% 
  dplyr::select(all_of(controls$covariates_list))
mi = mice::mice(pei_center)
pei_center = mice::complete(mi) %>% 
  as_tibble()


#' 
#' ## Correlation
#' 
#' Positive correlation means that the covariable is higher in Lausanne compared to Bern, while negative correlation means the opposite.

pei_center = pei_center %>%
  dplyr::select(all_of(controls$covariates_list)) %>% 
  dplyr::mutate(across(where(is.factor), function(x) as.numeric(as.factor(x)))) %>%
  dplyr::mutate(across(where(is.character), function(x) as.numeric(as.factor(x)))) %>% 
  dplyr::rename(y=center_label) 
pei_center %>% 
  cor() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Covariable") %>%
  as_tibble() %>%
  filter(Covariable!="y") %>% 
  select(Covariable,y) %>% 
  arrange(y) %>% 
  mutate(rank=min_rank(y)) %>% 
  filter(rank %in% c(1:10,27:36)) %>% 
  print(n=20)

#' ## BART
#' 

# format
X_train = pei_center %>% dplyr::select(-y) %>% as.data.frame()
Y_train =  pei_center %>% dplyr::pull(y)
# Fit the BART model
bart_model = BART::wbart(x.train = X_train, y.train = Y_train)
# BART accuracy
y_pred = predict(bart_model, newdata = X_train)
y_pred_mean = apply(y_pred,2,mean)
rmse_within = sqrt(mean((Y_train - y_pred_mean)^2))
# Get variable importance
var_importance = bart_model$varcount.mean
var_importance = var_importance[!grepl("pfas",names(var_importance))]
# Sort and display variable importance
var_importance = sort(var_importance, decreasing = TRUE)
# Top 10 excluding sociodemo
var_importance_top10 = var_importance[c(1,2,3,5,6,7,8,9,11,12)]
var_importance_top10
save(var_importance,file=file.path(controls$savepoint,"bart_center_label.Rdata"))


#' The importance of each variable is an indicator of the differences across groups. 

#' ## Mediation analysis
#'

pei2 = shesp_3 %>% 
  dplyr::select(all_of(controls$covariates_list),pfas_cluster) %>% 
  dplyr::mutate(y=as.integer(pfas_cluster %in% c("2_Intermediate","3_High")),
                x=center_label,
                .before = age) 


# Mediator models
med_list = list()
med_list[[1]] = f1 = bf(x ~ cluster, family = categorical)
med_list[[2]] = bf(ex1_perfume_freq  ~ x, family = cumulative)
med_list[[3]] = bf(ex1_fish_freq  ~ x, family = cumulative)
med_list[[4]] = bf(ex1_pangasius  ~ x, family=bernoulli)
med_list[[5]] = bf(ex1_salmon ~ x, family=bernoulli)
med_list[[6]] = bf(ex1_deo_roll_freq ~ x, family=cumulative)
med_list[[7]] = bf(ex1_shrimp ~ x, family=bernoulli)
med_list[[8]] = bf(ex1_handcream_freq ~ x, family=cumulative)
med_list[[9]] = bf(ex1_hot_meals_package_freq ~ x, family=cumulative)
med_list[[10]] = bf(ex2_prof_cleaning_agents ~ x, family=bernoulli)

m_y = bf(y ~ x + cluster + ex1_perfume_freq + ex1_fish_freq + ex1_pangasius + ex1_salmon + ex1_deo_roll_freq +
           ex1_shrimp + ex1_handcream_freq + ex1_hot_meals_package_freq + ex2_prof_cleaning_agents, family=bernoulli)


for(i in 1:10) {
  fit = brm(m_y + f1 + set_rescor(FALSE), data = pei2)
  out = bayestestR::mediation(fit)
}

# Fit the joint model
fit = brm(m_z1 + m_z2 + m_z3 + m_z4 + m_z5 + m_z6 + m_z7 + m_z8 + m_z9 + m_z10 + m_y + set_rescor(FALSE), 
          data = pei_center2)
saveRDS(fit,file=file.path(controls$savepoint,"mediation_fit.rds"))

# --- Extract Posterior Draws ---
mediation(fit)


###############################################
# mediation analysis https://easystats.github.io/bayestestR/articles/mediation.html

pacman::p_load(bayestestR,mediation,brms)

# Fit Bayesian mediation model in brms
f1 <- bf(x ~ cluster, family = bernoulli)
f2 <- bf(y ~ x + cluster, family = bernoulli)
m2 <- brm(f1 + f2 + set_rescor(FALSE), data = pei2, refresh = 0)

mediation(m2)
