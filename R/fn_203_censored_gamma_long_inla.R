#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 202 - Hurdle gamma regression on long data (all PFASs together)
#' author: Julien Riou
#' date: 2024-10-03
#' ---

fn_203_censored_gamma_long_inla = function(dat_in, dependent, censor=0.1, covariates, multivariable, adj=NULL, partial=NULL, savemod=FALSE, ...) {
  
  # format data ----
  y_name = paste0(dependent,"_quant")
  dat_ = dat_in %>% 
    dplyr::select(did, weight,all_of(y_name), all_of(covariates))
  
  # generate priors ----
  # weakly informative prior normal(0,5) on the log scale
  prior_ = prior("normal(0,5)", class="b") + 
    prior("normal(0,5)", class="b", dpar="hu") +
    prior("normal(0,5)", class="sd") + 
    prior("normal(0,5)", class="sd", dpar="hu") 
  
  # identify adjustment covariates ----
  if(multivariable) {
    if(is.null(adj)) {
      adj = covariates
    }
    adj_label = tolower(paste(unname(unlist(controls$labs[adj])),collapse=", "))
  } else {
    adj_label = ""
  }
  
  ## long format
  dat_long = dat_ %>% 
    tidyr::pivot_longer(cols=starts_with("pfas"),names_to = "pfas",values_to="y") %>% 
    dplyr::relocate(pfas, y, .after="did") %>% 
    dplyr::mutate(pfas=gsub("pfas_","",pfas),
                  pfas=gsub("_quant","",pfas),
                  pfas=gsub("pfos_total","pfostotal",pfas),
                  pfas_group=as.numeric(as.factor(pfas)),
                  pfas_group_int=pfas_group,
                  y=ifelse(y<censor,NA,y),
                  dummy=ifelse(is.na(y),1,0))  
  
  # formula ----
  ## set formula for first model
  if(multivariable) {
    rhs_ = c(covariates[1],adj)
    list_cov = rhs_[!duplicated(rhs_)]
    rhs_ = paste(paste0(list_cov,"+ f(pfas_group,",list_cov,",model='iid')"), collapse=" + ")
  } else {
    rhs_ = paste0(covariates[1]," + f(pfas_group,",covariates[1],",model='iid')")
  } 
  focus_cov = covariates[1]
  ## formula format
  form_ = as.formula(paste0("y ~ -1 + f(pfas_group_int,model='iid') + ", rhs_))
  print(form_)
  
  # compile and fit ----
  ## empty res 
  res_ = list()
  if(is.null(impute) |  sum(is.na(dat_[covariates[1]]))==0) {
    mod_ = INLA::inla(y~age,
                      family="gaussian",
                      data=dat_long,
                      control.family = list(dummy = list(type = "left", value = 0.1)),
                      control.compute = list(waic = TRUE)
    )
    mod_ = INLA::inla(y~ -1 + f(pfas_group_int,model='iid') + age + f(pfas_group, age, model = "iid"),
                      family="gamma",
                      data=dat_long,
                      control.family = list(dummy = list(type = "left", value = 0.1)),
                      control.compute = list(waic = TRUE)
    )
    mod_ = INLA::inla(form_,
                      family="gamma",
                      data=dat_long,
                      control.family = list(dummy = list(type = "left", value = 0.1)),
                      control.compute = list(waic = TRUE)
    )
      
  } else {
    ### brm_multiple automatically pools across multiple imputed datasets created with mice https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html
    ### we can't do the imputation during model fitting proposed by Paul Bürkner in the above link as it is not supported for hurdle_gamma models
    mod_ = suppressMessages(
      brms::brm_multiple(form_,
                         family = hurdle_gamma(),
                         data = dat_long_imp, 
                         prior = prior_,
                         chains = 4,
                         cores = 4,
                         iter = 2000,
                         warmup = 1000,
                         control = list(adapt_delta = 0.95)) 
    )
  }
  if(savemod) saveRDS(mod_, file=file.path(controls$savepoint,paste0("4_univariable_pfas_",covariates[1],"_",format(Sys.time(), "%Y%m%d_%H%M%S"),".rds")))
  ## extract posterior
  res_[[1]] = fn_202_hurdle_gamma_extract(mod_,focus_cov)
  ## loop on other covariates 
  if(is.null(partial)) partial = 2:length(covariates)
  for(i in partial) {
    ## set formula
    if(multivariable) {
      rhs_ = c(covariates[i],adj)
      list_cov = rhs_[!duplicated(rhs_)]
      rhs_ = paste(paste0(list_cov," + (",list_cov,"|pfas)"), collapse=" + ")
    } else {
      rhs_ = paste0(covariates[i]," + (",covariates[i],"|pfas)")
    } 
    focus_cov = covariates[i]
    ## formula format
    form_ = bf(paste0("y ~ ", rhs_), 
               paste0("hu ~ ", rhs_))
    print(form_)
    ## fit (function update avoids recompiling)
    if(is.null(impute) |  sum(is.na(dat_[covariates[i]]))==0) {
      mod_ = suppressMessages(
        update(mod_,
               formula. = form_,
               newdata = dat_long)
      )
    } else {
      mod_ = suppressMessages(
        update(mod_,
               formula. = form_,
               newdata = dat_long_imp)
      )
    }
    if(savemod) saveRDS(mod_, file=file.path(controls$savepoint,paste0("4_univariable_pfas_",covariates[i],".rds")))
    ## extract results
    res_[[i]] = fn_202_hurdle_gamma_extract(mod_,focus_cov)
  }
  
  # bind and format results ----
  pars_ = res_ %>% 
    do.call("rbind",.) %>% 
    dplyr::mutate(exp_beta=exp(Estimate),
                  exp_beta_lob=exp(`l-95% CI`),
                  exp_beta_upb=exp(`u-95% CI`),
                  type=ifelse(type=="OR","Odds ratio of detection","Fold-change"),
                  dependent_variable=pfas,
                  adjusted=ifelse(multivariable,"Adjusted","Unadjusted"), 
                  adjustment_variables=adj_label,
                  prior=prior_[1,1])
  
  return(pars_)
}

fn_202_hurdle_gamma_extract = function(mod, focus) {
  ## extract posterior summary for random effects (gamma part)
  samples_fix_1 = as_draws_df(mod, variable = paste0("^b_",focus,"*"), regex=TRUE) %>% 
    as_tibble() %>% 
    dplyr::select(-.chain,-.iteration,-.draw)
  samples_ran_1 = as_draws_df(mod, variable = paste0("^r_pfas\\[.*",focus,".*\\]$"), regex = TRUE) %>% 
    as_tibble() %>% 
    dplyr::select(-.chain,-.iteration,-.draw)
  ran_1 = samples_ran_1 %>% 
    dplyr::mutate(across(everything(), ~ . + dplyr::pull(samples_fix_1))) %>% 
    tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
    dplyr::group_by(pfas) %>% 
    dplyr::summarise(Estimate = mean(value), 
                     Est.Error=sd(value), 
                     `l-95% CI`=quantile(value, 0.025),
                     `u-95% CI`=quantile(value, 0.975),
                     .groups = "drop") %>% 
    dplyr::mutate(pfas=gsub("r_pfas\\[","",pfas),
                  pfas=gsub(paste0(",",focus,"\\]"),"",pfas)) %>% 
    dplyr::mutate(parameter=focus,type="FC",.before="pfas")
  ## extract posterior summary for random effects (hurdle part)
  samples_fix_2 = as_draws_df(mod, variable = paste0("^b_hu_",focus,"*"), regex=TRUE) %>% 
    as_tibble() %>% 
    dplyr::select(-.chain,-.iteration,-.draw)
  samples_ran_2 = as_draws_df(mod, variable = paste0("^r_pfas__hu\\[.*",focus,".*\\]$"), regex = TRUE) %>% 
    as_tibble() %>% 
    dplyr::select(-.chain,-.iteration,-.draw)
  ran_2 = samples_ran_2 %>% 
    dplyr::mutate(across(everything(), ~ . + dplyr::pull(samples_fix_2))) %>% 
    tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
    dplyr::group_by(pfas) %>% 
    dplyr::summarise(Estimate =mean(-value), # take minus value because hurdle is parameterized using the probability of zero and not the probability of non-zero
                     Est.Error=sd(value), 
                     `l-95% CI`=quantile(-value, 0.025),
                     `u-95% CI`=quantile(-value, 0.975),
                     .groups = "drop") %>% 
    dplyr::mutate(pfas=gsub("r_pfas__hu\\[","",pfas),
                  pfas=gsub(paste0(",",focus,"\\]"),"",pfas)) %>% 
    dplyr::mutate(parameter=paste0("hu_",focus),type="OR",.before="pfas")
  ## extract posterior summaries for fixed effects
  fix_1 = samples_fix_1 %>% 
    tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
    dplyr::summarise(Estimate =mean(value), # take minus value because hurdle is parameterized using the probability of zero and not the probability of non-zero
                     Est.Error=sd(value), 
                     `l-95% CI`=quantile(value, 0.025),
                     `u-95% CI`=quantile(value, 0.975)) %>% 
    dplyr::mutate(parameter=focus,type="FC",pfas="average_pfas",.before="Estimate") %>% 
    dplyr::mutate(Rhat= brms::rhat(mod,names(samples_fix_1)[1]))
  fix_2 = samples_fix_2 %>% 
    tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
    dplyr::summarise(Estimate =mean(-value), # take minus value because hurdle is parameterized using the probability of zero and not the probability of non-zero
                     Est.Error=sd(value), 
                     `l-95% CI`=quantile(-value, 0.025),
                     `u-95% CI`=quantile(-value, 0.975)) %>% 
    dplyr::mutate(parameter=paste0("hu_",focus),type="OR",pfas="average_pfas",.before="Estimate") %>% 
    dplyr::mutate(Rhat= brms::rhat(mod,names(samples_fix_2)[1]))
  ## put together
  rr = dplyr::bind_rows(fix_1,ran_1,fix_2,ran_2)
  return(rr)
}

