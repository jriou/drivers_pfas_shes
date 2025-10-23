#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 202 - Hurdle gamma regression on long data (all PFASs together)
#' author: Julien Riou
#' date: 2024-10-03
#' ---

fn_202_hurdle_gamma_long = function(dat_in, dependent, covariates, multivariable, adj=NULL, impute=NULL, partial=NULL, savemod=FALSE, ...) {
  
  # format data ----
  y_name = paste0(dependent,"_quant")
  dat_ = dat_in %>% 
    dplyr::select(did, weight,all_of(y_name), all_of(covariates), all_of(adj))
  
  # generate priors ----
  # weakly informative prior normal(0,1) on the log scale
  prior_ = prior("normal(0,1)", class="b") + 
    prior("normal(0,1)", class="b", dpar="hu") +
    prior("normal(0,1)", class="sd") + 
    prior("normal(0,1)", class="sd", dpar="hu") 
  
  # identify adjustment covariates ----
  if(multivariable) {
    if(is.null(adj)) {
      adj = covariates
    }
    adj_label = tolower(paste(unname(unlist(controls$labs[adj])),collapse=", "))
  } else {
    adj_label = ""
  }
  
  # apply multiple imputation ----
  if(!is.null(impute)) {
    m = 5
    predmat = mice::quickpred(dat_,
                              include = impute, # multiple imputation based on variables listed in `impute`
                              mincor = 1, # exclude all other covariates
                              exclude = "y")
    dat_imp = mice::mice(dat_, 
                         m = m, 
                         predictorMatrix = predmat,
                         print = FALSE)
    ## long format for each imputed dataset
    dat_long_imp = list()
    for(i in 1:m) {
      tmp_ = mice::complete(dat_imp,i)
      dat_long_imp[[i]] = tmp_ %>% 
        tidyr::pivot_longer(cols=starts_with("pfas"),names_to = "pfas",values_to="y") %>% 
        dplyr::relocate(pfas, y, .after="did") %>% 
        dplyr::mutate(y=ifelse(is.na(y),0,y)) %>% 
        dplyr::mutate(pfas=gsub("pfas_","",pfas),
                      pfas=gsub("_quant","",pfas),
                      pfas=gsub("pfos_total","pfostotal",pfas)) 
    }
  }
  ## long format
  dat_long = dat_ %>% 
    tidyr::pivot_longer(cols=starts_with("pfas"),names_to = "pfas",values_to="y") %>% 
    dplyr::relocate(pfas, y, .after="did") %>% 
    dplyr::mutate(y=ifelse(is.na(y),0,y)) %>% 
    dplyr::mutate(pfas=gsub("pfas_","",pfas),
                  pfas=gsub("_quant","",pfas),
                  pfas=gsub("pfos_total","pfostotal",pfas)) 
  
  # formula ----
  ## set formula for first model
  if(multivariable) {
    rhs_ = c(covariates[1],adj)
    list_cov = rhs_[!duplicated(rhs_)]
    rhs_ = paste(paste0(list_cov," + (-1+",list_cov,"|pfas)"), collapse=" + ")
    rhs_ = paste0("(1|pfas) + ",rhs_)
  } else {
    rhs_ = paste0(covariates[1]," + (",covariates[1],"|pfas)")
  } 
  focus_cov = covariates[1]
  ## formula format
  form_ = bf(paste0("y ~ ", rhs_), 
             paste0("hu ~ ", rhs_))
  print(form_)
  
  # compile and fit ----
  ## empty res 
  res_ = list()
  if(is.null(impute) |  sum(is.na(dat_[covariates[1]]))==0) {
    mod_ = suppressMessages(
      brms::brm(form_,
                family = hurdle_gamma(),
                data = dat_long, 
                prior = prior_,
                chains = 4,
                cores = 4,
                iter = 2000,
                warmup = 1000,
                control = list(adapt_delta = 0.95)) 
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
  if(savemod) saveRDS(mod_, file=file.path(controls$savepoint,paste0("4_",
                                                                     ifelse(multivariable,"multivariable_","univariable_"),
                                                                     paste(adj,collapse="_"),"_pfas_",covariates[1],".rds")))
  ## extract posterior
  res_[[1]] = fn_202_hurdle_gamma_extract(mod_,focus_cov)
  ## loop on other covariates 
  if(is.null(partial)) partial = 2:length(covariates)
  for(i in partial) {
    ## set formula
    if(multivariable) {
      rhs_ = c(covariates[i],adj)
      list_cov = rhs_[!duplicated(rhs_)]
      rhs_ = paste(paste0(list_cov," + (-1+",list_cov,"|pfas)"), collapse=" + ")
      rhs_ = paste0("(1|pfas) + ",rhs_)
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
    if(savemod) saveRDS(mod_, file=file.path(controls$savepoint,paste0("4_",
                                                                       ifelse(multivariable,"multivariable_","univariable_"),
                                                                       paste(adj,collapse="_"),"_pfas_",covariates[i],".rds")))
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
  k = ncol(samples_fix_1)
  samples_ran_1 = as_draws_df(mod, variable = paste0("^r_pfas\\[.*",focus,".*\\]$"), regex = TRUE) %>% 
    as_tibble() %>% 
    dplyr::select(-.chain,-.iteration,-.draw)
  p = ncol(samples_ran_1)/k
  if(k==1) {
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
  } else {
    ran_1 = NULL
    for(i in 1:k) {
      p_select = (1:p)+(i-1)*p
      ran_1 = samples_ran_1 %>% 
        dplyr::select(all_of(p_select)) %>% 
        dplyr::mutate(across(everything(), ~ . + dplyr::pull(samples_fix_1[,i]))) %>% 
        tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
        dplyr::group_by(pfas) %>% 
        dplyr::summarise(Estimate = mean(value), 
                         Est.Error=sd(value), 
                         `l-95% CI`=quantile(value, 0.025),
                         `u-95% CI`=quantile(value, 0.975),
                         .groups = "drop") %>% 
        dplyr::mutate(pfas=gsub("r_pfas\\[","",pfas),
                      pfas=gsub(paste0(",",focus,"\\]"),"",pfas)) %>% 
        dplyr::mutate(parameter=focus,type="FC",.before="pfas") %>% 
        dplyr::bind_rows(ran_1)
    }
  }
  ## extract posterior summary for random effects (hurdle part)
  samples_fix_2 = as_draws_df(mod, variable = paste0("^b_hu_",focus,"*"), regex=TRUE) %>% 
    as_tibble() %>% 
    dplyr::select(-.chain,-.iteration,-.draw)
  samples_ran_2 = as_draws_df(mod, variable = paste0("^r_pfas__hu\\[.*",focus,".*\\]$"), regex = TRUE) %>% 
    as_tibble() %>% 
    dplyr::select(-.chain,-.iteration,-.draw)
  if(k==1) {
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
  } else {
    ran_2 = NULL
    for(i in 1:k) {
      p_select = (1:p)+(i-1)*p
      ran_2 = samples_ran_2 %>% 
        dplyr::select(all_of(p_select)) %>% 
        dplyr::mutate(across(everything(), ~ . + dplyr::pull(samples_fix_2[,i]))) %>% 
        tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
        dplyr::group_by(pfas) %>% 
        dplyr::summarise(Estimate =mean(-value), # take minus value because hurdle is parameterized using the probability of zero and not the probability of non-zero
                         Est.Error=sd(value), 
                         `l-95% CI`=quantile(-value, 0.025),
                         `u-95% CI`=quantile(-value, 0.975),
                         .groups = "drop") %>% 
        dplyr::mutate(pfas=gsub("r_pfas__hu\\[","",pfas),
                      pfas=gsub(paste0(",",focus,"\\]"),"",pfas)) %>% 
        dplyr::mutate(parameter=paste0("hu_",focus),type="OR",.before="pfas") %>% 
        dplyr::bind_rows(ran_2)
    }
  }
  ## extract posterior summaries for fixed effects
  fix_1 = samples_fix_1 %>% 
    tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
    dplyr::group_by(pfas) %>% 
    dplyr::summarise(Estimate =mean(value), # take minus value because hurdle is parameterized using the probability of zero and not the probability of non-zero
                     Est.Error=sd(value), 
                     `l-95% CI`=quantile(value, 0.025),
                     `u-95% CI`=quantile(value, 0.975)) %>% 
    dplyr::mutate(parameter=gsub("^b_","",pfas),
                  type="FC",
                  pfas="average_pfas",
                  .before="Estimate") %>% 
    dplyr::mutate(Rhat= brms::rhat(mod,names(samples_fix_1)[1]))
  fix_2 = samples_fix_2 %>% 
    tidyr::pivot_longer(everything(), names_to = "pfas", values_to = "value") %>%
    dplyr::group_by(pfas) %>% 
    dplyr::summarise(Estimate =mean(-value), # take minus value because hurdle is parameterized using the probability of zero and not the probability of non-zero
                     Est.Error=sd(value), 
                     `l-95% CI`=quantile(-value, 0.025),
                     `u-95% CI`=quantile(-value, 0.975)) %>% 
    dplyr::mutate(parameter=gsub("^b_hu_","",pfas),
                  type="OR",
                  pfas="average_pfas",
                  .before="Estimate") %>% 
    dplyr::mutate(Rhat= brms::rhat(mod,names(samples_fix_2)[1]))
  ## put together
  rr = dplyr::bind_rows(fix_1,ran_1,fix_2,ran_2) %>% 
    dplyr::mutate(temp1=str_extract(pfas,pattern="(?<=,)(.*?)(?=]$)"),
                  pfas=str_extract(pfas,pattern="^[^,]+"),
                  parameter=ifelse(pfas=="average_pfas",parameter,temp1)) %>% 
    dplyr::select(-temp1)
  return(rr)
}

fn_202_hurdle_gamma_format_labels = function(tt) {
  nn = tt %>%
    dplyr::select(pfas)
  nn2 = tibble(pfas=controls$pfas_list_retained,
               labels=controls$pfas_labels_retained) %>% 
    dplyr::mutate(pfas=gsub("pfas_","",pfas),
                  pfas=gsub("_quant","",pfas),
                  pfas=gsub("pfos_total","pfostotal",pfas)) %>% 
    dplyr::left_join(nn,by = join_by(pfas))
  tt = tt %>% 
    left_join(nn2, by=join_by("pfas"=="pfas"),relationship = "many-to-many")
  return(tt)
}


fn_202_hurdle_gamma_loopread = function(var_string) {
  all_files = list.files(
    path = controls$savepoint, 
    pattern = "^4_multivariable_age_pfas", 
    full.names = TRUE)
  out_ = NULL
  for(i in 1:length(all_files)) {
    tmp = read_rds(all_files[7])
  }
  return(tt)
}
