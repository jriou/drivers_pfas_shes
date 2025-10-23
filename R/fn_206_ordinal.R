#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 206 - Ordinal regression on PFAS clusters
#' author: Julien Riou
#' date: 2025-05-01
#' ---

fn_206_ordinal = function(dat_in, covariates, multivariable, regularization=FALSE, adj=NULL, multiple_runs=TRUE, impute=NULL, monotonic=NULL, ...) {
  
  # format data ----
  dat_ = dat_in %>% 
    dplyr::rename(y=pfas_cluster) %>% 
    dplyr::mutate(y=factor(y,ordered=TRUE)) %>% 
    dplyr::select(y,all_of(covariates),all_of(impute),all_of(adj))
  
  # generate priors ----
  # weakly informative prior normal(0,1) on the log scale
  prior_ = prior("normal(0,1)", class="b")
  # horseshoe prior for regularization (see Piironen and Vehtari 2017, we expect 20% of parameters away from 0)
  if(multivariable & regularization) {
    prior_ = prior(horseshoe(df = 1, par_ratio = 0.2), class="b") +
      prior(horseshoe(df = 1, scale_global = 0.2), class="b", dpar="hu")
  }
  
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
    predmat = mice::quickpred(dat_,
                              include = impute, # multiple imputation based on variables listed in `impute`
                              mincor = 1, # exclude all other covariates
                              exclude = "y")
    dat_imp = mice::mice(dat_, 
                         m = 5, 
                         predictorMatrix = predmat,
                         print = FALSE)
  }
  
  # formula ----
  ## set formula for first model
  if(multivariable) {
    rhs_ = c(covariates[1],adj)
    rhs_ = paste(rhs_[!duplicated(rhs_)], collapse=" + ")
  } else {
    rhs_ = covariates[1]
  } 
  ## implement monotonic effects as in https://cran.r-project.org/web/packages/brms/vignettes/brms_monotonic.html
  if(!is.null(monotonic)) {
    for(i in monotonic) {
      rhs_ = gsub(i,paste0("mo(",i,")"),rhs_)
    }
  }
  ## formula format
  form_ = bf(paste0("y ~ ", rhs_))
  print(form_)
  
  # compile and fit ----
  ## empty res 
  res_ = list()
  if(is.null(impute) |  sum(is.na(dat_[covariates[1]]))==0) {
    ### automatically removes incomplete observations
    mod_ = suppressMessages(
      brms::brm(form_,
                family = cumulative(link = "logit"),
                data = dat_, 
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
                         family = cumulative(link = "logit"),
                         data = dat_imp, 
                         prior = prior_,
                         chains = 4,
                         cores = 4,
                         iter = 2000,
                         warmup = 1000,
                         control = list(adapt_delta = 0.95)) 
    )
  }
  res_[[1]] = summary(mod_)$fixed %>% 
    tibble::rownames_to_column(var="parameter") %>% 
    tibble::as_tibble()
  if(multiple_runs) res_[[1]] = res_[[1]] %>% dplyr::filter(grepl(covariates[1],parameter))
  
  ## loop on other covariates 
  if(multiple_runs) {
    for(i in 2:length(covariates)) {
      ## set formula
      if(multivariable) {
        rhs_ = c(covariates[i],adj)
        rhs_ = paste(rhs_[!duplicated(rhs_)], collapse=" + ")
      } else {
        rhs_ = covariates[i]
      } 
      ## implement monotonic effects as in https://cran.r-project.org/web/packages/brms/vignettes/brms_monotonic.html
      if(!is.null(monotonic)) {
        for(j in monotonic) {
          rhs_ = gsub(j,paste0("mo(",j,")"),rhs_)
        }
      }
      ## formula format
      form_ = bf(paste0("y ~ ", rhs_))
      print(form_)
      ## fit (function update avoids recompiling)
      if(is.null(impute) |  sum(is.na(dat_[covariates[1]]))==0) {
        mod_ = suppressMessages(
          update(mod_,
                 formula. = form_,
                 newdata = dat_)
        )
      } else {
        mod_ = suppressMessages(
          update(mod_,
                 formula. = form_,
                 newdata = dat_imp)
        )
      }
      ## extract results
      res_[[i]] = summary(mod_)$fixed %>% 
        tibble::rownames_to_column(var="parameter") %>% 
        tibble::as_tibble()
      if(multiple_runs) res_[[i]] = res_[[i]] %>% dplyr::filter(grepl(covariates[i],parameter))
      
      if(!is.null(monotonic)) {
        brms::conditional_effects(mod_,variable = "simo", regex = TRUE)
      }
    }
  }
  
  # bind and format results ----
  pars_ = res_ %>% 
    do.call("rbind",.) %>% 
    dplyr::mutate(type=if_else(grepl("hu_",parameter),"Odds ratio of detection","Fold-change"),
                  exp_beta=if_else(grepl("hu_",parameter),exp(-Estimate),exp(Estimate)),
                  exp_beta_lob=if_else(grepl("hu_",parameter),exp(-`u-95% CI`),exp(`l-95% CI`)),
                  exp_beta_upb=if_else(grepl("hu_",parameter),exp(-`l-95% CI`),exp(`u-95% CI`)),
                  parameter2=gsub("hu_","",parameter),
                  # dependent_variable=dependent,
                  adjusted=ifelse(multivariable,"Adjusted","Unadjusted"), 
                  adjustment_variables=adj_label,
                  prior=prior_[1,1])
  
  return(pars_)
}
