#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 205 - Variable selection with BART
#' author: Julien Riou
#' date: 2024-12-18
#' ---

fn_205_bart_selection = function(dat_in, dependent, covariates, impute, ...) {
  
  # Log-transform, normalize and replace zeros by LOD/sqrt(2)
  if(dependent!="pfas_cluster") {
    y_name = paste0(dependent,"_quant")
    dat_ = dat_in %>% 
      dplyr::rename(y=all_of(y_name)) %>% 
      dplyr::mutate(y=if_else(is.na(y),0.1/sqrt(2),y),
                    y=log(y),
                    y=(y-mean(y))/sd(y)) %>% 
      dplyr::select(y,all_of(covariates))
  } else {
    y_name = dependent
    dat_ = dat_in %>% 
      dplyr::rename(y=all_of(y_name)) %>% 
      dplyr::mutate(y=as.numeric(as.factor(y))) %>% 
      dplyr::select(y,all_of(covariates))
  }
  
  # 1 imputation
  m = 1
  predmat = mice::quickpred(dat_,
                            include = impute, # multiple imputation based on variables listed in `impute`
                            mincor = 1, # exclude all other covariates
                            exclude = "y")
  dat_imp = mice::mice(dat_, 
                       m = m, 
                       predictorMatrix = predmat,
                       print = FALSE)
  dat_imp = mice::complete(dat_imp,1) %>% 
    dplyr::mutate(across(where(is.factor), function(x) as.numeric(as.factor(x)))) %>%
    dplyr::mutate(across(where(is.character), function(x) as.numeric(as.factor(x)))) 
  
  # format
  X_train = dat_imp %>% dplyr::select(-y) %>% as.data.frame()
  Y_train =  dat_ %>% dplyr::pull(y)
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
  
  return(list(var_importance = var_importance, rmse=rmse_within))
  if(FALSE) {
    print(paste("Within-sample RMSE:", round(rmse_within, 2)))
    print(paste("Approximate R2:", round(1-rmse_within^2/var(Y_train), 2)))
    print(var_importance)
    print(var_importance/sum(var_importance))
    
    n_order = names(var_importance)
    tibble(var=names(var_importance),
           imp=var_importance) %>% 
      ggplot() +
      geom_col(aes(x=var,y=imp)) +
      coord_flip() +
      scale_x_discrete(limits=rev(n_order))
  }
  
}


fn_205_bart_selection_long = function(dat_in, exclude="age_group", ...) {
  
  require(BART)
  # Log-transform, normalize and replace zeros by LOD/sqrt(2)
  dat_form = pei %>% 
    tidyr::pivot_longer(cols=starts_with("pfas"),names_to = "pfas") %>% 
    dplyr::relocate(pfas, value, .after="did") %>% 
    dplyr::mutate(pfas=gsub("pfas_","",pfas),
                  pfas=gsub("pfos_total","pfostotal",pfas)) %>% 
    tidyr::separate(pfas,sep = "_",into=c("pfas","detection")) %>% 
    tidyr::pivot_wider(names_from="detection",values_from=value) %>% 
    dplyr::relocate(nonzero, quant, .after="pfas") %>% 
    dplyr::select(-nonzero,-age_group) %>% 
    dplyr::mutate(across(where(is.factor)& !c("did"), function(x) as.numeric(as.factor(x)))) %>% 
    dplyr::mutate(across(where(is.character)& !c("did"), function(x) as.numeric(as.factor(x)))) %>% 
    # normalize and replace NA by LOD/sqrt(2)
    dplyr::group_by(pfas) %>% 
    dplyr::mutate(quant=ifelse(is.na(quant),0.1/sqrt(2),quant),
                  quant=log(quant),
                  quant=(quant-mean(quant))/sd(quant)) %>% 
    dplyr::ungroup()
  
  X_train = dat_form %>% dplyr::select(-did,-quant) %>% as.data.frame()
  Y_train =  dat_form %>% dplyr::pull(quant)
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
  
  return(list(var_importance = var_importance, rmse=rmse_within))
  if(FALSE) {
    print(paste("Within-sample RMSE:", round(rmse_within, 2)))
    print(paste("Approximate R2:", round(1-rmse_within^2/var(Y_train), 2)))
    print(var_importance)
    print(var_importance/sum(var_importance))
    
    n_order = names(var_importance)
    tibble(var=names(var_importance),
           imp=var_importance) %>% 
      ggplot() +
      geom_col(aes(x=var,y=imp)) +
      coord_flip() +
      scale_x_discrete(limits=rev(n_order))
  }
  
}


fn_205_bart_selection_supp = function(dat_in, dependent, covariates, ...) {
  
  # Log-transform, normalize and replace zeros by LOD/sqrt(2)
  y_name = dependent
  dat_ = dat_in %>% 
    dplyr::rename(y=all_of(y_name)) %>% 
    dplyr::mutate(y=as.numeric(as.factor(y))) %>% 
    dplyr::select(y,all_of(covariates))
  
  # format
  X_train = dat_ %>% dplyr::select(-y) %>% as.data.frame()
  Y_train =  dat_ %>% dplyr::pull(y)
  # Fit the BART model
  bart_model = BART::wbart(x.train = X_train, y.train = Y_train)
  # Get variable importance
  var_importance = bart_model$varcount.mean
  var_importance = var_importance[!grepl("pfas",names(var_importance))]
  # Sort and display variable importance
  var_importance = names(sort(var_importance, decreasing = TRUE))
  
  return(var_importance)
}
