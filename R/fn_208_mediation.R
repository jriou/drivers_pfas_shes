#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 208 - Mediation analysis
#' author: Julien Riou
#' date: 2025-06-16
#' ---

fn_208_mediation = function(dat_in, ovar_, sel_, adj=NULL) {
  out_ = NULL
  for(svar_ in sel_[1:10]) {
    for1_ = as.formula(paste(svar_,ovar_,sep=" ~ "))
    print(for1_)
    nlevels_ = length(unique(dat_imp2[,svar_][[1]]))
    if(nlevels_==2) {
      mediator_model = glm(for1_,family = binomial(link = "logit"),data = dat_in)
    } else if(nlevels_>2) {
      mediator_model = glm(for1_,family = gaussian(link = "identity"),data = dat_in)
    }
    if(is.null(adj)) {
      for2_ = as.formula(paste("y ~",svar_,"+",ovar_))
      
    } else {
      for2_ = as.formula(paste("y ~",svar_,"+",ovar_,"+",paste(adj,collapse="+")))
    }
    print(for2_)
    outcome_model = glm(for2_, family = binomial(link = "logit"), data = dat_in)
    med_ = mediation::mediate(mediator_model, 
                              outcome_model,
                              treat = ovar_, 
                              mediator = svar_,
                              covariates = varsocd.)
    medsum_ = summary(med_) 
    out_ = bind_rows(out_,
                     tibble(treatment=ovar_,
                            mediator=svar_,
                            adj=adj,
                            ACME=medsum_$d.avg,
                            ACME_lob=medsum_$d.avg.ci[1],
                            ACME_upb=medsum_$d.avg.ci[2],
                            ADE= medsum_$z.avg,
                            ADE_lob=medsum_$z.avg.ci[1],
                            ADE_upb=medsum_$z.avg.ci[2],
                            proportion_mediated = medsum_$n.avg,
                            proportion_mediated_lob =  medsum_$n.avg.ci[1],
                            proportion_mediated_upb =  medsum_$n.avg.ci[2]))
  }
  return(out_)
}
