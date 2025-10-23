#' ---
#' title: PFAS2
#' subtitle: Main 
#' author: jriou
#' date: 2024-06-06
#' output:
#'  html_document:
#'    toc: true
#'    toc_float: true
#'    number_sections: true
#'    code_folding: hide
#'    theme: cosmo
#'    highlight: pygments
#'    fig_width: 10
#'    fig_height: 8
#'    fig_caption: true
#' bibliography: misc/PFAS.bib  
#' ---


# Set-up ------------------------------------------------------------------

controls = list(
  render=TRUE,
  analysis_date = "2025-09-04",
  ffq_grouping = FALSE, # if TRUE, food frequency questions will be grouped in 3 groups (less than weekly, at least weekly, daily and more)
  data_path = file.path("L:/UNISANTE_DESS/S_MC_SHES/PILOT/14-Data & Stats/3- Data/FINAL DATA for analysis/Final DATA"),
  run_datamanagement=TRUE,
  fit_models=TRUE
)
source("R/setup.R")
if(controls$run_datamanagement) {
  controls = fn_000_additional_controls(controls)
}

if(controls$run_datamanagement) {

# Block 1: data preparation -----------------------------------------------

  ## 1.1: load and assemble data
  shesp_1 = fn_001_load()
  
  ## 1.2: apply exclusion criteria
  shesp_2 = fn_002_exclusions(shesp_1)
  
  ## 1.3: create new variables
  shesp_3 = fn_003_new_variables(shesp_2)
  
  ## 1.4: compute sampling weigths
  shesp_3 = fn_004_ponderation(shesp_3)
  controls$survey_design = survey::svydesign(ids = ~1, weights = ~weight, data = shesp_3)
  
  ## 1.5: clustering of food frequency questionnaire
  shesp_3 = fn_005_ffq_clustering(shesp_3)
  
  ## 1.6: clustering of PFASs
  shesp_3 = fn_006_pfas_clustering(shesp_3)
  
# Block 2: description ----------------------------------------------------
  
  ## 2.1: description
  rmarkdown::render("1_description.R", 
                    params=list(analysis_date=controls$analysis_date),
                    output_dir=controls$savepoint)
  
  ## 2.2: selection of variables
  source("1_description.R") # selection added to `controls` from within the description script, need to rerun and save
  write_rds(controls,file=file.path(controls$savepoint,"controls.rds")) # save updated controls
  
} else {
  controls = read_rds(file=file.path(controls$savepoint,"controls.rds")) 
}

# Block 3: inference ------------------------------------------------------

## 3.1: inference
rmarkdown::render("2_determinants_pfas.R",
                  params=list(analysis_date=controls$analysis_date),
                  output_dir=controls$savepoint)

rmarkdown::render("4_determinants_pfas_multilevel.R",
                  params=list(analysis_date=controls$analysis_date),
                  output_dir=controls$savepoint)

rmarkdown::render("8_supplementary.R",
                  output_format = "word_document",
                  output_dir=controls$savepoint)
