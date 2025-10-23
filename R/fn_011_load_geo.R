#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 011 - Load geographic data
#' author: Julien Riou
#' date: 2023-09-26
#' ---

fn_011_load_geo = function() {
  
  # shape file from https://www.swisstopo.admin.ch/fr/repertoire-officiel-des-localites
  shap = sf::read_sf("../data/ortschaftenverzeichnis_plz_2056.shp/AMTOVZ_SHP_LV95/AMTOVZ_ZIP.shp") %>% 
    dplyr::group_by(ZIP4) %>% 
    dplyr::summarise(geometry = st_union(geometry), .groups = "drop")
  shse = readRDS("../data/se.Rds")
  # ZIP code per participant (separated from other data)
  spat = readr::read_csv2("../data/ParticipantManagemen-Did_DATA_2025-05-27_0848.csv") %>% 
    dplyr::select(did, ZIP4=ec_q9_zip) %>% 
    dplyr::mutate(ZIP4=as.character(ZIP4)) 
  
  return(list(did_zip=spat,shapefile=shap, see=shse))
}