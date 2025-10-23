#' ---
#' title: Unified data-management in SHeS-pilot
#' subtitle: Initialize R session
#' author: Julien Riou
#' date: 2023-09-26
#' ---



# libs --------------------------------------------------------------------

library(pacman)
pacman::p_load(reshape2,
               tidyverse,
               gt,
               gtsummary,
               gWQS,
               brms,
               GGally,
               mice,
               survey,
               cluster,
               BART,
               ggridges,
               factoextra,
               sf)


# options -----------------------------------------------------------------

options(width=120)
options(mc.cores = 4)

# set paths ---------------------------------------------------------------

path_script = 
  file.path("R/")

# source functions --------------------------------------------------------

all_fns = 
  dir("R/", pattern="fn", full.names = TRUE) %>% 
  file.path()
lapply(all_fns, function(x) source(x,echo=FALSE))


# create savepoint --------------------------------------------------------

controls$savepoint = 
  paste0("savepoints/savepoint_",controls$analysis_date) %>% 
  file.path()
dir.create(controls$savepoint, showWarnings = FALSE)


# esthetics ---------------------------------------------------------------

theme_set(theme_bw())
col1 = "firebrick"
col2 = "dodgerblue"
col3 = "chartreuse"
col4 = "goldenrod"

# theme_gtsummary_journal(journal = "jama")


# short custom functions --------------------------------------------------

inv_logit = function(x) exp(x)/(1+exp(x))

frq_table = function(x,by=NULL, cap=NULL, missing="always", labels=TRUE, ...) {
  if(labels) {
    out = suppressWarnings(
      gtsummary::tbl_summary(
        x,
        by=by,
        label=controls$labs,
        type = all_dichotomous() ~ "categorical",
        statistic = list(all_continuous() ~ "{median} ({min}, {max})", all_categorical() ~ "{n} ({p}%)"),
        missing=missing)
    ) } else {
      out = suppressWarnings(
        gtsummary::tbl_summary(
          x,
          by=by,
          # label=controls$labs,
          type = all_dichotomous() ~ "categorical",
          statistic = list(all_continuous() ~ "{median} ({min}, {max})", all_categorical() ~ "{n} ({p}%)"),
          missing=missing)
      )
    }
  if(!is.null(cap)) {
    cap = paste0("**Table.** ",cap)
    out = out %>% 
      gtsummary::modify_caption(cap)
  }
  return(out)
}

qsum = function(a,b,c,digits) {
  out = 
    paste0(
      formatC(a,format="f",big.mark=",",digits=digits),
      " (",
      formatC(b,format="f",big.mark=",",digits=digits),
      " to ",
      formatC(c,format="f",big.mark=",",digits=digits),
      ")")
  return(out)
}

qsumaway = function(a,b,c,d,digits) {
  out = 
    paste0(
      d,
      formatC(a,format="f",big.mark=",",digits=digits),
      " (",
      formatC(b,format="f",big.mark=",",digits=digits),
      " to ",
      formatC(c,format="f",big.mark=",",digits=digits),
      ")",
      d)
  return(out)
}

qsum_p = function(a,b) {
  out = 
    paste0(
      formatC(a,format="f",digits=0),
      " (",
      formatC(b,format="f",digits=0),
      "%)")
  return(out)
}
convert_crlf <- function (infile) { 
  print(infile) 
  txt <- readLines(infile) 
  f <- file(infile, open = "wb") 
  cat(txt, file = f, sep = "\n") 
  close(f) 
}