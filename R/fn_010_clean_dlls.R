#' ---
#' title: PFAS exposure and immunity in SHeS-pilot
#' subtitle: 010 - Clean memory and DLLs 
#' author: Julien Riou
#' date: 2023-09-26
#' ---

fn_010_clean_dlls = function() {
  
  # clean memory
  gc()
  
  # clean dlls created by stan compilation (see https://stackoverflow.com/questions/24832030/exceeded-maximum-number-of-dlls-in-r/24834763#24834763)
  dso_filenames <- dir(tempdir(), pattern=.Platform$dynlib.ext)
  filenames  <- dir(tempdir())
  for (i in seq(dso_filenames))
    dyn.unload(file.path(tempdir(), dso_filenames[i]))
  for (i in seq(filenames))
    if (file.exists(file.path(tempdir(), filenames[i])) & nchar(filenames[i]) < 42) # some files w/ long filenames that didn't like to be removeed
      file.remove(file.path(tempdir(), filenames[i]))
  
}