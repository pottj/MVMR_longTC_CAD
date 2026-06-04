#' ---
#' title: "Summarize results"
#' subtitle: "MVMR - TC on CAD - sex-stratified"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # Load data ####
#' ***
myFiles = list.files(path = "../results/MVMR/",pattern = "MVMR")

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  loaded = load(paste0("../results/MVMR/",myFiles[i]))
  data = get(loaded)
  
  dummy = gsub(".RData","",myFiles[i])
  dummy = unlist(strsplit(dummy,"_"))
  data[,exposure_study := dummy[3]]
  data[,exposure_method := dummy[4]]
  if(length(dummy)==5) data[,setting2 := dummy[5]]
  data
  
}
myTab = rbindlist(dumTab1, fill=T)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
