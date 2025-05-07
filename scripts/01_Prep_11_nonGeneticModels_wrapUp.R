#' ---
#' title: "Test non-genetic model - wrap up"
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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
myFiles = list.files(path = "../results/",pattern = "NonGeneticAssociations")
myFiles

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  load(paste0("../results/",myFiles[i]))
  res
}
myTab = rbindlist(dumTab1, use.names = T,fill=T)

#' # Check effect of age ####
#' ***
myTab[estimate == "mu_age" & comment == "OK" & setting == "all", c(1,3:6,8:9)]

#' Effect is negative for men and post-menopausal women, and positive for women and pre-menopausal women(that is expected). Effect size in post-menopausal women is lower than in men (-0.005 vs. -0.018 in RIRI model, -0.010 vs. -0.024 in RI model).  
#' 
myTab[estimate == "mu_age" & comment == "OK" & setting == "subset", c(1,3:6,8:9)]

#' Effect is still negative for men, but there are mixed results for the post-menopausal women: positive in the RIRI model (but just significant), and negative in the RI model (similar effect size and significance as before). The effect is still positive for women and pre-menopausal women. 
#' 
myTab[estimate == "sigma_age" & comment == "OK" & setting == "all", c(1,3:6,8:9)]

#' Effect is always negative, but stronger in men (e.g. -0.011 vs -0.004 in R3 model). The effect in premenopausal women is not significant (-0.002, p=0.092 in RI model).
#' 
myTab[estimate == "sigma_age" & comment == "OK" & setting == "subset", c(1,3:6,8:9)]

#' In the statin-free subset, age still has a negative effect on the TC variability in men and post-menopausal women, but not women or pre-menopausal women. The effect in the pre-menopausal women is significant (0.031, p=2.63e-30 in RI model). 

#' # Check effect of lipid lowering medication ####
#' ***
myTab[estimate == "mu_lipLowMed" & comment == "OK",c(1,3:6,8:10)]

#' Effect is always negative (that is expected), and of similar size. 
#' 
myTab[estimate == "sigma_lipLowMed" & comment == "OK",c(1,3:6,8:10)]

#' Effect is always positive (that is expected: introducing treatment does lower the values, and there is a high variability from before and during treatment). The effect size is not so similar, with
#'  
#' - men: between 0.07 and 0.15
#' - women: between 0.22 and 0.35
#' - post-menopausal women: between 0.32 and 0.39
#' - pre-menopausal women: 0.63
#' 
#' Not yet clear why though. I thought that in men lipid medication and age work in the same direction, hence I expected higher variability there. 
#' 

#' # Save data ####
#' ***
save(myTab, file = "../results/01_Prep_13_models.RData")

WriteXLS(x=myTab, 
         ExcelFileName="../results/01_Prep_13_models.xlsx", 
         SheetNames="NonGeneticModels", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
