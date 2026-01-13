#' ---
#' title: "MVMR UKB using trajGWAS estimates per sex"
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
#' UKB sample sizes per setting
load("../results/UKB/01_descriptiveTable.RData")
tab1 = copy(res)

#' TrajGWAS results
load("../results/UKB/05_potentialIVs_trajGWAS_Wald.RData")

#' CAD data 
load("../results/MVMR/00_CAD_sexStrat.RData")

#' # MVMR ####
#' ***
#' Here, I want to test the age-groups against each other: young vs old, pre vs post menopausal
#' 
ToDoList = data.table(setting1 = rep(c("young","pre"),2),
                      setting2 = rep(c("old","post"),2),
                      IV_selection = rep(c("TrajGWAS_set","GLGC_set"),each=2),
                      sex = c("men","women","men","women"))
ToDoList

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  
  # get exposure data
  dataX = copy(exposureData)
  dataX = dataX[setting %in% c(myRow$setting1,myRow$setting2),]
  dataX[,IVs := get(myRow$IV_selection)]
  mySNPList = dataX[IVs ==T,unique(rsID)]
  dataX = dataX[rsID %in% mySNPList,]
  
  # get outcome data
  dataY = copy(CAD)
  dataY = dataY[rsid_ukb %in% dataX$rsID,]
  dataX = dataX[rsID %in% dataY$rsid_ukb,]
  
  # make sure that order is the same
  stopifnot(dataX$rsID == dataY$rsid_ukb)
  stopifnot(dataX$EA == dataY$reference_allele)
  
  # create MVMR object
  beta_matrix1 = as.matrix(dataX[setting==myRow$setting1,c(9,12)])
  beta_matrix2 = as.matrix(dataX[setting==myRow$setting2,c(9,12)])
  beta_matrix = cbind(beta_matrix1,beta_matrix2)

  se_matrix1 = as.matrix(dataX[setting==myRow$setting1,c(10,13)])
  se_matrix2 = as.matrix(dataX[setting==myRow$setting2,c(10,13)])
  se_matrix = cbind(se_matrix1,se_matrix2)
  
  myExposures1 = paste(myRow$setting1,c("mean","variability"),sep = " - ")
  myExposures2 = paste(myRow$setting2,c("mean","variability"),sep = " - ")
  
  if(myRow$sex == "men"){
    myBY = dataY$male_beta
    myBYse = dataY$male_se
    myOutcome = "CAD_male"
  }else{
    myBY = dataY$female_beta
    myBYse = dataY$female_se
    myOutcome = "CAD_female"
  }
  mvmr_obj = mr_mvinput(bx = beta_matrix,
                        bxse = se_matrix,
                        by = myBY, 
                        byse = myBYse,
                        exposure = c(myExposures1,myExposures2),
                        outcome = myOutcome)
  
  n1 = tab1[setting == myRow$setting1,sampleSize]
  n2 = tab1[setting == myRow$setting2,sampleSize]
  
  res0 = mr_mvivw(mvmr_obj,nx = c(n1,n1,n2,n2)) 
  
  res1 = data.table(MR_model = rep("multivariate",4),
                    exposure = res0@Exposure,
                    outcome = rep(res0@Outcome,4),
                    NR_SNPs_total = rep(res0@SNPs,4),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8 & setting == myRow$setting1])[1],
                                     dim(dataX[var_pval<1e-5 & setting == myRow$setting1])[1],
                                     dim(dataX[mean_pval<5e-8 & setting == myRow$setting2])[1],
                                     dim(dataX[var_pval<1e-5 & setting == myRow$setting2])[1]),
                    beta_IVW = c(res0@Estimate),
                    SE_IVW = c(res0@StdError),
                    pval_IVW = c(res0@Pvalue),
                    condFstat = c(res0@CondFstat),
                    HeteroStat = rep(res0@Heter.Stat[1],4),
                    HeteroStat_pval = rep(res0@Heter.Stat[2],4),
                    IV_selection = myRow$IV_selection,
                    setting = myRow$sex)
  res1
  
}

MR_results = rbindlist(dumTab)
save(MR_results, file = "../results/MVMR/05_MVMR_UKB_trajGWAS_ageGroups.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
