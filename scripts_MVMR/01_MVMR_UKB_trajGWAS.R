#' ---
#' title: "MVMR UKB using trajGWAS estimates"
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
#' ## Men - trajGWAS SNPs ####
#' 
mySettings = c("men","old","young")
dumTab1 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  
  # get exposure data
  dataX = copy(exposureData)
  dataX = dataX[setting == mySettings[i]]
  dataX = dataX[TrajGWAS_set == T,]
  
  # get outcome data
  dataY = copy(CAD)
  dataY = dataY[rsid_ukb %in% dataX$rsID,]
  dataX = dataX[rsID %in% dataY$rsid_ukb,]
  
  # make sure that order is the same
  stopifnot(dataX$rsID == dataY$rsid_ukb)
  stopifnot(dataX$EA == dataY$reference_allele)
  
  # create MVMR object
  mvmr_obj = mr_mvinput(bx = as.matrix(dataX[,c(9,12)]),
                        bxse = as.matrix(dataX[,c(10,13)]),
                        by = dataY$male_beta, 
                        byse = dataY$male_se,
                        exposure = c("mean","variability"),
                        outcome = "CAD_males")
  
  res0 = mr_mvivw(mvmr_obj,nx = tab1[setting == mySettings[i],sampleSize]) 
  
  res1 = data.table(MR_model = rep("multivariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = rep(res0@Outcome,2),
                    NR_SNPs_total = rep(res0@SNPs,2),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-5])[1]),
                    beta_IVW = c(res0@Estimate),
                    SE_IVW = c(res0@StdError),
                    pval_IVW = c(res0@Pvalue),
                    condFstat = c(res0@CondFstat),
                    HeteroStat = rep(res0@Heter.Stat[1],2),
                    HeteroStat_pval = rep(res0@Heter.Stat[2],2))
  res1
 
  # create MR objects
  filt1 = dataX$mean_pval < 5e-8
  filt2 = dataX$var_pval <1e-5
  
  mr_obj1 = mr_input(bx = as.matrix(dataX[filt1,c(9,12)])[,1],
                     bxse = as.matrix(dataX[filt1,c(10,13)])[,1],
                     by = dataY$male_beta[filt1], 
                     byse = dataY$male_se[filt1],
                     exposure = "mean",
                     outcome = "CAD_males") 
  mr_obj2 = mr_input(bx = as.matrix(dataX[filt2,c(9,12)])[,2],
                     bxse = as.matrix(dataX[filt2,c(10,13)])[,2],
                     by = dataY$male_beta[filt2], 
                     byse = dataY$male_se[filt2],
                     exposure = "variability",
                     outcome = "CAD_males") 
  
  res2 = mr_ivw(mr_obj1)
  res3 = mr_ivw(mr_obj2)
  
  res4 = data.table(MR_model = rep("univariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = c(res2@Outcome,res3@Outcome),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-5])[1]),
                    beta_IVW = c(res2@Estimate,res3@Estimate),
                    SE_IVW = c(res2@StdError,res3@StdError),
                    pval_IVW = c(res2@Pvalue,res3@Pvalue),
                    condFstat = c(res2@Fstat,res3@Fstat),
                    HeteroStat = c(res2@Heter.Stat[1],res3@Heter.Stat[1]),
                    HeteroStat_pval = c(res2@Heter.Stat[2],res3@Heter.Stat[2]))
  res4
  
  # combine and return
  res = rbind(res1,res4,fill=T)
  res[,setting := mySettings[i]]
  res
}

MR_results_male_TrajGWAS = rbindlist(dumTab1)

#' ## Men - GLGC SNPs ####
#' 
dumTab1 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  
  # get exposure data
  dataX = copy(exposureData)
  dataX = dataX[setting == mySettings[i]]
  dataX = dataX[GLGC_set == T,]
  
  # get outcome data
  dataY = copy(CAD)
  dataY = dataY[rsid_ukb %in% dataX$rsID,]
  dataX = dataX[rsID %in% dataY$rsid_ukb,]
  
  # make sure that order is the same
  stopifnot(dataX$rsID == dataY$rsid_ukb)
  stopifnot(dataX$EA == dataY$reference_allele)
  
  # create MVMR object
  mvmr_obj = mr_mvinput(bx = as.matrix(dataX[,c(9,12)]),
                        bxse = as.matrix(dataX[,c(10,13)]),
                        by = dataY$male_beta, 
                        byse = dataY$male_se,
                        exposure = c("mean","variability"),
                        outcome = "CAD_males")
  
  res0 = mr_mvivw(mvmr_obj,nx = tab1[setting == mySettings[i],sampleSize]) 
  
  res1 = data.table(MR_model = rep("multivariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = rep(res0@Outcome,2),
                    NR_SNPs_total = rep(res0@SNPs,2),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-5])[1]),
                    beta_IVW = c(res0@Estimate),
                    SE_IVW = c(res0@StdError),
                    pval_IVW = c(res0@Pvalue),
                    condFstat = c(res0@CondFstat),
                    HeteroStat = rep(res0@Heter.Stat[1],2),
                    HeteroStat_pval = rep(res0@Heter.Stat[2],2))
  res1
  
  # create MR objects
  filt1 = dataX$mean_pval < 5e-8
  filt2 = dataX$var_pval <1e-5
  
  mr_obj1 = mr_input(bx = as.matrix(dataX[filt1,c(9,12)])[,1],
                     bxse = as.matrix(dataX[filt1,c(10,13)])[,1],
                     by = dataY$male_beta[filt1], 
                     byse = dataY$male_se[filt1],
                     exposure = "mean",
                     outcome = "CAD_males") 
  mr_obj2 = mr_input(bx = as.matrix(dataX[filt2,c(9,12)])[,2],
                     bxse = as.matrix(dataX[filt2,c(10,13)])[,2],
                     by = dataY$male_beta[filt2], 
                     byse = dataY$male_se[filt2],
                     exposure = "variability",
                     outcome = "CAD_males") 
  
  res2 = mr_ivw(mr_obj1)
  res3 = mr_ivw(mr_obj2)
  
  res4 = data.table(MR_model = rep("univariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = c(res2@Outcome,res3@Outcome),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-5])[1]),
                    beta_IVW = c(res2@Estimate,res3@Estimate),
                    SE_IVW = c(res2@StdError,res3@StdError),
                    pval_IVW = c(res2@Pvalue,res3@Pvalue),
                    condFstat = c(res2@Fstat,res3@Fstat),
                    HeteroStat = c(res2@Heter.Stat[1],res3@Heter.Stat[1]),
                    HeteroStat_pval = c(res2@Heter.Stat[2],res3@Heter.Stat[2]))
  res4
  
  # combine and return
  res = rbind(res1,res4,fill=T)
  res[,setting := mySettings[i]]
  res
}

MR_results_male_GLGC = rbindlist(dumTab1)

#' ## Women - trajGWAS SNPs ####
#' 
mySettings = c("women","post","pre")
dumTab1 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  
  # get exposure data
  dataX = copy(exposureData)
  dataX = dataX[setting == mySettings[i]]
  dataX = dataX[TrajGWAS_set == T,]
  
  # get outcome data
  dataY = copy(CAD)
  dataY = dataY[rsid_ukb %in% dataX$rsID,]
  dataX = dataX[rsID %in% dataY$rsid_ukb,]
  
  # make sure that order is the same
  stopifnot(dataX$rsID == dataY$rsid_ukb)
  stopifnot(dataX$EA == dataY$reference_allele)
  
  # create MVMR object
  mvmr_obj = mr_mvinput(bx = as.matrix(dataX[,c(9,12)]),
                        bxse = as.matrix(dataX[,c(10,13)]),
                        by = dataY$female_beta, 
                        byse = dataY$female_se,
                        exposure = c("mean","variability"),
                        outcome = "CAD_females")
  
  res0 = mr_mvivw(mvmr_obj,nx = tab1[setting == mySettings[i],sampleSize]) 
  
  res1 = data.table(MR_model = rep("multivariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = rep(res0@Outcome,2),
                    NR_SNPs_total = rep(res0@SNPs,2),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-5])[1]),
                    beta_IVW = c(res0@Estimate),
                    SE_IVW = c(res0@StdError),
                    pval_IVW = c(res0@Pvalue),
                    condFstat = c(res0@CondFstat),
                    HeteroStat = rep(res0@Heter.Stat[1],2),
                    HeteroStat_pval = rep(res0@Heter.Stat[2],2))
  res1
  
  # create MR objects
  filt1 = dataX$mean_pval < 5e-8
  filt2 = dataX$var_pval <1e-5
  
  mr_obj1 = mr_input(bx = as.matrix(dataX[filt1,c(9,12)])[,1],
                     bxse = as.matrix(dataX[filt1,c(10,13)])[,1],
                     by = dataY$female_beta[filt1], 
                     byse = dataY$female_se[filt1],
                     exposure = "mean",
                     outcome = "CAD_females") 
  mr_obj2 = mr_input(bx = as.matrix(dataX[filt2,c(9,12)])[,2],
                     bxse = as.matrix(dataX[filt2,c(10,13)])[,2],
                     by = dataY$female_beta[filt2], 
                     byse = dataY$female_se[filt2],
                     exposure = "variability",
                     outcome = "CAD_females") 
  
  res2 = mr_ivw(mr_obj1)
  res3 = mr_ivw(mr_obj2)
  
  res4 = data.table(MR_model = rep("univariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = c(res2@Outcome,res3@Outcome),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-5])[1]),
                    beta_IVW = c(res2@Estimate,res3@Estimate),
                    SE_IVW = c(res2@StdError,res3@StdError),
                    pval_IVW = c(res2@Pvalue,res3@Pvalue),
                    condFstat = c(res2@Fstat,res3@Fstat),
                    HeteroStat = c(res2@Heter.Stat[1],res3@Heter.Stat[1]),
                    HeteroStat_pval = c(res2@Heter.Stat[2],res3@Heter.Stat[2]))
  res4
  
  # combine and return
  res = rbind(res1,res4,fill=T)
  res[,setting := mySettings[i]]
  res
}

MR_results_female_TrajGWAS = rbindlist(dumTab1)

#' ## Women - GLGC SNPs ####
#' 
dumTab1 = foreach(i = 1:length(mySettings))%do%{
  #i=3
  
  # get exposure data
  dataX = copy(exposureData)
  dataX = dataX[setting == mySettings[i]]
  dataX = dataX[GLGC_set == T,]
  
  # get outcome data
  dataY = copy(CAD)
  dataY = dataY[rsid_ukb %in% dataX$rsID,]
  dataX = dataX[rsID %in% dataY$rsid_ukb,]
  
  # make sure that order is the same
  stopifnot(dataX$rsID == dataY$rsid_ukb)
  stopifnot(dataX$EA == dataY$reference_allele)
  
  # create MVMR object
  mvmr_obj = mr_mvinput(bx = as.matrix(dataX[,c(9,12)]),
                        bxse = as.matrix(dataX[,c(10,13)]),
                        by = dataY$female_beta, 
                        byse = dataY$female_se,
                        exposure = c("mean","variability"),
                        outcome = "CAD_females")
  
  res0 = mr_mvivw(mvmr_obj,nx = tab1[setting == mySettings[i],sampleSize]) 
  
  res1 = data.table(MR_model = rep("multivariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = rep(res0@Outcome,2),
                    NR_SNPs_total = rep(res0@SNPs,2),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-5])[1]),
                    beta_IVW = c(res0@Estimate),
                    SE_IVW = c(res0@StdError),
                    pval_IVW = c(res0@Pvalue),
                    condFstat = c(res0@CondFstat),
                    HeteroStat = rep(res0@Heter.Stat[1],2),
                    HeteroStat_pval = rep(res0@Heter.Stat[2],2))
  res1
  
  # create MR objects
  filt1 = dataX$mean_pval < 5e-8
  filt2 = dataX$var_pval <1e-5
  if(i==3) filt2 = dataX$var_pval <1e-2
  
  mr_obj1 = mr_input(bx = as.matrix(dataX[filt1,c(9,12)])[,1],
                     bxse = as.matrix(dataX[filt1,c(10,13)])[,1],
                     by = dataY$female_beta[filt1], 
                     byse = dataY$female_se[filt1],
                     exposure = "mean",
                     outcome = "CAD_females") 
  mr_obj2 = mr_input(bx = as.matrix(dataX[filt2,c(9,12)])[,2],
                     bxse = as.matrix(dataX[filt2,c(10,13)])[,2],
                     by = dataY$female_beta[filt2], 
                     byse = dataY$female_se[filt2],
                     exposure = "variability",
                     outcome = "CAD_females") 
  
  res2 = mr_ivw(mr_obj1)
  res3 = mr_ivw(mr_obj2)
  
  res4 = data.table(MR_model = rep("univariate",2),
                    exposure = c("TC mean","TC variability"),
                    outcome = c(res2@Outcome,res3@Outcome),
                    NR_SNPs_type = c(dim(dataX[mean_pval<5e-8])[1],dim(dataX[var_pval<1e-2])[1]),
                    beta_IVW = c(res2@Estimate,res3@Estimate),
                    SE_IVW = c(res2@StdError,res3@StdError),
                    pval_IVW = c(res2@Pvalue,res3@Pvalue),
                    condFstat = c(res2@Fstat,res3@Fstat),
                    HeteroStat = c(res2@Heter.Stat[1],res3@Heter.Stat[1]),
                    HeteroStat_pval = c(res2@Heter.Stat[2],res3@Heter.Stat[2]))
  res4
  
  # combine and return
  res = rbind(res1,res4,fill=T)
  res[,setting := mySettings[i]]
  res
}

MR_results_female_GLGC = rbindlist(dumTab1)

#' ## Combine the results ####
#' 
MR_results_female_TrajGWAS[,IV_selection := "trajGWAS"]
MR_results_male_TrajGWAS[,IV_selection := "trajGWAS"]
MR_results_female_GLGC[,IV_selection := "GLGC"]
MR_results_male_GLGC[,IV_selection := "GLGC"]

MR_results = rbind(MR_results_female_TrajGWAS,MR_results_male_TrajGWAS,MR_results_female_GLGC,MR_results_male_GLGC, fill = T)
save(MR_results, file = "../results/MVMR/01_MVMR_UKB_trajGWAS.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
