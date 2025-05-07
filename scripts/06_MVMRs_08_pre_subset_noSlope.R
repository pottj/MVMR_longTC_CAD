#' ---
#' title: "MVMR - pre-menopausal women - subset - noSlope"
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
source("../helperfunctions/MVMR_jp_POPS.R")

#' # Get data ####
#' ***
#' ## SNP info
#' Get high quality SNPs only
load("../results/02_SNP_05_LD_women.RData")
LDTab[,SNP2 := as.character(SNP2)]
load("../results/02_SNP_03_SNPListWomen.RData")

#' ## Exposure
load("../results/05_GAMLSS_08_pre_subset_noSlope.RData")
length(unique(myAssocs_X$SNP))

#' Check the GX associations
myAssocs_X[,table(pval_mean==0)]
myAssocs_X[,table(pval_mean<5e-8,model)]
myAssocs_X[,table(pval_var<5e-8,model)]

myAssocs_X[,table(pval_mean<5e-8,pval_var<5e-8, model)]
myAssocs_X[,cor.test(beta_mean,beta_var),by=model]

#' Transform into long format
data_long1 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("beta_mean", "beta_var"),
                  variable.name="type",
                  value.name="beta")
data_long2 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("SE_mean", "SE_var"),
                  variable.name="type",
                  value.name="SE")
data_long3 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("tval_mean", "tval_var"),
                  variable.name="type",
                  value.name="tval")
data_long4 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:3],
                  measure.vars=c("pval_mean", "pval_var"),
                  variable.name="type",
                  value.name="pval")

myAssocs_X_long = cbind(data_long1,data_long2[,5],
                        data_long3[,5],data_long4[,5])
myAssocs_X_long[,type := gsub("beta_","",type)]
myAssocs_X_long = myAssocs_X_long[!is.na(beta),]

#' ## Outcome
data_long5 = melt(SNPList,
                  id.vars=names(SNPList)[1:3],
                  measure.vars=c("Aragam_beta", "UKB_beta"),
                  variable.name="type",
                  value.name="beta")
data_long6 = melt(SNPList,
                  id.vars=names(SNPList)[1:3],
                  measure.vars=c("Aragam_SE", "UKB_SE"),
                  variable.name="type",
                  value.name="SE")
data_long7 = melt(SNPList,
                  id.vars=names(SNPList)[1:3],
                  measure.vars=c("Aragam_tval", "UKB_tval"),
                  variable.name="type",
                  value.name="tval")
data_long8 = melt(SNPList,
                  id.vars=names(SNPList)[1:3],
                  measure.vars=c("Aragam_pval", "UKB_pval"),
                  variable.name="type",
                  value.name="pval")
data_long9 = melt(SNPList,
                  id.vars=names(SNPList)[1:3],
                  measure.vars=c("Aragam_sampleSize", "UKB_sampleSize"),
                  variable.name="type",
                  value.name="sampleSize")

myAssocs_Y_long = cbind(data_long9[,c(1,4,5)],data_long5[,5],data_long6[,5],
                        data_long7[,5],data_long8[,5])
myAssocs_Y_long[,type := gsub("_sampleSize","",type)]
myAssocs_Y_long = myAssocs_Y_long[!is.na(beta),]

load(paste0("../results/02_SNP_09_CAD_SumStats.RData"))
myAssocs_Y_UKB = myAssocs_Y_UKB[SNP %in% myAssocs_Y_long$rsID & model == "premenopausal"]
stopifnot(myAssocs_Y_long[type=="UKB",rsID] == myAssocs_Y_UKB$SNP)
myAssocs_Y_long[type=="UKB",beta := myAssocs_Y_UKB$beta_mean]
myAssocs_Y_long[type=="UKB",SE := myAssocs_Y_UKB$SE_mean]
myAssocs_Y_long[type=="UKB",tval := myAssocs_Y_UKB$tval_mean]
myAssocs_Y_long[type=="UKB",pval := myAssocs_Y_UKB$pval_mean]
myAssocs_Y_long[type=="UKB",sampleSize := myAssocs_Y_UKB$sampleSize]

#' ## save as temporary files
save(myAssocs_X_long,myAssocs_Y_long, file = paste0("../temp/06_MVMRInput_pre_subset_noSlope.RData"))

#' # Do MVMR ####
#' ***
#' I want 
#' 
#' - all SNPs (including all non-significant ones)
#' - genome-wide significant SNPs
#' 
myExposures = unique(myAssocs_X_long$model)
mySampleSize = unique(myAssocs_X$sampleSize)
myOutcomes = unique(myAssocs_Y_long$type)
myFlag = unique(myAssocs_X$model)

names(myAssocs_Y_long) = c("SNP", "phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean" )
myAssocs_X_long[,phenotype := paste0("TC_",myExposures)]
myAssocs_X_long[,dumID := myFlag]

dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
  #k=1
  myOutcome = myOutcomes[k]
  myAssocs_Y2 = copy(myAssocs_Y_long)
  myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
  
  message("Working on exposure ",myExposures," and outcome ",myOutcome," ...")
  
  # do MVMRs
  MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long,
                       data_outcome = myAssocs_Y2,
                       exposure_name = paste0("TC_",myExposures), 
                       outcome_name = myOutcome,
                       flag = myFlag,
                       GX_pval_treshold = 1,
                       getPlot = F,
                       corTab = LDTab,
                       corTab_threshold = 0.1,sampleSize_GX = mySampleSize,
                       random = F,getCondF = T,getUni = T)
  
  MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long,
                       data_outcome = myAssocs_Y2,
                       exposure_name = paste0("TC_",myExposures), 
                       outcome_name = myOutcome,
                       flag = myFlag,
                       GX_pval_treshold = 5e-8,
                       getPlot = F,
                       corTab = LDTab,
                       corTab_threshold = 0.1,sampleSize_GX = mySampleSize,
                       random = F,getCondF = T,getUni = T)
  
  MVMR0[,threshold := "all_SNPs"]
  MVMR2[,threshold := "gw_SNPs"]
  MVMR = rbind(MVMR0,MVMR2,fill=T)
  MVMR
  
}
MVMR_results = rbindlist(dumTab3,fill = T)

save(MVMR_results,file = paste0("../results/06_MVMR_08_pre_subset_noSlope.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
