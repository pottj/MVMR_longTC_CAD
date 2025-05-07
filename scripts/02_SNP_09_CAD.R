#' ---
#' title: "Get association of SNPs on CAD"
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
#' Here, I want to estimate the SNP effects in a **glm** model for CAD. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()

#' # Load and prep UKB data ####
#' ***
#' Load data (genetic data & phenotype data)
load(paste0(data_QC,"/01_Prep_01_UKB_filtered_CAD.RData"))
load(paste0(data_QC,"/02_SNP_08_UKB_CAD.RData"))

stopifnot(pvar$ID == colnames(geno_mat))
stopifnot(psam$FID == rownames(geno_mat))
pvar[,rsID := ID]

#' # Get effects ####
#' ***
data0 = copy(myTab)
names(data0)[10:19] = gsub("_","",names(data0)[10:19])

# prepare loop over SNPs
mySNPs = pvar$ID
data0 = data0[ID %in% psam$FID,]
matched = match(data0$ID,psam$FID)
table(is.na(matched))

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
#counter = seq(1,450,25)

dumTab2 = foreach(i = 1:length(mySNPs))%dorng%{
#dumTab2 = foreach(i = 1:length(mySNPs))%do%{
  #i=1
  #if(i %in% counter) message("Working on SNP NR ",i)
  mySNP = geno_mat[,i]
  
  data2 = copy(data0)
  data2[,myG := mySNP[matched]]
  
  mod1 = glm(CAD ~ myG + age + lipidLow_init + Array + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             data = data2,sub = sex==1,family="binomial")
  
  mod2 = glm(CAD ~ myG + age + lipidLow_init + Array + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             data = data2,sub = sex==0,family="binomial")
  
  mod3 = glm(CAD ~ myG + age + lipidLow_init + Array + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             data = data2,sub = group=="postmenopausal",family="binomial")
  
  mod4 = glm(CAD ~ myG + age + lipidLow_init + Array + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             data = data2,sub = group=="premenopausal",family="binomial")
  
  dummy1 = summary(mod1)$coef
  dummy1 = dummy1[grepl("myG",rownames(dummy1)),]
  dummy2 = summary(mod2)$coef
  dummy2 = dummy2[grepl("myG",rownames(dummy2)),]
  dummy3 = summary(mod3)$coef
  dummy3 = dummy3[grepl("myG",rownames(dummy3)),]
  dummy4 = summary(mod4)$coef
  dummy4 = dummy4[grepl("myG",rownames(dummy4)),]
  dummy = rbind(dummy1,dummy2,dummy3,dummy4)
  
  n1 = length(mod1$fitted.values)
  n2 = length(mod2$fitted.values)
  n3 = length(mod3$fitted.values)
  n4 = length(mod4$fitted.values)
  
  res1 = data.table(SNP = rep(pvar[i,ID],4),
                    model = c("men","women","postmenopausal","premenopausal"),
                    sampleSize = c(n1,n2,n3,n4),
                    beta_mean = dummy[,1],
                    SE_mean =   dummy[,2],
                    tval_mean = dummy[,3],
                    pval_mean = dummy[,4])
  res1
}
myAssocs_Y_UKB = rbindlist(dumTab2)
myAssocs_Y_UKB

save(myAssocs_Y_UKB,file=paste0("../results/02_SNP_09_CAD_SumStats.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
