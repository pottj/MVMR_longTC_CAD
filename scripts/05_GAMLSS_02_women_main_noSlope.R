#' ---
#' title: "Get association of SNPs on TC (women - MAIN)"
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
#' Here, I want to estimate the SNP effects in a **gamlssIA** model for TC. 
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
load(paste0(data_QC,"/01_Prep_02_UKB_GP_allSamples.RData"))
load(paste0(data_QC,"/02_SNP_05_UKB_TC_GLGC_women_merged.RData"))

stopifnot(pvar$ID == colnames(geno_mat))
stopifnot(psam$FID == rownames(geno_mat))
pvar[,rsID := ID]

#' # Get effects ####
#' ***
matched = match(myTab6$BSU_ID,myTab7$ID)
myTab_long = cbind(myTab6,myTab7[matched,c(2,10:19,26)])
myTab_long[sex==0,sex:=2]
myTab_long[,exposure_type := NULL]

data1 = copy(myTab_long)
data1 = data1[sex==2,]
setnames(data1,"exposure_value","TC")
setnames(data1,"BSU_ID","ID")
names(data1)[13:22] = gsub("_","",names(data1)[13:22])

# prepare loop over SNPs
mySNPs = pvar$ID
data1 = data1[ID %in% psam$FID,]
matched = match(data1$ID,psam$FID)
table(is.na(matched))

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab2 = foreach(i = 1:length(mySNPs))%dorng%{
  #dumTab2 = foreach(i = 1:length(mySNPs))%do%{
  #i=1
  mySNP = geno_mat[,i]
  
  data2 = copy(data1)
  data2[,myG := mySNP[matched]]
  
  mod2 = gamlss(TC ~ myG + exposure_age + lipLowMed + Array +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + random(x = as.factor(ID)),   
                sigma.formula = ~ myG + exposure_age + lipLowMed + Array, 
                data = na.omit(data2), family = "NO")
  
  dummy2 = summary(mod2)
  dummy2 = dummy2[grepl("myG",rownames(dummy2)),]
  
  n1= length(unique(data2[,ID]))
  
  res1 = data.table(SNP = rep(pvar[i,ID],1),
                    model = c("women_main_noSlope"),
                    sampleSize = c(n1),
                    beta_mean = c(dummy2[1,1]),
                    SE_mean =   c(dummy2[1,2]),
                    tval_mean = c(dummy2[1,3]),
                    pval_mean = c(dummy2[1,4]),
                    beta_var = c(dummy2[2,1]),
                    SE_var =   c(dummy2[2,2]),
                    tval_var = c(dummy2[2,3]),
                    pval_var = c(dummy2[2,4]))
  res1
}
myAssocs_X = rbindlist(dumTab2)
myAssocs_X

save(myAssocs_X,file=paste0("../results/05_GAMLSS_02_women_main_noSlope.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
