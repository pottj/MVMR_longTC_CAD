#' ---
#' title: "Fix gamlss"
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

#' # Loop 1: get gene dosage matrix ####
#' ***
load("../results/UKB/05_potentialIVs_trajGWAS_Wald.RData")
load("../results/UKB/06_potentialIVs_GAMLSS.RData")

data_QC2 = gsub("~","../../",data_QC)
load(paste0(data_QC2,"UKB/temp_pgen/UKB_merged.RData"))

#' # Check previous results ####
#' ***
#' Because of some out of memory problem, the first script did not save all SNPs. Here, I just want to rerun those missing SNP - setting combinations. 
#' 
exposureData[,dumID := paste0(setting,rsID, sep="::")]
myAssocs[,dumID := paste0(setting,rsID, sep="::")]
exposureData2 = copy(exposureData)
exposureData2 = exposureData2[!is.element(dumID,myAssocs$dumID)]

#' # Loop 2: gamlss per SNP ####
#' ***
mySettings = unique(exposureData2$setting)
mySettings

registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab3 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  message("Working on setting ",mySettings[i]," (",i," of ",length(mySettings),")")
  
  # load phenotype data
  data1 = fread(paste0(data_QC,"/UKB/01_Prep_TrajGWAS_",mySettings[i],".csv"))
  
  # get genetic data for SNP subset
  exposureData3 = copy(exposureData2)
  exposureData3 = exposureData3[setting == mySettings[i],]
  
  dumTab2 = foreach(k = 1:dim(exposureData3)[1])%dorng%{
    #k=1
    mySNPID = exposureData3[k,rsID]
    myRow = pvar[ID == mySNPID,]
    filt = colnames(geno_mat) == mySNPID
    mySNP = geno_mat[,filt]
    
    # add SNP vector
    matched = match(data1$ID,psam$IID)
    data1[,myG := mySNP[matched]]
    
    # get some SNP info
    myRow[,EAF := sum(data1$myG)/(2*dim(data1)[1])]
    myRow[,MAF := EAF]
    myRow[EAF>0.5,MAF := 1-EAF]
    
    # run gamlss 
    mod2 = gamlss(TC ~ myG + age + age2 + statin +
                    PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10 + random(x = as.factor(ID)),   
                  sigma.formula = ~ myG + age + age2 + statin + 
                    PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10, 
                  data = na.omit(data1), family = "NO")
    
    dummy2 = summary(mod2)
    dummy2 = dummy2[grepl("myG",rownames(dummy2)),]
    
    n1= length(unique(data1[,ID]))
    myRow2 = copy(exposureData3)
    myRow2 = myRow2[setting == mySettings[i] & rsID ==myRow$ID]
    myRow2[,EAF := 1-myRow$EAF]
    myRow2[,mean_beta := dummy2[1,1] *(-1)]
    myRow2[,var_beta := dummy2[2,1] *(-1)]
    myRow2[,mean_SE := dummy2[1,2]]
    myRow2[,var_SE := dummy2[2,2]]
    myRow2[,mean_pval := dummy2[1,4]]
    myRow2[,var_pval := dummy2[2,4]]
    myRow2[,sampleSize := n1]
    myRow2[,method := "GAMLSS"]
    myRow2[,dumID := NULL]
    myRow2
  }
  myAssocs_3 = rbindlist(dumTab2)
  
  # Save temporarily
  message("Saving results for setting ",mySettings[i]," (",i," of ",length(mySettings),")")
  save(myAssocs_3, file = paste0("../temp/UKB_GAMLSS_",mySettings[i],"2.RData"))
  
  # return data table
  myAssocs_3
}
  
myAssocs_2 = rbindlist(dumTab3)
myAssocs[,dumID := NULL]
myAssocs = rbind(myAssocs,myAssocs_2)
setorder(myAssocs,chr,pos_b37)

save(myAssocs,file=paste0("../results/UKB/06_potentialIVs_GAMLSS_bugfixed.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
