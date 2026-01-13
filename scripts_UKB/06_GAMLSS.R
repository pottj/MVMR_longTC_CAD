#' ---
#' title: "Run gamlss"
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
load("../results/UKB/05_candidateSNPlist_pgen.RData")
mySNPs = copy(pvar)
myCHR = mySNPs[,unique(CHR)]

data_QC2 = gsub("~","../../",data_QC)
psam1 = fread(paste0(data_QC2,'UKB/temp_pgen/UKB_chr1.psam'))
myGenoMat = matrix(data = NA, ncol = 1,nrow = dim(psam1)[1])

dumTab1 = foreach(j = 1:length(myCHR))%do%{
  #j=1
  pvar1 = NewPvar(paste0(data_QC2, 'UKB/temp_pgen/UKB_chr',myCHR[j],'.pvar'))
  pgen = NewPgen(paste0(data_QC2,'UKB/temp_pgen/UKB_chr',myCHR[j],'.pgen'), pvar=pvar1)
  pvar = fread(paste0(data_QC2,'UKB/temp_pgen/UKB_chr',myCHR[j],'.pvar'))
  psam = fread(paste0(data_QC2,'UKB/temp_pgen/UKB_chr',myCHR[j],'.psam'))
  
  setnames(pvar,"#CHROM","CHR")
  pvar[,dumID := paste(CHR,POS,ALT,REF,sep=":")]
  
  myNRs = 1:dim(pvar)[1]
  filt = pvar$dumID %in% mySNPs$dumID
  pvar = pvar[filt,]
  
  geno_mat <- ReadList(pgen, variant_subset = myNRs[filt] , meanimpute=F)
  dim(geno_mat)
  
  # if(j>=19){
  #   filt3 = psam$IID %in% psam1$IID
  #   psam = psam[filt3,]
  #   geno_mat = geno_mat[filt3,]
  #   
  # }
  # checks
  stopifnot(pvar$ID == mySNPs[CHR==myCHR[j],ID])
  stopifnot(pvar$POS == mySNPs[CHR==myCHR[j],POS])
  stopifnot(pvar$ALT == mySNPs[CHR==myCHR[j],ALT])
  stopifnot(pvar$REF == mySNPs[CHR==myCHR[j],REF])
  
  colnames(geno_mat) = pvar[,ID]
  rownames(geno_mat) = psam$IID
  
  myGenoMat = cbind(myGenoMat,geno_mat)
  pvar
}
pvar = rbindlist(dumTab1)
head(pvar)

dim(myGenoMat)
myGenoMat[1:10,1:10]
geno_mat = myGenoMat[,-1]
dim(geno_mat)
geno_mat[1:10,1:10]

pvar[,EAF := colSums(geno_mat)/(2*dim(psam)[1])]
pvar[,MAF := EAF]
pvar[EAF>0.5,MAF := 1-EAF]

save(geno_mat,pvar,psam,file = paste0(data_QC2,"UKB/temp_pgen/UKB_merged.RData"))

#' # Loop 2: gamlss per SNP ####
#' ***
mySubpops = list.dirs(trajGWAS_UKB_results_WaldTest)
mySubpops = mySubpops[-1]
mySettings = gsub(".*//","",mySubpops)

exposureData2 = copy(exposureData)
registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

dumTab3 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  message("Working on setting ",mySettings[i]," (",i," of ",length(mySettings),")")
  
  # load phenotype data
  data1 = fread(paste0(data_QC,"/UKB/01_Prep_TrajGWAS_",mySettings[i],".csv"))
  
  dumTab2 = foreach(k = 1:dim(pvar)[1])%dorng%{
    #k=1
    mySNP = geno_mat[,k]
    myRow = pvar[k,]
  
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
    myRow2 = copy(exposureData2)
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
    myRow2
  }
  myAssocs_3 = rbindlist(dumTab2)
  
  # Save temporarily
  message("Saving results for setting ",mySettings[i]," (",i," of ",length(mySettings),")")
  save(myAssocs_3, file = paste0("../temp/UKB_GAMLSS_",mySettings[i]))
  
  # return data table
  myAssocs_3
}
  
myAssocs = rbindlist(dumTab3)
myAssocs

save(myAssocs,file=paste0("../results/UKB/06_potentialIVs_GAMLSS.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
