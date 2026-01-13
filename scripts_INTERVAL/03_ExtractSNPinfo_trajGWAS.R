#' ---
#' title: "Get associated SNPs"
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

#' # Load data ####
#' ***
load(file = "../results/UKB/03_GLGC_candidates.RData")
load(file = "../results/UKB/03_TrajGWAS_candidates.RData")

mySubpops = list.dirs(trajGWAS_INTERVAL_results_WaldTest)
mySubpops = mySubpops[-1]

dumTab1 = foreach(i = 1:length(mySubpops))%do%{
  #i=1
  myPath1 = mySubpops[i]
  myFiles = list.files(path = myPath1, pattern = "INTERVAL_chr")
  myCHR = gsub("INTERVAL_chr","",myFiles)
  myCHR = gsub(".pval.txt","",myCHR)
  myCHR = as.numeric(myCHR)
  myFiles = myFiles[order(myCHR)]
  myCHR = myCHR[order(myCHR)]
  
  dumTab2 = foreach(j = 1:length(myCHR))%do%{
    #j=1
    data2 = fread(paste0(myPath1,"/",myFiles[j]))
    data2
  }
  data1 = rbindlist(dumTab2)
  mySetting = gsub(".*//","",myPath1)
  data1[,setting := mySetting]
  
  data1[snpid %in% TrajGWAS_indep$snpid, TrajGWAS := T]
  data1[is.na(TrajGWAS), TrajGWAS := F]
  data1[snpid %in% TrajGWAS_indep[setting == mySetting,snpid], TrajGWAS_set := T]
  data1[is.na(TrajGWAS_set), TrajGWAS_set := F]
  data1[snpid %in% TC_data$rsID, GLGC := T]
  data1[is.na(GLGC), GLGC := F]
  if(mySetting %in% c("men","old","young"))  data1[snpid %in% TC_data[setting == "MALE", rsID], GLGC_set := T]
  if(mySetting %in% c("women","post","pre"))  data1[snpid %in% TC_data[setting == "FEMALE", rsID], GLGC_set := T]
  data1[is.na(GLGC_set), GLGC_set := F]
  
  save(data1,file = paste0(trajGWAS_INTERVAL_results_WaldTest, "/SumStats_",mySetting,".RData"))
  data1
}

exposureData = rbindlist(dumTab1)
exposureData[taupval<betapval,.N,by = setting]
exposureData[taupval>betapval,.N,by = setting]
exposureData[TrajGWAS_set==T,table(betapval<5e-8, taupval<1e-5,setting)]

setorder(exposureData,setting,chr,pos)

#' # Change column names
#' ***
#' According to github issue #38, the allele 2 is the effect allele. But this is not necessarly the minor allele. 
#' 
#' TrajGWAS allele 1 = OA
#' TrajGWAS allele 2 = EA
#' pvar REF = OA
#' pvar ALT = EA
#' 
data_QC2 = gsub("~","../../",data_QC)
test = exposureData[,.N,by = snpid]
test[,table(N)]

SNPList = copy(exposureData)
SNPList = SNPList[!duplicated(snpid)]
setorder(SNPList,chr,pos)
SNPList = SNPList[,c(1:6)]
SNPList[,dumID := paste(chr,pos,allele1,allele2,sep=":")]

dumTab6 = foreach(i = 1:22)%do%{
  #i=1
  pvar1 = NewPvar(paste0(data_QC2, 'INTERVAL/temp_pgen/INTERVAL_chr',i,'.pvar'))
  pgen = NewPgen(paste0(data_QC2,'INTERVAL/temp_pgen/INTERVAL_chr',i,'.pgen'), pvar=pvar1)
  pvar = fread(paste0(data_QC2,'INTERVAL/temp_pgen/INTERVAL_chr',i,'.pvar'))
  psam = fread(paste0(data_QC2,'INTERVAL/temp_pgen/INTERVAL_chr',i,'.psam'))
  
  setnames(pvar,"#CHROM","CHR")
  pvar[,dumID := paste(CHR,POS,ALT,REF,sep=":")]
  
  myNRs = 1:dim(pvar)[1]
  filt = pvar$dumID %in% SNPList$dumID
  pvar = pvar[filt,]
  
  geno_mat <- ReadList(pgen, variant_subset = myNRs[filt] , meanimpute=F)
  dim(geno_mat)
  
  # checks
  stopifnot(pvar$ID == SNPList[chr==i,snpid])
  stopifnot(pvar$POS == SNPList[chr==i,pos])
  stopifnot(pvar$ALT == SNPList[chr==i,allele1])
  stopifnot(pvar$REF == SNPList[chr==i,allele2])
  
  colnames(geno_mat) = pvar[,ID]
  rownames(geno_mat) = psam$IID
  
  pvar[,EAF := colSums(geno_mat)/(2*dim(psam)[1])]
  pvar[,MAF := EAF]
  pvar[EAF>0.5,MAF := 1-EAF]
  
  # REF = other allele, ALT = effect allele
  pvar[,EA := ALT]
  pvar[,OA := REF]
  pvar
  
}
pvar = rbindlist(dumTab6)
plot(pvar$MAF,SNPList$maf)
abline(0,1)
plot(pvar$EAF,SNPList$maf)
abline(0,1)

table(pvar$EA == SNPList$allele2)
pvar[,minorAllele := EA]
pvar[EAF>0.5,minorAllele := OA]
table(pvar$minorAllele == SNPList$allele2)
pvar[minorAllele != SNPList$allele2]
save(pvar,file = "../results/INTERVAL/03_candidateSNPlist_pgen.RData")

matched = match(exposureData$snpid,pvar$ID)
pvar = pvar[matched,]
stopifnot(pvar$ID == exposureData$snpid)
stopifnot(pvar$POS == exposureData[,pos])
stopifnot(pvar$ALT == exposureData[,allele1])
stopifnot(pvar$REF == exposureData[,allele2])
exposureData[,table(allele2 == pvar$minorAllele)]
exposureData[,EA := allele2]
exposureData[,OA := allele1]
exposureData[,EAF := maf]
exposureData[EA != pvar$minorAllele,EAF := 1-maf]
exposureData[,minorAllele := pvar$minorAllele]

names(exposureData)
exposureData = exposureData[,c(3,1,2,20,19,22, 6,21, 8:18)]
names(exposureData) = c("rsID","chr","pos_b37", "OA","EA","MA","MAF","EAF",
                        "mean_beta","mean_SE","mean_pval","var_beta","var_SE","var_pval",
                        "setting","TrajGWAS","TrajGWAS_set","GLGC","GLGC_set"    )

save(exposureData,file = "../results/INTERVAL/03_potentialIVs_trajGWAS_Wald.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
