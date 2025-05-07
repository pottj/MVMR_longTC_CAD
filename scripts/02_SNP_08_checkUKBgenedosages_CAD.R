#' ---
#' title: "Check UKB genetic data for CAD analyses"
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
#' Here, I want to load the genetic data. 
#' 
#' In my PLINK2 call, I already filter for imputation quality (mach r2 >0.8) and MAF (>0.01). In the SNP selection script, I already selected the variants I want to test and for the samples with longitudinal data.
#' 
#' Here, I just load the pgen data and store in an format for R for the later GLM runs. Then, I check per chromosome the pairwise LD.
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()

#' # Load genetic data ####
#' ***
data_QC2 = gsub("~","../../",data_QC)
pvar1 = NewPvar(paste0(data_QC2, '02_SNP_06_UKB_CAD.pvar'))
pgen = NewPgen(paste0(data_QC2,'02_SNP_06_UKB_CAD.pgen'), pvar=pvar1)
pvar = fread(paste0(data_QC2,'02_SNP_06_UKB_CAD.pvar'))
psam = fread(paste0(data_QC2,'02_SNP_06_UKB_CAD.psam'))

load("../results/02_SNP_06_SNPList_CAD.RData")
table(is.element(pvar$ID,SNPList$rsID))
table(is.element(SNPList$rsID,pvar$ID))
table(duplicated(pvar$ID))

#' There is one triallelic SNP: rs10127939. I will keep that variant with the matching alleles in the GLGC data set
#' 
pvar[ID == "rs10127939"]
SNPList[rsID == "rs10127939"]
filt = pvar$ID == "rs10127939" & pvar$REF == "C"
setnames(pvar,"#CHROM","CHR")

myNRs = 1:dim(pvar)[1]
geno_mat <- ReadList(pgen, variant_subset = myNRs , meanimpute=F)
dim(geno_mat)
colnames(geno_mat) = pvar[,ID]
rownames(geno_mat) = psam$IID
pvar = pvar[!filt,]
geno_mat = geno_mat[,!filt]

# filt samples
load(paste0(data_QC,"/01_Prep_01_UKB_filtered_CAD.RData"))
myTab[sex==0,sex:=2]
matched = match(psam$IID,myTab$ID)
myTab = myTab[matched,]
filt = !is.na(myTab$ID) & myTab$sex == psam$SEX
table(filt)
psam = psam[filt,]
geno_mat = geno_mat[filt,]
setnames(psam, "#FID","FID")

pvar[,EAF := colSums(geno_mat)/(2*dim(psam)[1])]
pvar[,MAF := EAF]
pvar[EAF>0.5,MAF := 1-EAF]
pvar[,table(MAF<0.01)]
pvar_AF = fread(paste0(data_QC, '/02_SNP_01_UKB_TC_GLGC_merged_AF.afreq'))
pvar_AF = pvar_AF[ID %in% pvar$ID]
filt = pvar_AF$ID == "rs10127939" & pvar_AF$REF == "C"
pvar_AF = pvar_AF[!filt,]
stopifnot(pvar$ID == pvar_AF$ID)
pvar[,EAF_fullUKB := pvar_AF$ALT_FREQS]
plot(pvar$EAF,pvar$EAF_fullUKB)
pvar[,EAF_fullUKB := NULL]

#' # LD checks ####
#' ***
pvar[,comment_LD := "LD OK"]
test = pvar[,.N,CHR]
test = test[N>1,]
myCHRs = test$CHR

LDTab = foreach(i = 1:length(myCHRs))%do%{
  #i=1
  myCHR = myCHRs[i]
  filt = pvar$CHR == myCHR
  dumMatrix = geno_mat[,filt]
  CorTab = cor(dumMatrix)^2
  heatmap(CorTab,Rowv = NA,Colv = NA, main =paste0("Chromosome ",myCHRs[i]))
  CorTab2 = as.data.table(CorTab)
  CorTab2[,SNP1 := rownames(CorTab)]
  CorTab_long = melt(CorTab2,id.vars="SNP1",measure.vars=rownames(CorTab),variable.name="SNP2")
  CorTab_long = CorTab_long[SNP1 != SNP2,]
  CorTab2_long = CorTab_long[value>0.05,]
  pvar[ID %in% CorTab2_long$SNP1 | ID %in% CorTab2_long$SNP2, comment_LD := "LD >0.05 with at least 1 SNP"]
  CorTab_long
}
LDTab = rbindlist(LDTab)
dev.off()
table(pvar$comment_LD)

#' Okay, there are many SNPs with no LD proxy, and some SNP which have at least one LD proxy with LD r2>0.05. 
#'
#' I will decide which SNP to use when I have the SNP associations. Then I can pick the best of these correlated variants. 
#' 
#' # Save data ####
#' ***
save(pvar, psam, geno_mat, file = paste0(data_QC,"02_SNP_08_UKB_CAD.RData"))
save(LDTab, file = paste0("../results/02_SNP_08_LD_CAD.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
