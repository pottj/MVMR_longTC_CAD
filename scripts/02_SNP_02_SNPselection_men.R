#' ---
#' title: "SNP selection in men"
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
#' 1) Load data from some consortium (GLGC) --> n0
#' 2) Filter SNPs                           --> n1
#'    - p<5x10^-8
#'    - MAF>=1%
#'    - autosomal
#'    - rsID available
#'    - match to hg38 possible
#' 3) Extract genotypes from study (UKBB)  --> n2
#' 4) Extract summary statistics from other study (Aragam et al.) and restrict to overlap --> n3
#' 5) Priority pruning based on p-values from consortium (GLGC) and position (1 MB window) --> n4 = n_final
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP150.GRCh38))
source("../helperfunctions/getSmallestDist.R")

#' # Load GLGC data ####
#' ***
GLGC_data = fread(SumStats_TC_men, header=TRUE, sep="\t")
n0 = dim(GLGC_data)[1]
head(GLGC_data)
table(GLGC_data$CHROM)

#' # Filter GLGC data ####
#' ***
#' - pvalue < 5e-8
#' - EAF < 0.99 and EAF > 0.01
#' - autosomal SNPs (chr<23)
#' - rsID available
#' - match to hg38 possible
#' 
GLGC_data = GLGC_data[pvalue_neg_log10_GC >=-log10(5e-8),]
GLGC_data = GLGC_data[POOLED_ALT_AF>0.01 & POOLED_ALT_AF<0.99,]
GLGC_data = GLGC_data[CHROM<=22,]
GLGC_data = GLGC_data[grepl("rs",rsID)]
GLGC_data[,dumID := paste(CHROM,POS_b37,sep=":")]
dupPos = GLGC_data[duplicated(dumID)]
GLGC_data = GLGC_data[!is.element(dumID,dupPos$dumID)]
GLGC_data[,dumID := NULL]

save(GLGC_data, file = "../temp/02_SNP_02_GLGC_men.RData")
load("../temp/02_SNP_02_GLGC_men.RData")

mySNPs = GLGC_data$rsID
snps = SNPlocs.Hsapiens.dbSNP150.GRCh38
myLift = snpsById(snps, mySNPs, ifnotfound="drop")

myPos = pos(myLift)
myIDs = myLift$RefSNP_id

matched = match(mySNPs,myIDs)
GLGC_data[,pos_b38 := myPos[matched]]
GLGC_data = GLGC_data[!is.na(pos_b38)]

n1 = dim(GLGC_data)[1]
save(GLGC_data, file = "../temp/02_SNP_02_GLGC_men.RData")
load("../temp/02_SNP_02_GLGC_men.RData")

#' Add ID as in UKB data
GLGC_data[,dumID1 := paste0(CHROM,":",POS_b37,":",ALT,":",REF)]
GLGC_data[,dumID2 := paste0(CHROM,":",POS_b37,":",REF,":",ALT)]

#' Add ID as in Aragam data
GLGC_data[,dumID3 := paste0(CHROM,":",POS_b37,"_",ALT,"_",REF)]
GLGC_data[,dumID4 := paste0(CHROM,":",POS_b37,"_",REF,"_",ALT)]

#' Add ID as in FinnGen + UKB meta data
GLGC_data[,dumID5 := paste0(CHROM,":",pos_b38,":",ALT,":",REF)]
GLGC_data[,dumID6 := paste0(CHROM,":",pos_b38,":",REF,":",ALT)]

#' # Extract SNPs in exposure data ####
#' ***
pvar = fread(paste0(data_QC, '/02_SNP_01_UKB_TC_GLGC_merged.pvar'))
pvar_AF = fread(paste0(data_QC, '/02_SNP_01_UKB_TC_GLGC_merged_AF.afreq'))

stopifnot(pvar$ID == pvar_AF$ID)
stopifnot(pvar$REF == pvar_AF$REF)
stopifnot(pvar$ALT == pvar_AF$ALT)
pvar[,ALT_FREQ := pvar_AF$ALT_FREQS]

pvar[,ID2 := paste0(`#CHROM`,":",POS,":",REF,":",ALT)]
pvar = pvar[ID2 %in% c(GLGC_data$dumID1,GLGC_data$dumID2),]
GLGC_data = GLGC_data[dumID1 %in% pvar$ID2 | dumID2 %in% pvar$ID2,]
n2 = dim(GLGC_data)[1]

save(GLGC_data, file = "../temp/02_SNP_02_GLGC_men.RData")
load("../temp/02_SNP_02_GLGC_men.RData")

#' # Extract SNPs in outcome data ####
#' ***
#' ## 2-sample MR: Aragam et al (2022)
Aragam = fread(SumStats_CAD_sex_stratified)
setnames(Aragam,"rs_number","markername")
table(is.element(GLGC_data$dumID3,Aragam$markername))
table(is.element(GLGC_data$dumID4,Aragam$markername))

Aragam = Aragam[markername %in% GLGC_data$dumID3 | markername %in% GLGC_data$dumID4,]
GLGC_data = GLGC_data[dumID3 %in% Aragam$markername | dumID4 %in% Aragam$markername,]
n3 = dim(GLGC_data)[1]

#' ## 1-sample MR: UKB (as reported in the Neale lab, as I need the sex-stratified data)
UKB_CAD = fread(SumStats_CAD_UKB_men)
table(is.element(GLGC_data$dumID1,UKB_CAD$variant))
table(is.element(GLGC_data$dumID2,UKB_CAD$variant))

UKB_CAD = UKB_CAD[variant %in% GLGC_data$dumID1 | variant %in% GLGC_data$dumID2,]
GLGC_data = GLGC_data[dumID1 %in% UKB_CAD$variant | dumID2 %in% UKB_CAD$variant,]
n3b = dim(GLGC_data)[1]

#' # Perform priority pruning ####
#' ***
setorder(GLGC_data,CHROM,POS_b37)
myCHRs = unique(GLGC_data$CHROM)

result.22 = foreach(s2 = myCHRs) %do% {
  # s2 = myCHRs[1]
  subdata2 = copy(GLGC_data)
  subdata2 = subdata2[CHROM == s2, ]
  setkey(subdata2, POS_b37)
  
  if(dim(subdata2)[1]<=1){
    subdata2[, keep := T]
    subdata2[, NR_SNPs := 0]
  }else{
    subdata2[, keep := NA]
    subdata2[, NR_SNPs := as.numeric(NA)]
    
    smallestDist = getSmallestDist(subdata2[, POS_b37])
    while(smallestDist < 1000000) {
      #minP = min(subdata2[is.na(keep), p])
      maxLogP = max(subdata2[is.na(keep), pvalue_neg_log10_GC])
      myPOS = subdata2[maxLogP == pvalue_neg_log10_GC & is.na(keep), POS_b37]
      if(length(myPOS)>1){
        myPOS = myPOS[1]
      }
      subdata2[POS_b37 == myPOS, keep := T]
      
      #filter for SNPs that can stay within the set (outside the +- 1MB range or keep==T)
      myFilt = (subdata2[, POS_b37] < (myPOS - 1000000)) | 
        (subdata2[, POS_b37] > (myPOS + 1000000)) | 
        subdata2[, keep] 
      myFilt[is.na(myFilt)] = FALSE
      subdata2 = subdata2[myFilt == TRUE, ]
      
      subdata2[POS_b37 == myPOS, NR_SNPs := sum(myFilt==F)]
      smallestDist = getSmallestDist(subdata2[, POS_b37])
    }
    
    #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
    subdata2[is.na(keep), NR_SNPs := 0]
    subdata2[is.na(keep), keep := TRUE]
  }
  
  subdata2
}
GLGC_data2 = rbindlist(result.22)
n4 = dim(GLGC_data2)[1]

#' # Compare Allele frequencies ####
#' ***
#' I will use the exposure data as reference: everything should have the same effect allele as the pgen file (I dont want to turn the gamlss results later on).
#' 
#' ## Checking UKB data
pvar = pvar[ID2 %in% c(GLGC_data2$dumID1,GLGC_data2$dumID2),]
head(pvar)
setnames(pvar, "#CHROM","CHR")

table(GLGC_data2$dumID1==pvar$ID2,GLGC_data2$dumID2==pvar$ID2)
table(GLGC_data2$ALT==pvar$ALT, GLGC_data2$ALT==pvar$REF)
plot(GLGC_data2$POOLED_ALT_AF, pvar$ALT_FREQ)
abline(0,1)

filt = GLGC_data2$ALT==pvar$REF
table(filt)
GLGC_data2[, effect_allele := ALT]
GLGC_data2[filt, effect_allele := REF]
GLGC_data2[, other_allele := REF]
GLGC_data2[filt, other_allele := ALT]
GLGC_data2[, EAF := POOLED_ALT_AF]
GLGC_data2[filt, EAF := 1-POOLED_ALT_AF]
GLGC_data2[, beta := EFFECT_SIZE]
GLGC_data2[filt, beta := EFFECT_SIZE *(-1)]

table(GLGC_data2$effect_allele==pvar$ALT, GLGC_data2$other_allele==pvar$REF)
plot(GLGC_data2$EAF, pvar$ALT_FREQ)
abline(0,1)

#' ## Checking Aragam data
Aragam = Aragam[markername %in% GLGC_data2$dumID3 | markername %in% GLGC_data2$dumID4,]
setorder(Aragam,CHR,BP)
Aragam[,EA := reference_allele]
Aragam[,OA := other_allele]
Aragam[,EAF := male_eaf]

table(GLGC_data2$effect_allele==Aragam$EA,GLGC_data2$other_allele==Aragam$OA)
table(GLGC_data2$effect_allele==Aragam$OA,GLGC_data2$other_allele==Aragam$OA)
filt = GLGC_data2$effect_allele==Aragam$OA
table(filt)
Aragam[, effect_allele := EA]
Aragam[filt, effect_allele := OA]
Aragam[, other_allele := OA]
Aragam[filt, other_allele := EA]
Aragam[filt, EAF := 1-EAF]
Aragam[filt, beta := beta *(-1)]

table(GLGC_data2$effect_allele==Aragam$effect_allele, GLGC_data2$other_allele==Aragam$other_allele)
plot(GLGC_data2$EAF, Aragam$EAF)
abline(0,1)

#' ## Checking FinnGen + UKB data
UKB_CAD = UKB_CAD[variant %in% GLGC_data2$dumID1 | variant %in% GLGC_data2$dumID2,]
UKB_CAD[,ALT := gsub(".*:","",variant)]
UKB_CAD[,REF := gsub(".*[0123456789]:","",variant)]
UKB_CAD[,REF := gsub(":.*","",REF)]

table(GLGC_data2$effect_allele==UKB_CAD$ALT,GLGC_data2$other_allele==UKB_CAD$REF)
filt = GLGC_data2$effect_allele==UKB_CAD$REF
table(filt)
UKB_CAD[, effect_allele := ALT]
UKB_CAD[filt, effect_allele := REF]
UKB_CAD[, other_allele := REF]
UKB_CAD[filt, other_allele := ALT]
UKB_CAD[, EAF := minor_AF]
filt2 = UKB_CAD$minor_allele == UKB_CAD$other_allele
UKB_CAD[filt2, EAF := 1-EAF]
UKB_CAD[, beta := beta]
UKB_CAD[filt, beta := beta *(-1)]

table(GLGC_data2$effect_allele==UKB_CAD$effect_allele, GLGC_data2$other_allele==UKB_CAD$other_allele)
plot(GLGC_data2$EAF, UKB_CAD$EAF)
abline(0,1)

#' # Combine and save data 
#' ***
#' I want some finale file including 
#' 
#' - SNP information (rsID, chr, pos b37, pos b38, effect allele, other allele)
#' - matching ID for UKB and UKB eaf (as in pgen files)
#' - GLGC summary statistics (eaf, beta, se, t, p, n)
#' - Aragam summary statistics (eaf, beta, se, t, p, n)
#' - UKB summary statistics (eaf, beta, se, t, p, n)
#' 
names(GLGC_data2)
names(pvar)
names(Aragam)
names(UKB_CAD)

SNPList = data.table(rsID = GLGC_data2$rsID,
                     chr = GLGC_data2$CHROM, 
                     pos_b37 = GLGC_data2$POS_b37,
                     pos_b38 = GLGC_data2$pos_b38, 
                     effect_allele = GLGC_data2$effect_allele, 
                     other_allele = GLGC_data2$other_allele)

SNPList[, UKB_TC_ID := pvar$ID]
SNPList[, UKB_TC_EAF := pvar$ALT_FREQ]

SNPList[, GLGC_EAF := GLGC_data2$EAF]
SNPList[, GLGC_beta := GLGC_data2$beta]
SNPList[, GLGC_SE := GLGC_data2$SE]
SNPList[, GLGC_tval := GLGC_data2$beta / GLGC_data2$SE]
SNPList[, GLGC_logp := GLGC_data2$pvalue_neg_log10_GC]
SNPList[, GLGC_sampleSize := GLGC_data2$N]

SNPList[, Aragam_ID := Aragam$markername]
SNPList[, Aragam_EAF := Aragam$EAF]
SNPList[, Aragam_beta := Aragam$beta]
SNPList[, Aragam_SE := Aragam$male_se]
SNPList[, Aragam_tval := Aragam$male_z]
SNPList[, Aragam_pval := Aragam$male_p_value]
SNPList[, Aragam_sampleSize := Aragam$male_n_samples]

SNPList[, UKB_ID := UKB_CAD$variant]
SNPList[, UKB_EAF := UKB_CAD$EAF]
SNPList[, UKB_beta := UKB_CAD$beta]
SNPList[, UKB_SE := UKB_CAD$se]
SNPList[, UKB_tval := UKB_CAD$tstat]
SNPList[, UKB_pval := UKB_CAD$pval]
SNPList[, UKB_sampleSize := UKB_CAD$n_complete_samples]

save(SNPList,file = "../results/02_SNP_02_SNPListMen.RData")

n0;n1;n2;n3b;n4

#' # Create PLINK call ####
#' ***
#' I want to create a plink call to extract those 239 SNPs from the larger sample set
write.table(unique(SNPList$rsID),file = paste0(data_QC,"/02_SNP_02_SNPListMen.txt"), 
            col.names = F, row.names = F, quote = F)

call1 = paste0("plink2", 
               " --pfile ",data_QC, "/02_SNP_01_UKB_TC_GLGC_merged",
               " --extract ",data_QC,"/02_SNP_02_SNPListMen.txt", 
               " --keep-fam ",data_QC,"/01_Prep_02_SampleList_TC_men.txt",
               " --make-pgen --out ",data_QC,"/02_SNP_02_UKB_TC_GLGC_men")
print(call1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
