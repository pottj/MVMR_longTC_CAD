#' ---
#' title: "Coloc between candidate SNP regions"
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
#' Step 1: get overlapping regions (pruning by position and GLGC p-value)
#' 
#' Step 2: load GLGC data
#' 
#' Step 3: run coloc per region
#' 
#' Step 4: save results
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()
source("../helperfunctions/getSmallestDist.R")

#' # Overlapping regions ####
#' ***
load("../results/02_SNP_02_SNPListMen.RData")
SNPList_men = copy(SNPList)

load("../results/02_SNP_03_SNPListWomen.RData")
SNPList_women = copy(SNPList)

SNPList_men[,sex := "male"]
SNPList_women[,sex := "female"]

SNPList = rbind(SNPList_men,SNPList_women)
setorder(SNPList,chr,pos_b37,-GLGC_logp)

# remove duplicate SNPs
SNPList = SNPList[!duplicated(rsID)]

myCHRs = unique(SNPList$chr)

result.22 = foreach(s2 = myCHRs) %do% {
  # s2 = myCHRs[1]
  subdata2 = copy(SNPList)
  subdata2 = subdata2[chr == s2, ]
  setkey(subdata2, pos_b37)
  
  if(dim(subdata2)[1]<=1){
    subdata2[, keep := T]
    subdata2[, NR_SNPs := 0]
  }else{
    subdata2[, keep := NA]
    subdata2[, NR_SNPs := as.numeric(NA)]
    
    smallestDist = getSmallestDist(subdata2[, pos_b37])
    while(smallestDist < 1000000) {
      #minP = min(subdata2[is.na(keep), p])
      maxLogP = max(subdata2[is.na(keep), GLGC_logp])
      myPOS = subdata2[maxLogP == GLGC_logp & is.na(keep), pos_b37]
      if(length(myPOS)>1){
        myPOS = myPOS[1]
      }
      subdata2[pos_b37 == myPOS, keep := T]
      
      #filter for SNPs that can stay within the set (outside the +- 1MB range or keep==T)
      myFilt = (subdata2[, pos_b37] < (myPOS - 1000000)) | 
        (subdata2[, pos_b37] > (myPOS + 1000000)) | 
        subdata2[, keep] 
      myFilt[is.na(myFilt)] = FALSE
      subdata2 = subdata2[myFilt == TRUE, ]
      
      subdata2[pos_b37 == myPOS, NR_SNPs := sum(myFilt==F)]
      smallestDist = getSmallestDist(subdata2[, pos_b37])
    }
    
    #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
    subdata2[is.na(keep), NR_SNPs := 0]
    subdata2[is.na(keep), keep := TRUE]
  }
  
  subdata2
}
SNPList2 = rbindlist(result.22)

#' Okay, there are 295 independent loci that I want to check for colocalization. My assumption: 
#' 
#' - if there was a SNP in the other sex within the region: potential coloc
#' - if there was no SNP in the other sex: potential sex-specific
#' 
#' But this could be wrong, because coloc does not need genome-wide significance (coloc even when the SNP is missing in the other sex).
#' 
#' # Load GLGC data ####
#' ***
GLGC_men = fread(SumStats_TC_men, header=TRUE, sep="\t")
GLGC_women = fread(SumStats_TC_women, header=TRUE, sep="\t")


#' # Colocalization ####
#' ***
dumTab1 = foreach(i = 1:dim(SNPList2)[1])%do%{
#dumTab1 = foreach(i = 1:10)%do%{
  #i=1 
  myRow = SNPList2[i,]
  message("Working on region around SNP ", myRow$rsID)
  
  # filter data
  data_1 = copy(GLGC_men)
  data_1 = data_1[CHROM == myRow$chr & POS_b37>= myRow$pos_b37-250000 & POS_b37<= myRow$pos_b37+250000]
  dim(data_1)
  
  data_2 = copy(GLGC_women)
  data_2 = data_2[CHROM == myRow$chr & POS_b37>= myRow$pos_b37-250000 & POS_b37<= myRow$pos_b37+250000]
  dim(data_2)
  
  # get overlap
  data_1[,dumID := paste(POS_b37,REF,ALT,sep="_")]
  data_2[,dumID := paste(POS_b37,REF,ALT,sep="_")]
  data_1 = data_1[dumID %in% data_2$dumID]
  data_2 = data_2[dumID %in% data_1$dumID]
  
  # exclude triallelic SNPs
  TriAllelicSNPs = data_1[duplicated(POS_b37),POS_b37]
  TriAllelicSNP_IDs = data_1[POS_b37 %in% TriAllelicSNPs,dumID]
  data_1 = data_1[!is.element(dumID,TriAllelicSNP_IDs)]
  data_2 = data_2[!is.element(dumID,TriAllelicSNP_IDs)]
  
  # order and check if same alleles were used
  setorder(data_1,POS_b37,REF)
  setorder(data_2,POS_b37,REF)
  stopifnot(data_1$dumID == data_2$dumID)
  stopifnot(data_1$REF == data_2$REF)
  stopifnot(data_1$ALT == data_2$ALT)
  #stopifnot(myRow$rsID %in% data_1$rsID)
  
  # check for NA in beta estimates (should have been filtered before!)
  filt<-is.na(data_1$EFFECT_SIZE) | is.na(data_2$EFFECT_SIZE)
  data_1<-data_1[!filt,]
  data_2<-data_2[!filt,]
  
  # add MAF
  data_1[,MAF := POOLED_ALT_AF]
  data_1[POOLED_ALT_AF>0.5,MAF := 1-POOLED_ALT_AF]
  data_2[,MAF := POOLED_ALT_AF]
  data_2[POOLED_ALT_AF>0.5,MAF := 1-POOLED_ALT_AF]
  
  # perform coloc
  my_res<- coloc::coloc.abf(dataset1=list(beta=data_1$EFFECT_SIZE,
                                          varbeta=(data_1$SE)^2,
                                          N=data_1$N,
                                          snp=data_1$dumID,
                                          MAF=data_1$MAF,
                                          type="quant"),
                            dataset2=list(beta=data_2$EFFECT_SIZE,
                                          varbeta=(data_2$SE)^2,
                                          N=data_2$N,
                                          snp=data_2$dumID,
                                          MAF=data_2$MAF,
                                          type="quant"))
  
  
  my_res2<-my_res[[1]]
  myRow[,nsnps := as.numeric(my_res2[1])]
  myRow[,PP.H0.abf := as.numeric(my_res2[2])]
  myRow[,PP.H1.abf := as.numeric(my_res2[3])]
  myRow[,PP.H2.abf := as.numeric(my_res2[4])]
  myRow[,PP.H3.abf := as.numeric(my_res2[5])]
  myRow[,PP.H4.abf := as.numeric(my_res2[6])]
  myRow
  }
SNPList3 = rbindlist(dumTab1)  

SNPList3[, table(PP.H4.abf>0.8)]
SNPList3[, table(PP.H3.abf>0.8)]
SNPList3[, table(PP.H2.abf>0.8)]
SNPList3[, table(PP.H1.abf>0.8)]
SNPList3[, table(PP.H0.abf>0.8)]

#' # Save ####
#' ***
save(SNPList3,file = "../results/02_SNP_06_SNPListColoc_complete.RData")

SNPList3 = SNPList3[PP.H4.abf>0.8]
save(SNPList3,file = "../results/02_SNP_06_SNPListColoc_filtered.RData")
write.table(unique(SNPList3$rsID),file = paste0(data_QC,"/02_SNP_06_SNPListColoc_filtered.txt"), 
            col.names = F, row.names = F, quote = F)

call1 = paste0("plink2", 
               " --pfile ",data_QC, "/02_SNP_01_UKB_TC_GLGC_merged",
               " --extract ",data_QC,"/02_SNP_06_SNPListColoc_filtered.txt", 
               " --keep-fam ",data_QC,"/01_Prep_01_SampleList_TC_GLGC.txt",
               " --make-pgen --out ",data_QC,"/02_SNP_06_UKB_TC_GLGC_coloc")
print(call1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
