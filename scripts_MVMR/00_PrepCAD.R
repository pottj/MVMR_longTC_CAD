#' ---
#' title: "MVMR get CAD data"
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
#' UKB TrajGWAS results
load("../results/UKB/05_potentialIVs_trajGWAS_Wald.RData")
data_UKB = copy(exposureData)

#' INTERVAL TrajGWAS results
load("../results/INTERVAL/03_potentialIVs_trajGWAS_Wald.RData")
data_INT = copy(exposureData)

#' Restrict to SNPs available in both studies
SNPList_UKB = copy(data_UKB)
SNPList_UKB = SNPList_UKB[!duplicated(rsID)]
setorder(SNPList_UKB,chr,pos_b37)

SNPList_INT = copy(data_INT)
SNPList_INT = SNPList_INT[!duplicated(rsID)]
setorder(SNPList_INT,chr,pos_b37)

SNPList = SNPList_UKB[,c(1:8)]
SNPList = SNPList[rsID %in% SNPList_INT$rsID,]

#' Create dummy ID for matching to CAD data
SNPList[,dumID1 := paste0(chr,":",pos_b37,"_",OA,"_",EA)]
SNPList[,dumID2 := paste0(chr,":",pos_b37,"_",EA,"_",OA)]

#' Load CAD summary statistics
CAD_data = fread(SumStats_CAD_sex_stratified)
head(CAD_data)

table(SNPList$rsID %in% CAD_data$rsid_ukb)
table(SNPList$dumID1 %in% CAD_data$rs_number, SNPList$dumID2 %in% CAD_data$rs_number)

CAD = copy(CAD_data)
CAD = CAD[rs_number %in% SNPList$dumID1 | rs_number %in% SNPList$dumID2,]
SNPList = SNPList[rsID %in% CAD$rsid_ukb,]
CAD[duplicated(rsid_ukb)]
setorder(CAD,CHR,BP)

stopifnot(CAD$rsid_ukb == SNPList$rsID)
stopifnot(CAD$BP == SNPList$pos_b37)

filt = CAD$reference_allele != SNPList$EA
table(filt)
plot(SNPList$EAF, CAD$eaf)

save(CAD, file = "../results/MVMR/00_CAD_sexStrat.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
