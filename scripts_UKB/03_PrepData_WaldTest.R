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
#' I want some template script to get the summary statistics of all sub-populations, and then filter for SNPs significant in at least one sub-population (either beta or tau). The resulting independent SNPs will then be used in the trajGWAS Wald test to estimate the SNP effect sizes. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
source("../helperfunctions/getSmallestDist.R")
source("../helperfunctions/positionBasedPruning.R")

#' # Load data ####
#' ***
mySubpops = list.dirs(trajGWAS_UKB_results)
mySubpops = mySubpops[-1]

dumTab1 = foreach(i = 1:length(mySubpops))%do%{
  #i=1
  myPath1 = mySubpops[i]
  myFiles = list.files(path = myPath1, pattern = "UKB_chr")
  myCHR = gsub("UKB_chr","",myFiles)
  myCHR = gsub(".pval.txt","",myCHR)
  myCHR = as.numeric(myCHR)
  myFiles = myFiles[order(myCHR)]
  myCHR = myCHR[order(myCHR)]

  dumTab2 = foreach(j = 1:length(myCHR))%do%{
    #j=1
    data2 = fread(paste0(myPath1,"/",myFiles[j]))
    data2[,jointpval2 := apply(.SD,1,harmonic.mean),.SDcols=c("betapval", "taupval")]
    data2
  }
  data1 = rbindlist(dumTab2)
  setting = gsub(".*//","",myPath1)
  data1[,setting := setting]

  # Save data for later use
  save(data1,file = paste0(trajGWAS_UKB_results, "/SumStats_",setting,".RData"))
  load(paste0(trajGWAS_UKB_results, "/SumStats_",setting,".RData"))
  
  # Create Manhattan plots
  data1[,pval := betapval]
  data1[betapval<1e-25,pval := 1e-25]
  filenm = paste0("../results/UKB/03_ManhattanPlots/",setting,"_beta.png")
  png(filename = filenm,width = 2250, height = 1125, res=200)
  manhattan(data1[-log10(pval)>2], chr="chr", bp="pos", snp="snpid", p="pval")
  dev.off()

  data1[,pval := taupval]
  data1[taupval<1e-25,pval := 1e-25]
  filenm = paste0("../results/UKB/03_ManhattanPlots/",setting,"_tau.png")
  png(filename = filenm,width = 2250, height = 1125, res=200)
  manhattan(data1[-log10(pval)>1], chr="chr", bp="pos", snp="snpid", p="pval" )
  dev.off()
  
  # calculate lambda 
  chisq_beta <- qchisq(1 - data1$betapval, 1)
  lambda_beta <- median(chisq_beta) / qchisq(0.5, 1)
  chisq_tau <- qchisq(1 - data1$taupval, 1)
  lambda_tau <- median(chisq_tau) / qchisq(0.5, 1)
  
  # QQ plots
  filenm = paste0("../results/UKB/03_QQPlots/",setting,"_beta.png")
  png(filename = filenm,width = 1125, height = 1125, res=200)
  qq(data1$betapval,main = paste0("UKB - ",setting," - beta, Inflation factor: ",round(lambda_beta,2)))
  dev.off()

  filenm = paste0("../results/UKB/03_QQPlots/",setting,"_tau.png")
  png(filename = filenm,width = 1125, height = 1125, res=200)
  qq(data1$taupval,main = paste0("UKB - ",setting," - tau, Inflation factor: ",round(lambda_tau,2)))
  dev.off()
  
  # return inflation factor and SNP size
  res = data.table(setting = setting,
                   lambda_beta = round(lambda_beta,5),
                   lambda_tau = round(lambda_tau,5),
                   nSNPs = dim(data1)[1],
                   nSNPs_MAF = dim(data1[maf>0.01,])[1])
  print(res)
  
  # Filter data for significant SNPs
  data1 = data1[betapval<5e-8 | taupval<1e-5,]
  data1[,pval := NULL]
  data1
  
}

TrajGWAS = rbindlist(dumTab1)
setorder(TrajGWAS,setting,chr,pos)
TrajGWAS = TrajGWAS[grepl("rs",snpid)]
TrajGWAS[,.N,by = c("setting")]
TrajGWAS[taupval<betapval,.N,by = setting]
TrajGWAS[taupval>betapval,.N,by = setting]

save(TrajGWAS,file = "../results/UKB/03_SumStatsFiltered.RData")
load("../results/UKB/03_SumStatsFiltered.RData")

#' # Position based pruning ####
#' ***
#' Now I perform position based pruning per setting and trait (beta/mean and tau/variability). I will do this per chromosome, e.g. select the best associated SNP of the chromosome (lowest p-value), and then remove all SNPs within +/- 500,000 kb. I will keep track of the number of SNPs within the window as information about the hits stability (more tagged SNPs higher stability/support). However, I will not filter for this information. 
#'  

mySubpops = unique(TrajGWAS$setting)

dumTab3 = foreach(i = 1:length(mySubpops))%do%{
  #i=1
  data3 = copy(TrajGWAS)
  data3 = data3[setting == mySubpops[i],]
  
  data4 = positionBasedPruning(data = data3, distance = 500000, 
                               col_chr = "chr",col_pos = "pos",col_prio = "jointpval2",
                               smallerIsBetter = T)
  data4
}

TrajGWAS_indep = rbindlist(dumTab3)
hist(TrajGWAS_indep$NR_SNPs)
TrajGWAS_indep[,.N,by = c("setting")]
TrajGWAS_indep[NR_SNPs>0,.N,by = c("setting")]
TrajGWAS_indep[NR_SNPs>4,.N,by = c("setting")]
TrajGWAS_indep[NR_SNPs>4 & jointpval2<5e-8,.N,by = c("setting")]
TrajGWAS_indep[NR_SNPs>4 & jointpval2<1e-6 & taupval<betapval,.N,by = c("setting")]
save(TrajGWAS_indep,file = "../results/UKB/03_independentIVs.RData")

#' # Colocalization within sub-population
#' ***
#' For 16 loci, I want to check if the beta and tau signals are one shared signal (colocalizing) or two independent signals. As I do not have the beta estimates (and their standard errors), I will use the p-values as input. Hence, I need the MAF and sample size  
#' 
coloc_candidates = copy(TrajGWAS_indep)
coloc_candidates = coloc_candidates[betapval<5e-8 & taupval<1e-5,]

mySubpops = unique(coloc_candidates$setting)
load("../results/UKB/01_descriptiveTable.RData")

ToDoList = res[setting %in% mySubpops,]
setorder(ToDoList,setting)
stopifnot(ToDoList$setting == mySubpops)

dumTab1 = foreach(i = 1:length(mySubpops))%do%{
  #i=2
  coloc_candidates2 = copy(coloc_candidates)
  coloc_candidates2 = coloc_candidates2[setting == mySubpops[i],]
  
  load(paste0(trajGWAS_UKB_results, "/SumStats_",mySubpops[i],".RData"))
  
  dumTab2 = foreach(j = 1:dim(coloc_candidates2)[1])%do%{
    #j=2
    myHit = copy(coloc_candidates2)
    myHit = myHit[j,]
    
    data2 = copy(data1)
    data2 = data2[chr == myHit$chr & pos > myHit$pos - 500000 & pos < myHit$pos + 500000,]
    data2 = data2[grepl("rs",snpid)]
    dupIDs = data2[duplicated(snpid),snpid]
    #data2[snpid %in% dupIDs,]
    data2 = data2[!is.element(snpid,dupIDs),]
    
    #data2[,meanlogP := log10(betapval) * betadir]
    #data2[,varlogP := log10(taupval) * taudir]
    #plot(data2$meanlogP,data2$varlogP)
    
    list1 = list(snp = data2$snpid,
                 position = data2$pos,
                 MAF= data2$maf, 
                 type = "quant",
                 #sdY = ToDoList[i,sd],
                 N = ToDoList[i,sampleSize], 
                 pvalues = data2$betapval)
    
    list2 = list(snp = data2$snpid,
                 position = data2$pos,
                 MAF  = data2$maf, 
                 type = "quant",
                 #sdY = ToDoList[i,sd],
                 N = ToDoList[i,sampleSize], 
                 pvalues = data2$taupval)
    
    #check_dataset(list1)
    #check_dataset(list2)
    
    my.res <- coloc.abf(dataset1=list1,
                        dataset2=list2)
    res = my.res$summary
    myHit[,nsnps := res[1]]
    myHit[,PP.H0.abf := res[2]]
    myHit[,PP.H1.abf := res[3]]
    myHit[,PP.H2.abf := res[4]]
    myHit[,PP.H3.abf := res[5]]
    myHit[,PP.H4.abf := res[6]]
    myHit
    
  }
  coloc_subPop = rbindlist(dumTab2)
  coloc_subPop
}
coloc_results = rbindlist(dumTab1)
coloc_results[PP.H4.abf<0.8]
save(coloc_results,file = "../results/UKB/03_coloc_within.RData")

#' Okay, there is one hit on chromosome 8 (near TRIB1AL), that does not perfectly co-localize, but it has still a PP4>50%. I will just treat them all as shared hits. 
#' 
TrajGWAS_indep[betapval<5e-8 & taupval<1e-5, comment:= "shared signal"]
TrajGWAS_indep[betapval<5e-8 & taupval>1e-5, comment:= "beta only"]
TrajGWAS_indep[betapval>5e-8 & taupval<1e-5, comment:= "tau only"]
TrajGWAS_indep[,table(comment)]
TrajGWAS_indep[NR_SNPs>0,table(comment)]
TrajGWAS_indep[NR_SNPs>0,table(comment,setting)]
TrajGWAS_indep = TrajGWAS_indep[NR_SNPs>0,]

length(unique(TrajGWAS_indep$snpid))
save(TrajGWAS_indep,file = "../results/UKB/03_TrajGWAS_candidates.RData")

#' # Compare with GLGC ####
#' ***
#' I load the GLGC sex-stratified data (Kanoni et al. 2022), and filter for
#' 
#' - genome-wide significance (p<5e-8)
#' - rsID available
#' - MAF>1%
#' - no duplicated chr:pos (in filtered data set)
#' - position based pruning as before
#' 
myFiles = c(SumStats_TC_men,SumStats_TC_women)
myFiles

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  mySetting = gsub(".meta.singlevar.results.gz","",myFiles[i])
  mySetting = gsub(".*_","",mySetting)
  message("Working on ",mySetting)
  
  data = fread(myFiles[i])
  message("     Initial SNP number: ",dim(data)[1])
  
  # some filtering
  data = data[pvalue_neg_log10_GC >=-log10(5e-8),]
  data = data[!is.na(rsID)]
  data = data[POOLED_ALT_AF >= 0.01]
  data = data[POOLED_ALT_AF <= 0.99]
  data[,dumID := paste(CHROM,POS_b37,sep=":")]
  dupPos = data[duplicated(dumID)]
  data = data[!is.element(dumID,dupPos$dumID)]
  data[,dumID := NULL]
  message("     SNP number after filtering: ",dim(data)[1])
  
  # add column for sex-setting
  data[,setting := mySetting]
  
  # pruning 
  data5 = positionBasedPruning(data = data, distance = 500000, 
                               col_chr = "CHROM",col_pos = "POS_b37",col_prio = "pvalue_neg_log10_GC",
                               smallerIsBetter = F)
  message("     SNP number after pruning: ",dim(data5)[1])
  
  # return data table
  data5
}
TC_data = rbindlist(dumTab1)
dim(TC_data)
setorder(TC_data,CHROM,POS_b37,setting)
TC_data[,.N,by = c("setting")]
TC_data[,length(unique(rsID))]

table(TC_data$rsID %in% TrajGWAS_indep$snpid)
table(TrajGWAS_indep$snpid %in% TC_data$rsID)
save(TC_data,file = "../results/UKB/03_GLGC_candidates.RData")

#' Okay, now I have a list of 246 unique SNPs from the TrajGWAS, and 554 SNPs from GLGC, which I will both test using trajGWAS and gamlss. 
#' 
candidateSNPs = c(TrajGWAS_indep[,snpid],TC_data$rsID)
candidateSNPs = candidateSNPs[!duplicated(candidateSNPs)]

write.table(candidateSNPs,
            file = paste0(data_QC,"/UKB/03_candidateSNPs.txt"), 
            col.names = F, row.names = F, quote = F)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
