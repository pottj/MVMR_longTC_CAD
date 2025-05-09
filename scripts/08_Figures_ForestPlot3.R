#' ---
#' title: "Figures - forest plots 3"
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
#' I want two Forest Plots here, comparing men and women in the no Slope model. I want all exposure types in one plot. Restrict to 
#' 
#' 1) Men and women (age groups)
#' 2) Co-localizing loci 
#' 3) 2-sample MVMR (Aragam outcome data)
#' 
#' Plot 1: all SNPs 
#' 
#' Plot 2: use genome-wide significant SNPs only
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # Load data ####
#' ***
myMVMR6_results = list.files(path = "../results/",pattern = "06_MVMR")
myMVMR6_results = myMVMR6_results[c(11:14)]
myMVMR_results = c(myMVMR6_results)

dumTab4 = foreach(i = 1:length(myMVMR_results))%do%{
  #i=1
  load(paste0("../results/",myMVMR_results[i]))
  names(MVMR_results)
  MVMR_results = MVMR_results[setting == "multivariate",]
  MVMR_results[outcome == "UKB", setting := "1-sample"]
  MVMR_results[outcome != "UKB", setting := "2-sample"]
  MVMR_results[outcome != "UKB", outcome := "Aragam_CAD"]
  MVMR_results[outcome == "UKB", outcome := "UKB_CAD"]
  
  MVMR_tab_wide = dcast(MVMR_results, ID + setting + exposure + outcome + threshold + NR_SNPs_total + HeteroStat + HeteroStat_pval ~ exposure_type, 
                        value.var=c("NR_SNPs_type","beta_IVW","SE_IVW","pval_IVW","condFstat"))
  names(MVMR_tab_wide) = gsub("IVW_","",names(MVMR_tab_wide))
  names(MVMR_tab_wide) = gsub("NR_SNPs_type","SNPs",names(MVMR_tab_wide))
  setnames(MVMR_tab_wide,"ID","flag")
  MVMR_tab_wide[,exposure := "TC"]
  
  MVMR_tab_wide
  
}
tab3 = rbindlist(dumTab4,fill=T)
names(tab3)
x1 = grep("mean",names(tab3))
x2 = grep("slope",names(tab3))
x3 = grep("var",names(tab3))
x = c(1:8,x1,x2,x3)
tab3 = tab3[,x,with=F]
tab3[,sex := gsub("_.*","",flag)]
tab3[,flag2 := gsub(".*_","",flag)]
tab3[flag2 == "coloc",flag2 := "main"]
tab3 = tab3[outcome == "Aragam_CAD"]
tab3[1:4,sex := "women"]
tab3[c(1,2,7,8),flag2 := "old"]
tab3[c(3:6),flag2 := "young"]

#' # Plot 1: all SNPs ####
#' ***
#' All SNP approach, 2-sample MVMR, all 12 settings
#' 
myTab = copy(tab3)
myTab = myTab[threshold == "all_SNPs",]

myTab2_SNPs = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
             measure.vars = c("SNPs_mean", "SNPs_var"))
myTab2_SNPs[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_SNPs,"value","SNPs")
myTab2_SNPs[,variable := gsub("SNPs_","",variable)]

myTab2_beta = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                   measure.vars = c("beta_mean", "beta_var"))
myTab2_beta[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_beta,"value","beta")
myTab2_beta[,variable := gsub("beta_","",variable)]

myTab2_SE = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                   measure.vars = c("SE_mean", "SE_var"))
myTab2_SE[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_SE,"value","SE")
myTab2_SE[,variable := gsub("SE_","",variable)]

myTab2_pval = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                 measure.vars = c("pval_mean", "pval_var"))
myTab2_pval[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_pval,"value","pval")
myTab2_pval[,variable := gsub("pval_","",variable)]

myTab2_condFstat = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                   measure.vars = c("condFstat_mean", "condFstat_var"))
myTab2_condFstat[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_condFstat,"value","condF")
myTab2_condFstat[,variable := gsub("condFstat_","",variable)]

myTab2 = cbind(myTab2_SNPs,myTab2_beta$beta,myTab2_SE$SE,myTab2_pval$pval,myTab2_condFstat$condF)

names(myTab2)[9:12] = c("beta","SE","pval","condF")
myTab2 = myTab2[!is.na(beta),]

myTab2[,lowerCI95 := beta-1.96*SE]
myTab2[,upperCI95 := beta+1.96*SE]

myTab2[,condF := round(condF,2)]
myTab2[,condF := as.character(condF)]

myTab2$` ` <- paste(rep(" ", 50), collapse = " ")
myTab2$`Estimate \n[95% CI]` <- ifelse(is.na(myTab2$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 myTab2$beta, myTab2$lowerCI95,myTab2$upperCI95))
setorder(myTab2, variable, sex, flag2)

#' Plot per exposure type
myTab2[,Setting := paste0("    ",sex, " - ",flag2)]
dummy = data.table(Setting = c("mean","variability"))
plotData = rbind(myTab2,dummy, fill=T)
plotData = plotData[c(9,1:4,10,5:8)]

plotData[is.na(flag),` ` := ""]
plotData[is.na(flag),`Estimate \n[95% CI]`:= ""]
plotData[is.na(flag),condF := ""]

dummy = plotData$Setting
dummy[c(1,6)] = "white"
dummy[c(2,7)] = "#BDD3E5"
dummy[c(3,8)] = "#BCF1FF"
dummy[c(4,9)] = "#B7FFB7"
dummy[c(5,10)] = "#FFECC9"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Causal effect of the of TC on CAD per exposure type (all candidate SNPs)")
setnames(plotData,"condF","cond \nF-stat")
setnames(plotData,"Setting", "Exposure type \n sex - age group")

p2<- forest(plotData[,c(17,15,16,12)],
            est = plotData$beta,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/ForestPlots/3_AllSNPs_2sample_MenWomenAgeGroups_noSlope_perXType.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' Plot per sample group
myTab2[,Setting := variable]
myTab2[Setting=="var",Setting := "variability"]
dummy = data.table(Setting = c("older men","younger men","post-menopausal women","pre-menopausal women"))
plotData = rbind(myTab2,dummy, fill=T)
plotData = plotData[c(9,1,5,10,2,6,11,3,7,12,4,8)]

plotData[is.na(flag),` ` := ""]
plotData[is.na(flag),`Estimate \n[95% CI]`:= ""]
plotData[is.na(flag),condF := ""]

dummy = plotData$Setting
dummy[c(1,4,7,10)] = "white"
dummy[c(2,3)] = "#BDD3E5"
dummy[c(5,6)] = "#BCF1FF"
dummy[c(8,9)] = "#B7FFB7"
dummy[c(11,12)] = "#FFECC9"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Causal effect of the of TC on CAD per sex and age group (all candidate SNPs)")
setnames(plotData,"condF","cond \nF-stat")
setnames(plotData,"Setting", "Sex and age group \n exposure type")

p2<- forest(plotData[,c(17,15,16,12)],
            est = plotData$beta,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/ForestPlots/3_AllSNPs_2sample_MenWomenAgeGroups_noSlope_perGroup.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # Plot 2: genome-wide significant SNPs ####
#' ***
#' Genome-wide significant SNPs approach, 2-sample MVMR, all 12 settings
#' 
myTab = copy(tab3)
myTab = myTab[threshold != "all_SNPs",]

myTab2_SNPs = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                   measure.vars = c("SNPs_mean", "SNPs_var"))
myTab2_SNPs[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_SNPs,"value","SNPs")
myTab2_SNPs[,variable := gsub("SNPs_","",variable)]

myTab2_beta = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                   measure.vars = c("beta_mean", "beta_var"))
myTab2_beta[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_beta,"value","beta")
myTab2_beta[,variable := gsub("beta_","",variable)]

myTab2_SE = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                 measure.vars = c("SE_mean", "SE_var"))
myTab2_SE[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_SE,"value","SE")
myTab2_SE[,variable := gsub("SE_","",variable)]

myTab2_pval = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                   measure.vars = c("pval_mean", "pval_var"))
myTab2_pval[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_pval,"value","pval")
myTab2_pval[,variable := gsub("pval_","",variable)]

myTab2_condFstat = melt(myTab, id.vars = names(myTab)[c(1,6:8,19,20)],
                        measure.vars = c("condFstat_mean", "condFstat_var"))
myTab2_condFstat[,flag := paste(sex,flag2,sep="_")]
setnames(myTab2_condFstat,"value","condF")
myTab2_condFstat[,variable := gsub("condFstat_","",variable)]

myTab2 = cbind(myTab2_SNPs,myTab2_beta$beta,myTab2_SE$SE,myTab2_pval$pval,myTab2_condFstat$condF)

names(myTab2)[9:12] = c("beta","SE","pval","condF")
myTab2 = myTab2[!is.na(beta),]

myTab2[,lowerCI95 := beta-1.96*SE]
myTab2[,upperCI95 := beta+1.96*SE]

myTab2[,condF := round(condF,2)]
myTab2[,condF := as.character(condF)]

myTab2$` ` <- paste(rep(" ", 50), collapse = " ")
myTab2$`Estimate \n[95% CI]` <- ifelse(is.na(myTab2$SE), "",
                                       sprintf("%.2f [%.2f, %.2f]",
                                               myTab2$beta, myTab2$lowerCI95,myTab2$upperCI95))
setorder(myTab2, variable, sex, flag2)

#' Plot per exposure type
myTab2[,Setting := paste0("    ",sex, " - ",flag2)]
dummy = data.table(Setting = c("mean","variability"))
plotData = rbind(myTab2,dummy, fill=T)
plotData = plotData[c(9,1:4,10,5:8)]

plotData[is.na(flag),` ` := ""]
plotData[is.na(flag),`Estimate \n[95% CI]`:= ""]
plotData[is.na(flag),condF := ""]

dummy = plotData$Setting
dummy[c(1,6)] = "white"
dummy[c(2,7)] = "#BDD3E5"
dummy[c(3,8)] = "#BCF1FF"
dummy[c(4,9)] = "#B7FFB7"
dummy[c(5,10)] = "#FFECC9"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Causal effect of the of TC on CAD per exposure type (gw. sig. SNPs)")
setnames(plotData,"condF","cond \nF-stat")
setnames(plotData,"Setting", "Exposure type \n sex - age group")

p2<- forest(plotData[,c(17,15,16,12)],
            est = plotData$beta,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/ForestPlots/3_GWSNPs_2sample_MenWomenAgeGroups_noSlope_perXType.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' Plot per sample group
myTab2[,Setting := variable]
myTab2[Setting=="var",Setting := "variability"]
dummy = data.table(Setting = c("older men","younger men","post-menopausal women","pre-menopausal women"))
plotData = rbind(myTab2,dummy, fill=T)
plotData = plotData[c(9,1,5,10,2,6,11,3,7,12,4,8)]

plotData[is.na(flag),` ` := ""]
plotData[is.na(flag),`Estimate \n[95% CI]`:= ""]
plotData[is.na(flag),condF := ""]

dummy = plotData$Setting
dummy[c(1,4,7,10)] = "white"
dummy[c(2,3)] = "#BDD3E5"
dummy[c(5,6)] = "#BCF1FF"
dummy[c(8,9)] = "#B7FFB7"
dummy[c(11,12)] = "#FFECC9"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Causal effect of the of TC on CAD per sex and age group (gw. sig. SNPs)")
setnames(plotData,"condF","cond \nF-stat")
setnames(plotData,"Setting", "Sex and age group \n exposure type")

p2<- forest(plotData[,c(17,15,16,12)],
            est = plotData$beta,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/ForestPlots/3_GWSNPs_2sample_MenWomenAgeGroups_noSlope_perGroup.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
