#' ---
#' title: "Figures - forest plots 4"
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
myMVMR_results = list.files(path = "../results/",pattern = "07_MVMR")

dumTab4 = foreach(i = 1:length(myMVMR_results))%do%{
  #i=1
  load(paste0("../results/",myMVMR_results[i]))
  names(MVMR_results)
  MVMR_results = MVMR_results[setting == "multivariate",]
  MVMR_results[outcome == "UKB", setting := "1-sample"]
  MVMR_results[outcome != "UKB", setting := "2-sample"]
  MVMR_results[outcome != "UKB", outcome := "Aragam_CAD"]
  MVMR_results[outcome == "UKB", outcome := "UKB_CAD"]
  
  # MVMR_tab_wide = dcast(MVMR_results, ID + setting + exposure + outcome + threshold + NR_SNPs_total + HeteroStat + HeteroStat_pval ~ exposure_type, 
  #                       value.var=c("NR_SNPs_type","beta_IVW","SE_IVW","pval_IVW","condFstat"))
  # names(MVMR_tab_wide) = gsub("IVW_","",names(MVMR_tab_wide))
  # names(MVMR_tab_wide) = gsub("NR_SNPs_type","SNPs",names(MVMR_tab_wide))
  # names(MVMR_tab_wide) = gsub("TC_women_","",names(MVMR_tab_wide))
  # names(MVMR_tab_wide) = gsub("TC_men_","",names(MVMR_tab_wide))
  # setnames(MVMR_tab_wide,"ID","flag")
  # if(i==1)MVMR_tab_wide[,exposure := "TC_women"]
  # if(i==2)MVMR_tab_wide[,exposure := "TC_men"]
  # 
  # MVMR_tab_wide
  MVMR_results
}
tab3 = rbindlist(dumTab4,fill=T)
names(tab3)

#' # Plot 1: all SNPs ####
#' ***
#' All SNP approach, 2-sample MVMR, all 12 settings
#' 
myTab2 = copy(tab3)
myTab2 = myTab2[threshold == "all_SNPs",]
myTab2 = myTab2[setting == "2-sample",]

names(myTab2)[7:10] = c("beta","SE","pval","condF")
myTab2 = myTab2[!is.na(beta),]

myTab2[,lowerCI95 := beta-1.96*SE]
myTab2[,upperCI95 := beta+1.96*SE]

myTab2[,condF := round(condF,2)]
myTab2[,condF := as.character(condF)]

myTab2$` ` <- paste(rep(" ", 50), collapse = " ")
myTab2$`Estimate \n[95% CI]` <- ifelse(is.na(myTab2$SE), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 myTab2$beta, myTab2$lowerCI95,myTab2$upperCI95))
setorder(myTab2, exposure_type)

#' Plot per exposure type
myTab2[,Setting := gsub("TC_women_","    ",exposure_type)]
myTab2[,Setting := gsub("TC_men_","    ",Setting)]
myTab2[,Setting := gsub(".*_","",Setting)]
myTab2[Setting == "var",Setting := "variability"]
myTab2[,Setting := paste("    ",Setting)]

dummy = data.table(Setting = c("men","women"))
plotData = rbind(myTab2,dummy, fill=T)
plotData = plotData[c(9,1:4,10,5:8)]

plotData[is.na(outcome),` ` := ""]
plotData[is.na(outcome),`Estimate \n[95% CI]`:= ""]
plotData[is.na(outcome),condF := ""]

dummy = plotData$Setting
dummy[c(1,6)] = "white"
dummy[c(2,3)] = "#BDD3E5"
dummy[c(4,5)] = "#BCF1FF"
dummy[c(7,8)] = "#B7FFB7"
dummy[c(9,10)] = "#FFECC9"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Causal effect of the of TC on CAD (all candidate SNPs)")
setnames(plotData,"condF","cond \nF-stat")
setnames(plotData,"Setting", "Sex \n type")

p2<- forest(plotData[,c(19,17,18,10)],
            est = plotData$beta,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/ForestPlots/4_AllSNPs_2sample_AgeGroups_noSlope.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # Plot 2: genome-wide significant SNPs ####
#' ***
#' Genome-wide significant SNPs approach, 2-sample MVMR, all 12 settings
#' 
myTab2 = copy(tab3)
myTab2 = myTab2[threshold != "all_SNPs",]
myTab2 = myTab2[setting == "2-sample",]

names(myTab2)[7:10] = c("beta","SE","pval","condF")
myTab2 = myTab2[!is.na(beta),]

myTab2[,lowerCI95 := beta-1.96*SE]
myTab2[,upperCI95 := beta+1.96*SE]

myTab2[,condF := round(condF,2)]
myTab2[,condF := as.character(condF)]

myTab2$` ` <- paste(rep(" ", 50), collapse = " ")
myTab2$`Estimate \n[95% CI]` <- ifelse(is.na(myTab2$SE), "",
                                       sprintf("%.2f [%.2f, %.2f]",
                                               myTab2$beta, myTab2$lowerCI95,myTab2$upperCI95))
setorder(myTab2, exposure_type)

#' Plot per exposure type
myTab2[,Setting := gsub("TC_women_","    ",exposure_type)]
myTab2[,Setting := gsub("TC_men_","    ",Setting)]
myTab2[,Setting := gsub(".*_","",Setting)]
myTab2[Setting == "var",Setting := "variability"]
myTab2[,Setting := paste("    ",Setting)]

dummy = data.table(Setting = c("men","women"))
plotData = rbind(myTab2,dummy, fill=T)
plotData = plotData[c(9,1:4,10,5:8)]

plotData[is.na(outcome),` ` := ""]
plotData[is.na(outcome),`Estimate \n[95% CI]`:= ""]
plotData[is.na(outcome),condF := ""]

dummy = plotData$Setting
dummy[c(1,6)] = "white"
dummy[c(2,3)] = "#BDD3E5"
dummy[c(4,5)] = "#BCF1FF"
dummy[c(7,8)] = "#B7FFB7"
dummy[c(9,10)] = "#FFECC9"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))

myXlab = paste0("Causal effect of the of TC on CAD (gw. sig. SNPs)")
setnames(plotData,"condF","cond \nF-stat")
setnames(plotData,"Setting", "Sex \n type")

p2<- forest(plotData[,c(19,17,18,10)],
            est = plotData$beta,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            title = myXlab,
            theme = tm1)

plot(p2)

filename = paste0("../results/_figures/ForestPlots/4_GWSNPs_2sample_AgeGroups_noSlope.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(p2)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
