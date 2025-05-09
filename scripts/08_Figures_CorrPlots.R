#' ---
#' title: "Figures - genetic correlation"
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
#' **Correlation Plot between mean - slope - variability in the different sample groups**
#' 
#' - check 1: main model vs no slope model per subgroup (men, women, pre-, post-menopausal)
#' - check 2: men vs women, and young vs old 
#' - check 3: young and old men vs pre- and post-menopausal women (for the poster at the Wellcome Conference)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # Load data ####
#' ***
myFiles1 = list.files(path ="../results/",pattern="03_GAMLSS")
myFiles1 = myFiles1[grepl("coloc",myFiles1)]
myFiles1 = myFiles1[c(1:4,8,9)]
myFiles2 = list.files(path ="../results/",pattern="05_GAMLSS")
myFiles2 = myFiles2[grepl("coloc",myFiles2)]

myFiles = c(myFiles1,myFiles2)
myFiles
#' I want the model in the following format: sex_age_model, e.g. men_young_Slope
#' 
myTraits = c("men_combined_Slope","women_combined_Slope","women_old_Slope","women_young_Slope",
             "men_young_Slope","men_old_Slope","men_combined_noSlope","women_combined_noSlope",
             "women_old_noSlope","women_young_noSlope","men_young_noSlope","men_old_noSlope")

dumTab = foreach(i=1:length(myFiles))%do%{
  #i=1
  load(paste0("../results/",myFiles[i]))
  myAssocs_X[,model := myTraits[i]]
  myAssocs_X
}
myAssocs = rbindlist(dumTab,fill=T)
myAssocs[,table(model)]

#' # Get long & wide format ####
#' ***
data_long1 = melt(myAssocs,
                  id.vars=names(myAssocs)[1:3],
                  measure.vars=c("beta_mean", "beta_slope","beta_var"),
                  variable.name="type",
                  value.name="beta")

data_long2 = melt(myAssocs,
                  id.vars=names(myAssocs)[1:3],
                  measure.vars=c("SE_mean", "SE_slope","SE_var"),
                  variable.name="type",
                  value.name="SE")

data_long3 = melt(myAssocs,
                  id.vars=names(myAssocs)[1:3],
                  measure.vars=c("tval_mean", "tval_slope","tval_var"),
                  variable.name="type",
                  value.name="tval")

data_long4 = melt(myAssocs,
                  id.vars=names(myAssocs)[1:3],
                  measure.vars=c("pval_mean", "pval_slope","pval_var"),
                  variable.name="type",
                  value.name="pval")

myAssocs_long = cbind(data_long1,data_long2[,5],data_long3[,5],data_long4[,5])
myAssocs_long[,type := gsub("beta_","",type)]
myAssocs_long = myAssocs_long[!is.na(beta),]

myAssocs_long[,dumID := paste(type,model,sep="_")]
myAssocs_wide_beta = dcast(myAssocs_long,SNP ~ dumID,value.var = c("beta"))
data_wide_matrix = as.matrix(myAssocs_wide_beta[,-1])

CorTab = cor(data_wide_matrix,use = "pairwise.complete.obs")
testRes = cor.mtest(data_wide_matrix, conf.level = 0.95)

save(myAssocs_long, file = "../temp/CorPlot_Associations_longFormat.RData")

#' # Check 1: within sample groups
#' 
#' ## Men - combined
filt = grepl("_men_comb",colnames(CorTab))
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_men_combined_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",5)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/Men_combined_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Men - young
filt = grepl("_men_young",colnames(CorTab))
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_men_young_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",5)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/Men_young_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Men - old
filt = grepl("_men_old",colnames(CorTab))
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_men_old_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",5)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/Men_young_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Women - combined
filt = grepl("_women_combined",colnames(CorTab))
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_women_combined_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",5)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/Women_combined_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Women - young
filt = grepl("_women_young",colnames(CorTab))
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_women_young_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",5)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/Women_young_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Women - old
filt = grepl("_women_old",colnames(CorTab))
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_women_old_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",5)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/Women_old_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' # Check 2: between sample groups
#' ***
#' 
#' ## Men vs women - combined
filt1 = grepl("_men_comb",colnames(CorTab))
filt2 = grepl("_women_comb",colnames(CorTab))
filt = filt1 | filt2
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_combined_"," ",rownames(CorTab2))
rownames(CorTab2) = gsub("_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",10)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/MenWomen_combined_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Men vs women - young
filt1 = grepl("_men_young",colnames(CorTab))
filt2 = grepl("_women_young",colnames(CorTab))
filt = filt1 | filt2
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_young_"," ",rownames(CorTab2))
rownames(CorTab2) = gsub("_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",10)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/MenWomen_young_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Men vs women - old
filt1 = grepl("_men_old",colnames(CorTab))
filt2 = grepl("_women_old",colnames(CorTab))
filt = filt1 | filt2
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_old_"," ",rownames(CorTab2))
rownames(CorTab2) = gsub("_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",10)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/MenWomen_old_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Young vs old - men
filt1 = grepl("_men_young",colnames(CorTab))
filt2 = grepl("_men_old",colnames(CorTab))
filt = filt1 | filt2
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_men_"," ",rownames(CorTab2))
rownames(CorTab2) = gsub("_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",10)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/YoungOld_men_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' ## Young vs old - women
filt1 = grepl("_women_young",colnames(CorTab))
filt2 = grepl("_women_old",colnames(CorTab))
filt = filt1 | filt2
CorTab2 = CorTab[filt,filt]
rownames(CorTab2) = gsub("_women_"," ",rownames(CorTab2))
rownames(CorTab2) = gsub("_"," ",rownames(CorTab2))
colnames(CorTab2) = rep("",10)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/Correlation/YoungOld_women_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#' # Check 3: poster
#' ***
CorTab3 = cor(data_wide_matrix[,c(4,22,6,24,10,28,12,30)],use = "pairwise.complete.obs")
testRes = cor.mtest(data_wide_matrix[,c(4,22,6,24,10,28,12,30)], conf.level = 0.95)
rownames(CorTab3) = gsub("_noSlope","",rownames(CorTab3))
rownames(CorTab3) = gsub("_"," ",rownames(CorTab3))
colnames(CorTab3) = rep("",8)

corrplot(CorTab3,p.mat = testRes$p,number.cex = 1,number.digits = 2,sig.level = 0.005,
         type = "upper", tl.col = "black", tl.srt = 45,addCoef.col = 'white',na.label = " ")

png(file=paste0("../results/_figures/Correlation/YoungOld_MenWomen_Corrplot.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab3,p.mat = testRes$p,number.cex = 1,number.digits = 2,sig.level = 0.005,
         type = "upper", tl.col = "black", tl.srt = 45,addCoef.col = 'white',na.label = " ")
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")