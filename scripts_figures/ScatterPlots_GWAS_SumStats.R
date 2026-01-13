load("../results/UKB/05_potentialIVs_trajGWAS_Wald.RData")
dataX_UKB_trajGWAS = copy(exposureData)
load("../results/UKB/06_potentialIVs_GAMLSS_bugfixed.RData")
dataX_UKB_GAMLSS = copy(myAssocs)

load("../results/INTERVAL/03_potentialIVs_trajGWAS_Wald.RData")
dataX_INTERVAL_trajGWAS = copy(exposureData)
load("../results/INTERVAL/04_potentialIVs_GAMLSS.RData")
dataX_INTERVAL_GAMLSS = copy(myAssocs)

setorder(dataX_UKB_GAMLSS,setting,chr,pos_b37)
setorder(dataX_UKB_trajGWAS,setting,chr,pos_b37)
table(dataX_UKB_GAMLSS$rsID == dataX_UKB_trajGWAS$rsID)
table(dataX_UKB_GAMLSS$setting == dataX_UKB_trajGWAS$setting)

setorder(dataX_INTERVAL_GAMLSS,setting,chr,pos_b37)
setorder(dataX_INTERVAL_trajGWAS,setting,chr,pos_b37)
table(dataX_INTERVAL_GAMLSS$rsID == dataX_INTERVAL_trajGWAS$rsID)
table(dataX_INTERVAL_GAMLSS$setting == dataX_INTERVAL_trajGWAS$setting)

mySettings = unique(dataX_trajGWAS$setting)
mySettings2 = c("men","old men","post-menopausal women","pre-menopausal women","women","young men")

library(ggplot2)
library(ggpubr)

# Scatter plot per study, setting and exposure type for the effect estimates

# same type and study, compare method
for(i in 1:6){
  #i=1
  data1 = copy(dataX_UKB_trajGWAS)
  data1 = data1[setting == mySettings[i]]
  data2 = copy(dataX_UKB_GAMLSS)
  data2 = data2[setting == mySettings[i]]
  stopifnot(data1$rsID == data2$rsID)
  data3 = copy(dataX_INTERVAL_trajGWAS)
  data3 = data3[setting == mySettings[i]]
  data4 = copy(dataX_INTERVAL_GAMLSS)
  data4 = data4[setting == mySettings[i]]
  stopifnot(data3$rsID == data4$rsID)
  
  plotData1 = data.table(snp = data1$rsID, 
                         est1 = data1$mean_beta, 
                         est2 = data2$mean_beta, 
                         est1_se = data1$mean_SE, 
                         est2_se = data2$mean_SE, 
                         type = "mean", 
                         study = "UKB")
  plotData2 = data.table(snp = data1$rsID, 
                         est1 = data1$var_beta, 
                         est2 = data2$var_beta, 
                         est1_se = data1$var_SE, 
                         est2_se = data2$var_SE, 
                         type = "variability", 
                         study = "UKB")
  plotData3 = data.table(snp = data3$rsID, 
                         est1 = data3$mean_beta, 
                         est2 = data4$mean_beta, 
                         est1_se = data3$mean_SE, 
                         est2_se = data4$mean_SE, 
                         type = "mean", 
                         study = "INTERVAL")
  plotData4 = data.table(snp = data3$rsID, 
                         est1 = data3$var_beta, 
                         est2 = data4$var_beta, 
                         est1_se = data3$var_SE, 
                         est2_se = data4$var_SE, 
                         type = "variability", 
                         study = "INTERVAL")
  
  plotData = rbind(plotData1,plotData2,plotData3,plotData4)
  plotData[,categ := paste(study, type,sep = " - ")]
  plotData[,Zscore1 := est1/est1_se]
  plotData[,Zscore2 := est2/est2_se]
  
  p1 <- ggplot(plotData,aes(x=est1,y=est2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab("Estimate using trajGWAS") + 
    ylab("Estimate using GAMLSS")+
    ggtitle(paste0("Scatter plot of SNP estimators in ",mySettings2[i]))+
    theme_bw()
  
  print(p1)
  
  p2 <- ggplot(plotData,aes(x=Zscore1,y=Zscore2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab("Estimate using trajGWAS") + 
    ylab("Estimate using GAMLSS")+
    ggtitle(paste0("Scatter plot of SNP Z-scores in ",mySettings2[i]))+
    theme_bw()
  
  print(p2)
  
  filename1 = paste0("../results/_figures/ScatterPlots/TrajGWAS_GAMLSS_",mySettings[i],"_estiamtes.png")
  png(filename = filename1,width = 1800, height = 1600, res=200)
  print(p1)
  dev.off()
  
  filename2 = paste0("../results/_figures/ScatterPlots/TrajGWAS_GAMLSS_",mySettings[i],"_Zscores.png")
  png(filename = filename2,width = 1800, height = 1600, res=200)
  print(p2)
  dev.off()
  
  print(summary(lm(est1 ~ est2 + study, data = plotData,sub = type=="mean")))
  print(summary(lm(est1 ~ est2 + study, data = plotData,sub = type=="variability")))
}

# same method and type, compare study
for(i in 1:6){
  #i=1
  data1 = copy(dataX_UKB_trajGWAS)
  data1 = data1[setting == mySettings[i]]
  data2 = copy(dataX_UKB_GAMLSS)
  data2 = data2[setting == mySettings[i]]
  stopifnot(data1$rsID == data2$rsID)
  data3 = copy(dataX_INTERVAL_trajGWAS)
  data3 = data3[setting == mySettings[i]]
  data4 = copy(dataX_INTERVAL_GAMLSS)
  data4 = data4[setting == mySettings[i]]
  stopifnot(data3$rsID == data4$rsID)
  
  data1 = data1[rsID %in% data3$rsID]
  data2 = data2[rsID %in% data4$rsID]
  data3 = data3[rsID %in% data1$rsID]
  data4 = data4[rsID %in% data2$rsID]
  
  plotData1 = data.table(snp = data1$rsID, 
                         est1 = data1$mean_beta, 
                         est2 = data3$mean_beta, 
                         est1_se = data1$mean_SE, 
                         est2_se = data3$mean_SE, 
                         type = "mean", 
                         method = "trajGWAS")
  plotData2 = data.table(snp = data1$rsID, 
                         est1 = data1$var_beta, 
                         est2 = data3$var_beta, 
                         est1_se = data1$var_SE, 
                         est2_se = data3$var_SE, 
                         type = "variability", 
                         method = "trajGWAS")
  plotData3 = data.table(snp = data2$rsID, 
                         est1 = data2$mean_beta, 
                         est2 = data4$mean_beta, 
                         est1_se = data2$mean_SE, 
                         est2_se = data4$mean_SE, 
                         type = "mean", 
                         method = "GAMLSS")
  plotData4 = data.table(snp = data2$rsID, 
                         est1 = data2$var_beta, 
                         est2 = data4$var_beta, 
                         est1_se = data2$var_SE, 
                         est2_se = data4$var_SE, 
                         type = "variability", 
                         method = "GAMLSS")
  
  plotData = rbind(plotData1,plotData2,plotData3,plotData4)
  plotData[,categ := paste(method, type,sep = " - ")]
  plotData[,Zscore1 := est1/est1_se]
  plotData[,Zscore2 := est2/est2_se]
  
  p1 <- ggplot(plotData,aes(x=est1,y=est2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab("Estimate using UKB data") + 
    ylab("Estimate using INTERVAL data")+
    ggtitle(paste0("Scatter plot of SNP estimators in ",mySettings2[i]))+
    theme_bw()
  
  print(p1)
  
  p2 <- ggplot(plotData,aes(x=Zscore1,y=Zscore2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab("Estimate using UKB data") + 
    ylab("Estimate using INTERVAL data")+
    ggtitle(paste0("Scatter plot of SNP Z-scores in ",mySettings2[i]))+
    theme_bw()
  
  print(p2)
  
  filename1 = paste0("../results/_figures/ScatterPlots/UKB_INTERVAL_",mySettings[i],"_estiamtes.png")
  png(filename = filename1,width = 1800, height = 1600, res=200)
  print(p1)
  dev.off()
  
  filename2 = paste0("../results/_figures/ScatterPlots/UKB_INTERVAL_",mySettings[i],"_Zscores.png")
  png(filename = filename2,width = 1800, height = 1600, res=200)
  print(p2)
  dev.off()
  
  print(summary(lm(est1 ~ est2 + method, data = plotData,sub = type=="mean")))
  print(summary(lm(est1 ~ est2 + method, data = plotData,sub = type=="variability")))
}

# same method and study, compare type
for(i in 1:6){
  #i=1
  data1 = copy(dataX_UKB_trajGWAS)
  data1 = data1[setting == mySettings[i]]
  data2 = copy(dataX_UKB_GAMLSS)
  data2 = data2[setting == mySettings[i]]
  stopifnot(data1$rsID == data2$rsID)
  data3 = copy(dataX_INTERVAL_trajGWAS)
  data3 = data3[setting == mySettings[i]]
  data4 = copy(dataX_INTERVAL_GAMLSS)
  data4 = data4[setting == mySettings[i]]
  stopifnot(data3$rsID == data4$rsID)
  
  plotData1 = data.table(snp = data1$rsID, 
                         est1 = data1$mean_beta, 
                         est2 = data1$var_beta, 
                         est1_se = data1$mean_SE, 
                         est2_se = data1$var_SE, 
                         study = "UKB", 
                         method = "trajGWAS")
  plotData2 = data.table(snp = data2$rsID, 
                         est1 = data2$mean_beta, 
                         est2 = data2$var_beta, 
                         est1_se = data2$mean_SE, 
                         est2_se = data2$var_SE, 
                         study = "UKB", 
                         method = "GAMLSS")
  plotData3 = data.table(snp = data3$rsID, 
                         est1 = data3$mean_beta, 
                         est2 = data3$var_beta, 
                         est1_se = data3$mean_SE, 
                         est2_se = data3$var_SE, 
                         study = "INTERVAL", 
                         method = "trajGWAS")
  plotData4 = data.table(snp = data4$rsID, 
                         est1 = data4$mean_beta, 
                         est2 = data4$var_beta, 
                         est1_se = data4$mean_SE, 
                         est2_se = data4$var_SE, 
                         study = "INTERVAL", 
                         method = "GAMLSS")
  
  plotData = rbind(plotData1,plotData2,plotData3,plotData4)
  plotData[,categ := paste(study, method,sep = " - ")]
  plotData[,Zscore1 := est1/est1_se]
  plotData[,Zscore2 := est2/est2_se]
  
  p1 <- ggplot(plotData,aes(x=est1,y=est2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab("Estimate for mean") + 
    ylab("Estimate for variability")+
    ggtitle(paste0("Scatter plot of SNP estimators in ",mySettings2[i]))+
    theme_bw()
  
  print(p1)
  
  p2 <- ggplot(plotData,aes(x=Zscore1,y=Zscore2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab("Estimate for mean") + 
    ylab("Estimate for variability")+
    ggtitle(paste0("Scatter plot of SNP Z-scores in ",mySettings2[i]))+
    theme_bw()
  
  print(p2)
  
  filename1 = paste0("../results/_figures/ScatterPlots/mean_variability_",mySettings[i],"_estiamtes.png")
  png(filename = filename1,width = 1800, height = 1600, res=200)
  print(p1)
  dev.off()
  
  filename2 = paste0("../results/_figures/ScatterPlots/mean_variability_",mySettings[i],"_Zscores.png")
  png(filename = filename2,width = 1800, height = 1600, res=200)
  print(p2)
  dev.off()
  
}

# Compare men and women
mySettings1 = c("men","old","young","old","post")
mySettings2 = c("women","post","pre","young","pre")
mySettings3 = c("men","old men","young men","old men", "post-menopausal women")
mySettings4 = c("women","post-menopausal women","pre-menopausal women","young men", "pre-menopausal women")


for(i in 1:5){
  #i=1
  data1 = copy(dataX_UKB_trajGWAS)
  data1 = data1[setting == mySettings1[i]]
  data2 = copy(dataX_UKB_trajGWAS)
  data2 = data2[setting == mySettings2[i]]
  stopifnot(data1$rsID == data2$rsID)
  data3 = copy(dataX_INTERVAL_trajGWAS)
  data3 = data3[setting == mySettings1[i]]
  data4 = copy(dataX_INTERVAL_trajGWAS)
  data4 = data4[setting == mySettings2[i]]
  stopifnot(data3$rsID == data4$rsID)
  
  plotData1 = data.table(snp = data1$rsID, 
                         est1 = data1$mean_beta, 
                         est2 = data2$mean_beta, 
                         est1_se = data1$mean_SE, 
                         est2_se = data2$mean_SE, 
                         study = "UKB", 
                         method = "trajGWAS",
                         type = "mean")
  plotData2 = data.table(snp = data1$rsID, 
                         est1 = data1$var_beta, 
                         est2 = data2$var_beta, 
                         est1_se = data1$var_SE, 
                         est2_se = data2$var_SE, 
                         study = "UKB", 
                         method = "trajGWAS",
                         type = "variability")
  plotData3 = data.table(snp = data3$rsID, 
                         est1 = data3$mean_beta, 
                         est2 = data4$mean_beta, 
                         est1_se = data3$mean_SE, 
                         est2_se = data4$mean_SE, 
                         study = "INTERVAL", 
                         method = "trajGWAS",
                         type = "mean")
  plotData4 = data.table(snp = data3$rsID, 
                         est1 = data3$var_beta, 
                         est2 = data4$var_beta, 
                         est1_se = data3$var_SE, 
                         est2_se = data4$var_SE, 
                         study = "INTERVAL", 
                         method = "trajGWAS",
                         type = "variability")
  

  plotData = rbind(plotData1,plotData2, plotData3,plotData4)
  plotData[,categ := paste(study, method,type, sep = " - ")]
  plotData[,Zscore1 := est1/est1_se]
  plotData[,Zscore2 := est2/est2_se]
  
  p1 <- ggplot(plotData,aes(x=est1,y=est2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab(paste0("Estimate in ",mySettings3[i])) + 
    ylab(paste0("Estimate in ",mySettings4[i])) + 
    ggtitle(paste0("Scatter plot of SNP estimators")) +
    theme_bw()
  
  print(p1)
  
  p2 <- ggplot(plotData,aes(x=Zscore1,y=Zscore2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab(paste0("Estimate in ",mySettings3[i])) + 
    ylab(paste0("Estimate in ",mySettings4[i])) + 
    ggtitle(paste0("Scatter plot of SNP Z-scores")) +
    theme_bw()
  
  print(p2)
  
  filename1 = paste0("../results/_figures/ScatterPlots/",mySettings1[i],"_",mySettings2[i],"_trajGWAS_estiamtes.png")
  png(filename = filename1,width = 1800, height = 1600, res=200)
  print(p1)
  dev.off()
  
  filename2 = paste0("../results/_figures/ScatterPlots/",mySettings1[i],"_",mySettings2[i],"_trajGWAS_Zscores.png")
  png(filename = filename2,width = 1800, height = 1600, res=200)
  print(p2)
  dev.off()
  
}

for(i in 1:5){
  #i=1
  data1 = copy(dataX_UKB_GAMLSS)
  data1 = data1[setting == mySettings1[i]]
  data2 = copy(dataX_UKB_GAMLSS)
  data2 = data2[setting == mySettings2[i]]
  stopifnot(data1$rsID == data2$rsID)
  data3 = copy(dataX_INTERVAL_GAMLSS)
  data3 = data3[setting == mySettings1[i]]
  data4 = copy(dataX_INTERVAL_GAMLSS)
  data4 = data4[setting == mySettings2[i]]
  stopifnot(data3$rsID == data4$rsID)
  
  plotData1 = data.table(snp = data1$rsID, 
                         est1 = data1$mean_beta, 
                         est2 = data2$mean_beta, 
                         est1_se = data1$mean_SE, 
                         est2_se = data2$mean_SE, 
                         study = "UKB", 
                         method = "GAMLSS",
                         type = "mean")
  plotData2 = data.table(snp = data1$rsID, 
                         est1 = data1$var_beta, 
                         est2 = data2$var_beta, 
                         est1_se = data1$var_SE, 
                         est2_se = data2$var_SE, 
                         study = "UKB", 
                         method = "GAMLSS",
                         type = "variability")
  plotData3 = data.table(snp = data3$rsID, 
                         est1 = data3$mean_beta, 
                         est2 = data4$mean_beta, 
                         est1_se = data3$mean_SE, 
                         est2_se = data4$mean_SE, 
                         study = "INTERVAL", 
                         method = "GAMLSS",
                         type = "mean")
  plotData4 = data.table(snp = data3$rsID, 
                         est1 = data3$var_beta, 
                         est2 = data4$var_beta, 
                         est1_se = data3$var_SE, 
                         est2_se = data4$var_SE, 
                         study = "INTERVAL", 
                         method = "GAMLSS",
                         type = "variability")
  
  
  plotData = rbind(plotData1,plotData2, plotData3,plotData4)
  plotData[,categ := paste(study, method,type, sep = " - ")]
  plotData[,Zscore1 := est1/est1_se]
  plotData[,Zscore2 := est2/est2_se]
  
  p1 <- ggplot(plotData,aes(x=est1,y=est2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab(paste0("Estimate in ",mySettings3[i])) + 
    ylab(paste0("Estimate in ",mySettings4[i])) + 
    ggtitle(paste0("Scatter plot of SNP estimators")) +
    theme_bw()
  
  print(p1)
  
  p2 <- ggplot(plotData,aes(x=Zscore1,y=Zscore2)) + 
    facet_wrap(~ categ,scales = "free") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth=1)+
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey", linewidth=1) + 
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) +
    stat_cor(method = "pearson") + 
    xlab(paste0("Estimate in ",mySettings3[i])) + 
    ylab(paste0("Estimate in ",mySettings4[i])) + 
    ggtitle(paste0("Scatter plot of SNP Z-scores")) +
    theme_bw()
  
  print(p2)
  
  filename1 = paste0("../results/_figures/ScatterPlots/",mySettings1[i],"_",mySettings2[i],"_GAMLSS_estiamtes.png")
  png(filename = filename1,width = 1800, height = 1600, res=200)
  print(p1)
  dev.off()
  
  filename2 = paste0("../results/_figures/ScatterPlots/",mySettings1[i],"_",mySettings2[i],"_GAMLSS_Zscores.png")
  png(filename = filename2,width = 1800, height = 1600, res=200)
  print(p2)
  dev.off()
  
}
