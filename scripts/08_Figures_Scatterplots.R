#' ---
#' title: "Figures - scatter plots"
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
#' **Scatter Plot between mean or variability effects in the different sample groups or models**
#' 
#' - check 1: main model vs no slope model per subgroup (men, women, pre-, post-menopausal)
#' - check 2: men vs women, and young vs old 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # Load data ####
#' ***
load("../temp/CorPlot_Associations_longFormat.RData")
myAssocs_wide_beta = dcast(myAssocs_long,SNP ~ dumID,value.var = c("beta"))
data_wide_matrix = as.matrix(myAssocs_wide_beta[,-1])

#' # Check 1: within sample groups
#' 
#' ## Men - combined
plotData = rbind(myAssocs_wide_beta[,c(1,2,3)],myAssocs_wide_beta[,c(1,20,21)],use.names=F)
names(plotData)[2] = "main"
names(plotData)[3] = "noSlope"
plotData[,type := "mean"]
plotData[duplicated(SNP),type := "variability"]
plotData[,table(type)]

myPlot1 = ggplot(plotData, aes(x=main, y=noSlope, color=type)) +
  facet_wrap(~type, nrow = 1,scales = "free",
             labeller = as_labeller(c("mean" = "SNP effect on mean", 
                                      "variability" = "SNP effect on variability") ) )+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")+
  labs(x="SNP effect in main model",y = "SNP effect in no-slope model") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/Men_combined_Scatterplot.png"),
    width=1500,height=900,res = 200)
myPlot1
dev.off()

#' ## Men - young
plotData = rbind(myAssocs_wide_beta[,c(1,6,7)],myAssocs_wide_beta[,c(1,24,25)],use.names=F)
names(plotData)[2] = "main"
names(plotData)[3] = "noSlope"
plotData[,type := "mean"]
plotData[duplicated(SNP),type := "variability"]
plotData[,table(type)]

myPlot1 = ggplot(plotData, aes(x=main, y=noSlope, color=type)) +
  facet_wrap(~type, nrow = 1,scales = "free",
             labeller = as_labeller(c("mean" = "SNP effect on mean", 
                                      "variability" = "SNP effect on variability") ) )+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")+
  labs(x="SNP effect in main model",y = "SNP effect in no-slope model") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/Men_young_Scatterplot.png"),
    width=1500,height=900,res = 200)
myPlot1
dev.off()

#' ## Men - old
plotData = rbind(myAssocs_wide_beta[,c(1,4,5)],myAssocs_wide_beta[,c(1,22,23)],use.names=F)
names(plotData)[2] = "main"
names(plotData)[3] = "noSlope"
plotData[,type := "mean"]
plotData[duplicated(SNP),type := "variability"]
plotData[,table(type)]

myPlot1 = ggplot(plotData, aes(x=main, y=noSlope, color=type)) +
  facet_wrap(~type, nrow = 1,scales = "free",
             labeller = as_labeller(c("mean" = "SNP effect on mean", 
                                      "variability" = "SNP effect on variability") ) )+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")+
  labs(x="SNP effect in main model",y = "SNP effect in no-slope model") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/Men_old_Scatterplot.png"),
    width=1500,height=900,res = 200)
myPlot1
dev.off()

#' ## Women - combined
plotData = rbind(myAssocs_wide_beta[,c(1,8,9)],myAssocs_wide_beta[,c(1,26,27)],use.names=F)
names(plotData)[2] = "main"
names(plotData)[3] = "noSlope"
plotData[,type := "mean"]
plotData[duplicated(SNP),type := "variability"]
plotData[,table(type)]

myPlot1 = ggplot(plotData, aes(x=main, y=noSlope, color=type)) +
  facet_wrap(~type, nrow = 1,scales = "free",
             labeller = as_labeller(c("mean" = "SNP effect on mean", 
                                      "variability" = "SNP effect on variability") ) )+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")+
  labs(x="SNP effect in main model",y = "SNP effect in no-slope model") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/Women_combined_Scatterplot.png"),
    width=1500,height=900,res = 200)
myPlot1
dev.off()

#' ## Women - young
plotData = rbind(myAssocs_wide_beta[,c(1,12,13)],myAssocs_wide_beta[,c(1,30,31)],use.names=F)
names(plotData)[2] = "main"
names(plotData)[3] = "noSlope"
plotData[,type := "mean"]
plotData[duplicated(SNP),type := "variability"]
plotData[,table(type)]

myPlot1 = ggplot(plotData, aes(x=main, y=noSlope, color=type)) +
  facet_wrap(~type, nrow = 1,scales = "free",
             labeller = as_labeller(c("mean" = "SNP effect on mean", 
                                      "variability" = "SNP effect on variability") ) )+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")+
  labs(x="SNP effect in main model",y = "SNP effect in no-slope model") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/Women_young_Scatterplot.png"),
    width=1500,height=900,res = 200)
myPlot1
dev.off()

#' ## Women - old
plotData = rbind(myAssocs_wide_beta[,c(1,10,11)],myAssocs_wide_beta[,c(1,28,29)],use.names=F)
names(plotData)[2] = "main"
names(plotData)[3] = "noSlope"
plotData[,type := "mean"]
plotData[duplicated(SNP),type := "variability"]
plotData[,table(type)]

myPlot1 = ggplot(plotData, aes(x=main, y=noSlope, color=type)) +
  facet_wrap(~type, nrow = 1,scales = "free",
             labeller = as_labeller(c("mean" = "SNP effect on mean", 
                                      "variability" = "SNP effect on variability") ) )+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_point(size=3)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")+
  labs(x="SNP effect in main model",y = "SNP effect in no-slope model") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/Women_old_Scatterplot.png"),
    width=1500,height=900,res = 200)
myPlot1
dev.off()

#' # Check 2: between sample groups
#' 
#' ## Men vs women - combined
plotData = cbind(myAssocs_long[model %in% c("men_combined_Slope","men_combined_noSlope"),c(1:8)],
                 myAssocs_long[model %in% c("women_combined_Slope","women_combined_noSlope"),c(5:8)])
names(plotData)[5:8] = paste0("men_",names(plotData)[5:8])
names(plotData)[9:12] = paste0("women_",names(plotData)[9:12])
plotData[,table(type,model)]
plotData[,model2:= gsub("men_combined","",model)]
plotData[model2=="_Slope",model2:= "main"]
plotData[model2=="_noSlope",model2:= "no slope"]

plotData[,type2 := paste0(type,"_",model2)]
plotData[,table(type2)]

myPlot1 = ggplot(plotData, aes(x=men_beta, y=women_beta, color=model2)) +
  facet_wrap(~type, nrow = 1,scales = "free")+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_errorbar(aes(ymin = women_beta- 1.96*women_SE, ymax = women_beta+ 1.96*women_SE),linewidth=1,alpha=0.5) +
  geom_errorbarh(aes(xmin = men_beta- 1.96*men_SE, xmax = men_beta + 1.96*men_SE),linewidth=1,alpha=0.5) +
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        #legend.position = "none"
  )+
  labs(x="SNP effect in men",y = "SNP effect in women",color = "GAMLSS \nmodel") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/MenWomen_combined_Scatterplot.png"),
    width=4000,height=2000,res = 200)
myPlot1
dev.off()

#' ## Men vs women - young
plotData = cbind(myAssocs_long[model %in% c("men_young_Slope","men_young_noSlope"),c(1:8)],
                 myAssocs_long[model %in% c("women_young_Slope","women_young_noSlope"),c(5:8)])
names(plotData)[5:8] = paste0("men_",names(plotData)[5:8])
names(plotData)[9:12] = paste0("women_",names(plotData)[9:12])
plotData[,table(type,model)]
plotData[,model2:= gsub("men_young","",model)]
plotData[model2=="_Slope",model2:= "main"]
plotData[model2=="_noSlope",model2:= "no slope"]

plotData[,type2 := paste0(type,"_",model2)]
plotData[,table(type2)]

myPlot1 = ggplot(plotData, aes(x=men_beta, y=women_beta, color=model2)) +
  facet_wrap(~type, nrow = 1,scales = "free")+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_errorbar(aes(ymin = women_beta- 1.96*women_SE, ymax = women_beta+ 1.96*women_SE),linewidth=1,alpha=0.5) +
  geom_errorbarh(aes(xmin = men_beta- 1.96*men_SE, xmax = men_beta + 1.96*men_SE),linewidth=1,alpha=0.5) +
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        #legend.position = "none"
  )+
  labs(x="SNP effect in men",y = "SNP effect in women",color = "GAMLSS \nmodel") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/MenWomen_young_Scatterplot.png"),
    width=4000,height=2000,res = 200)
myPlot1
dev.off()

#' ## Men vs women - old
plotData = cbind(myAssocs_long[model %in% c("men_old_Slope","men_old_noSlope"),c(1:8)],
                 myAssocs_long[model %in% c("women_old_Slope","women_old_noSlope"),c(5:8)])
names(plotData)[5:8] = paste0("men_",names(plotData)[5:8])
names(plotData)[9:12] = paste0("women_",names(plotData)[9:12])
plotData[,table(type,model)]
plotData[,model2:= gsub("men_old","",model)]
plotData[model2=="_Slope",model2:= "main"]
plotData[model2=="_noSlope",model2:= "no slope"]

plotData[,type2 := paste0(type,"_",model2)]
plotData[,table(type2)]

myPlot1 = ggplot(plotData, aes(x=men_beta, y=women_beta, color=model2)) +
  facet_wrap(~type, nrow = 1,scales = "free")+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_errorbar(aes(ymin = women_beta- 1.96*women_SE, ymax = women_beta+ 1.96*women_SE),linewidth=1,alpha=0.5) +
  geom_errorbarh(aes(xmin = men_beta- 1.96*men_SE, xmax = men_beta + 1.96*men_SE),linewidth=1,alpha=0.5) +
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        #legend.position = "none"
  )+
  labs(x="SNP effect in men",y = "SNP effect in women",color = "GAMLSS \nmodel") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/MenWomen_old_Scatterplot.png"),
    width=4000,height=2000,res = 200)
myPlot1
dev.off()

#' ## Young vs old - men
plotData = cbind(myAssocs_long[model %in% c("men_young_Slope","men_young_noSlope"),c(1:8)],
                 myAssocs_long[model %in% c("men_old_Slope","men_old_noSlope"),c(5:8)])
names(plotData)[5:8] = paste0("young_",names(plotData)[5:8])
names(plotData)[9:12] = paste0("old_",names(plotData)[9:12])
plotData[,table(type,model)]
plotData[,model2:= gsub("men_young","",model)]
plotData[model2=="_Slope",model2:= "main"]
plotData[model2=="_noSlope",model2:= "no slope"]

plotData[,type2 := paste0(type,"_",model2)]
plotData[,table(type2)]

myPlot1 = ggplot(plotData, aes(x=young_beta, y=old_beta, color=model2)) +
  facet_wrap(~type, nrow = 1,scales = "free")+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_errorbar(aes(ymin = old_beta- 1.96*old_SE, ymax = old_beta+ 1.96*old_SE),linewidth=1,alpha=0.5) +
  geom_errorbarh(aes(xmin = young_beta- 1.96*young_SE, xmax = young_beta + 1.96*young_SE),linewidth=1,alpha=0.5) +
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        #legend.position = "none"
  )+
  labs(x="SNP effect in young men",y = "SNP effect in old men",color = "GAMLSS \nmodel") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/YoungOld_men_Scatterplot.png"),
    width=4000,height=2000,res = 200)
myPlot1
dev.off()

#' ## Young vs old - women
plotData = cbind(myAssocs_long[model %in% c("women_young_Slope","women_young_noSlope"),c(1:8)],
                 myAssocs_long[model %in% c("women_old_Slope","women_old_noSlope"),c(5:8)])
names(plotData)[5:8] = paste0("young_",names(plotData)[5:8])
names(plotData)[9:12] = paste0("old_",names(plotData)[9:12])
plotData[,table(type,model)]
plotData[,model2:= gsub("women_young","",model)]
plotData[model2=="_Slope",model2:= "main"]
plotData[model2=="_noSlope",model2:= "no slope"]

plotData[,type2 := paste0(type,"_",model2)]
plotData[,table(type2)]

myPlot1 = ggplot(plotData, aes(x=young_beta, y=old_beta, color=model2)) +
  facet_wrap(~type, nrow = 1,scales = "free")+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", linewidth=1.15)+
  geom_abline(intercept = 0,slope=1,color="grey", linetype="dashed", linewidth=1.15)+
  geom_errorbar(aes(ymin = old_beta- 1.96*old_SE, ymax = old_beta+ 1.96*old_SE),linewidth=1,alpha=0.5) +
  geom_errorbarh(aes(xmin = young_beta- 1.96*young_SE, xmax = young_beta + 1.96*young_SE),linewidth=1,alpha=0.5) +
  geom_point(size=2.5)+ 
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 14, face = "bold"),
        #legend.position = "none"
  )+
  labs(x="SNP effect in young women",y = "SNP effect in old women",color = "GAMLSS \nmodel") 

myPlot1
png(file=paste0("../results/_figures/Scatterplots/YoungOld_men_Scatterplot.png"),
    width=4000,height=2000,res = 200)
myPlot1
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
