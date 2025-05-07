#' ---
#' title: "Figures - boxplots"
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
#' Get all relevant trajectories of TC levels in: 
#' 
#' - all samples
#' - men only
#' - women only
#' - pre-menopausal women only
#' - post-menopausal women only
#' 
#' Data is either the full data set (n=73,778) or the statin-free subset (n=51,190). 
#' 
#' For the poster at the Wellcome Conference (Women's Health: Genes, Data and Advancing Approaches, 27/01/2025 – 29/01/2025), I also generate some boxplots for young and old men.   
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()

#' # Full data set ####
#' ***
load(paste0(data_QC,"/01_Prep_02_UKB_GP_allSamples.RData"))
names(myTab6)
setnames(myTab6,"exposure_value","TC")
setnames(myTab6,"exposure_age","age")

matched = match(myTab6$BSU_ID,myTab7$ID)
myTab6[,sex := myTab7[matched,sex]]
myTab6[sex==0,sex:=2]
myTab6[sex==2,sex2 := "female"]
myTab6[sex==1,sex2 := "male"]

myTab6[age<=45, age3 := "20-45"]
myTab6[age<=50 & age>45, age3 := "46-50"]
myTab6[age<=55 & age>50, age3 := "51-55"]
myTab6[age<=60 & age>55, age3 := "56-60"]
myTab6[age<=65 & age>60, age3 := "61-65"]
myTab6[age<=70 & age>65, age3 := "66-70"]
myTab6[age<=75 & age>70, age3 := "71-75"]
myTab6[age<=80 & age>75, age3 := "76-80"]

#' ## All samples 
ggp5 = ggplot(myTab6, aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (all samples)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp5
ggsave(plot = ggp5, filename = paste0("../results/_figures/Boxplots/AllSamples_all.png"), 
       height = 7, width = 14)

#' ## Men 
ggp6 = ggplot(myTab6[sex==1], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("steelblue"),
                    labels = c("male"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (men)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp6
ggsave(plot = ggp6, filename = paste0("../results/_figures/Boxplots/AllSamples_men.png"), 
       height = 7, width = 14)

#' ## Women 
ggp7 = ggplot(myTab6[sex==2], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred"),
                    labels = c("female"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (women)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp7
ggsave(plot = ggp7, filename = paste0("../results/_figures/Boxplots/AllSamples_women.png"), 
       height = 7, width = 14)

#' ## Post-menopausal women
ggp8 = ggplot(myTab6[flag_post==T & age>=51], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkgreen"),
                    labels = c("female"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (post-menopausal women, age >= 51 years)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp8
ggsave(plot = ggp8, filename = paste0("../results/_figures/Boxplots/AllSamples_women_post.png"), 
       height = 7, width = 14)

#' ## Pre-menopausal women 
ggp9 = ggplot(myTab6[flag_pre==T & age<=60,], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("orange"),
                    labels = c("female"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (pre-menopausal women, age <=60 years)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp9
ggsave(plot = ggp9, filename = paste0("../results/_figures/Boxplots/AllSamples_women_pre.png"), 
       height = 7, width = 14)

#' ## All samples by group
myTab6[sex==1, group := "men"]
myTab6[flag_post==T & age>=51, group := "post-menopausal women"]
myTab6[flag_pre==T & age<=60, group := "pre-menopausal women"]

ggp10 = ggplot(myTab6[!is.na(group)], aes(x=age3, y=TC, fill=group)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("steelblue","darkgreen","orange"),
                    labels = c( "male","female \n(post-menopausal)","female \n(pre-menopausal)"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (all samples)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp10
ggsave(plot = ggp10, filename = paste0("../results/_figures/Boxplots/AllSamples_allByGroup.png"), 
       height = 7, width = 14)

#' # Statin-free subset  ####
#' ***
load(paste0(data_QC,"/01_Prep_02_UKB_GP_noLipidLoweringSubSet.RData"))
names(myTab6)
setnames(myTab6,"exposure_value","TC")
setnames(myTab6,"exposure_age","age")

matched = match(myTab6$BSU_ID,myTab7$ID)
myTab6[,sex := myTab7[matched,sex]]
myTab6[sex==0,sex:=2]
myTab6[sex==2,sex2 := "female"]
myTab6[sex==1,sex2 := "male"]

myTab6[age<=45, age3 := "20-45"]
myTab6[age<=50 & age>45, age3 := "46-50"]
myTab6[age<=55 & age>50, age3 := "51-55"]
myTab6[age<=60 & age>55, age3 := "56-60"]
myTab6[age<=65 & age>60, age3 := "61-65"]
myTab6[age<=70 & age>65, age3 := "66-70"]
myTab6[age<=75 & age>70, age3 := "71-75"]
myTab6[age<=80 & age>75, age3 := "76-80"]

#' ## All samples 
ggp5 = ggplot(myTab6, aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (all samples)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp5
ggsave(plot = ggp5, filename = paste0("../results/_figures/Boxplots/Subset_all.png"), 
       height = 7, width = 14)

#' ## Men 
ggp6 = ggplot(myTab6[sex==1], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("steelblue"),
                    labels = c("male"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (men)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp6
ggsave(plot = ggp6, filename = paste0("../results/_figures/Boxplots/Subset_men.png"), 
       height = 7, width = 14)

#' ## Women 
ggp7 = ggplot(myTab6[sex==2], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkred"),
                    labels = c("female"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (women)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp7
ggsave(plot = ggp7, filename = paste0("../results/_figures/Boxplots/Subset_women.png"), 
       height = 7, width = 14)

#' ## Post-menopausal women
ggp8 = ggplot(myTab6[flag_post==T & age>=51], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("darkgreen"),
                    labels = c("female"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (post-menopausal women, age >= 51 years)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp8
ggsave(plot = ggp8, filename = paste0("../results/_figures/Boxplots/Subset_women_post.png"), 
       height = 7, width = 14)

#' ## Pre-menopausal women 
ggp9 = ggplot(myTab6[flag_pre==T & age<=60,], aes(x=age3, y=TC, fill=sex2)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("orange"),
                    labels = c("female"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (pre-menopausal women, age <=60 years)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp9
ggsave(plot = ggp9, filename = paste0("../results/_figures/Boxplots/Subset_women_pre.png"), 
       height = 7, width = 14)

#' ## All samples by group
myTab6[sex==1, group := "men"]
myTab6[flag_post==T & age>=51, group := "post-menopausal women"]
myTab6[flag_pre==T & age<=60, group := "pre-menopausal women"]

ggp10 = ggplot(myTab6[!is.na(group)], aes(x=age3, y=TC, fill=group)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("steelblue","darkgreen","orange"),
                    labels = c( "male","post-menopausal women","pre-menopausal women"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (all samples)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp10
ggsave(plot = ggp10, filename = paste0("../results/_figures/Boxplots/Subset_allByGroup.png"), 
       height = 7, width = 14)

#' # Poster plots ####
#' ***
load(paste0(data_QC,"/01_Prep_0X_UKB_GP_MiddleAgedMen.RData"))
myTab7_MY = copy(myTab7)
myTab6_MY = copy(myTab6)
load(paste0(data_QC,"/01_Prep_0X_UKB_GP_OldMen.RData"))
myTab7_MO = copy(myTab7)
myTab6_MO = copy(myTab6)

load(paste0(data_QC,"/01_Prep_02_UKB_GP_noLipidLoweringSubSet.RData"))
names(myTab6)
setnames(myTab6,"exposure_value","TC")
setnames(myTab6,"exposure_age","age")

matched = match(myTab6$BSU_ID,myTab7$ID)
myTab6[,sex := myTab7[matched,sex]]
myTab6[sex==0,sex:=2]
myTab6[sex==2,sex2 := "female"]
myTab6[sex==1,sex2 := "male"]

myTab_WY = copy(myTab6)
myTab_WY = myTab_WY[flag_pre == T,]

myTab_WO = copy(myTab6)
myTab_WO = myTab_WO[flag_post == T,]

myTab_MY = myTab6_MY[dumID %in% myTab6$dumID,]
myTab_MO = myTab6_MO[dumID %in% myTab6$dumID,]

myTab_WO[,group:= "post-menopausal"]
myTab_WY[,group:= "pre-menopausal"]
myTab_MO[,group:= "men (>50y)"]
myTab_MY[,group:= "men (<60y)"]

setnames(myTab_MO,"exposure_value","TC")
setnames(myTab_MY,"exposure_value","TC")
setnames(myTab_MO,"exposure_age","age")
setnames(myTab_MY,"exposure_age","age")
myTab = rbind(myTab_MO,myTab_MY,myTab_WO,myTab_WY, fill=T,use.names=T)
myTab[,table(lipLowMed)]

names(myTab)
myTab = myTab[,c(5,2,4,8,13)]

myTab[,min(age)]
myTab[age<=45, age3 := "40-45"]
myTab[age<=50 & age>45, age3 := "46-50"]
myTab[age<=55 & age>50, age3 := "51-55"]
myTab[age<=60 & age>55, age3 := "56-60"]
myTab[age<=65 & age>60, age3 := "61-65"]
myTab[age<=70 & age>65, age3 := "66-70"]

myTab[,age2 := ceiling(age)]
myTab[,table(age2)]

#' ## All samples 
ggp11 = ggplot(myTab[group %in% c("pre-menopausal","men (<60y)")], aes(x=as.factor(age3), y=TC, fill=group)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("#BCF1FF","#FFECC9"),
                    labels = c("men (<60y)", "women, \npre-menopausal"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels in younger people (>51y) without lipid-lowering medication")) +
  theme_classic() +
  ylim(1,11) #+
#theme(legend.position = "none")
ggp11
ggsave(plot = ggp11, filename = paste0("../results/_figures/Boxplots/Subset_young.png"), 
       height = 7, width = 14)

ggp12 = ggplot(myTab[group %in% c("post-menopausal","men (>50y)")], aes(x=as.factor(age3), y=TC, fill=group)) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values = c("#BDD3E5","#B7FFB7"),
                    labels = c("men (>51y)", "women, \npost-menopausal"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels in older people (>51y) without lipid-lowering medication")) +
  theme_classic() +
  ylim(1,11)# +
#theme(legend.position = "none")
ggp12

ggsave(plot = ggp12, filename = paste0("../results/_figures/Boxplots/Subset_old.png"), 
       height = 7, width = 14)


#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

