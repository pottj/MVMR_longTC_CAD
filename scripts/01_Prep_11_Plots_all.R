#' ---
#' title: "Get exposure data - plots"
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
.libPaths()

#' # Load data  ####
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

#' # Trajectories ####
#' ***
#' ## All samples 
ggp0  = ggplot(myTab6, aes(x=age, y=TC, group=BSU_ID, col=sex2, 
                                    fill=sex2, shape = as.factor(lipLowMed))) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(color="black")+
  labs(x="Age (years)", y="TC (mmol/L)", title = "All samples") +
  scale_shape_manual(values = c(21,24),
                     labels = c("no treatment","statin treatment"))+
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("female", "male"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("female", "male"))+
  theme_classic() +
  theme(legend.position = "none")
ggp0
ggsave(plot = ggp0, filename = paste0("../results/_figures/Trajectories/AllSamples_all.png"), 
       height = 7, width = 14)

#' ## Men 
ggp1  = ggplot(myTab6[sex==1,], aes(x=age, y=TC, group=BSU_ID, col=sex2, 
                                    fill=sex2, shape = as.factor(lipLowMed))) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(color="black")+
  labs(x="Age (years)", y="TC (mmol/L)", "Men") +
  scale_shape_manual(values = c(21,24),
                     labels = c("no treatment","statin treatment"))+
  scale_fill_manual(values = c("steelblue"),
                    labels = c("male"))+
  scale_colour_manual(values = c("steelblue"),
                      labels = c("male"))+
  theme_classic() +
  theme(legend.position = "none")
ggp1
ggsave(plot = ggp1, filename = paste0("../results/_figures/Trajectories/AllSamples_men.png"), 
       height = 7, width = 14)

#' ## Women 
ggp2  = ggplot(myTab6[sex==2,], aes(x=age, y=TC, group=BSU_ID, col=sex2, 
                                    fill=sex2, shape = as.factor(lipLowMed))) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(color="black")+
  labs(x="Age (years)", y="TC (mmol/L)",title = "Women") +
  scale_shape_manual(values = c(21,24),
                     labels = c("no treatment","statin treatment"))+
  scale_fill_manual(values = c("darkred"),
                    labels = c("female"))+
  scale_colour_manual(values = c("darkred"),
                      labels = c("female"))+
  theme_classic() +
  theme(legend.position = "none")
ggp2
ggsave(plot = ggp2, filename = paste0("../results/_figures/Trajectories/AllSamples_women.png"), 
       height = 7, width = 14)

#' ## Post-menopausal women
ggp3  = ggplot(myTab6[flag_post==T & age>=51,], aes(x=age, y=TC, group=BSU_ID, col=sex2, 
                                    fill=sex2, shape = as.factor(lipLowMed))) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(color="black")+
  labs(x="Age (years)", y="TC (mmol/L)",title = "Post-menopausal women") +
  scale_shape_manual(values = c(21,24),
                     labels = c("no treatment","statin treatment"))+
  scale_fill_manual(values = c("darkgreen"),
                    labels = c("female"))+
  scale_colour_manual(values = c("darkgreen"),
                      labels = c("female"))+
  theme_classic() +
  theme(legend.position = "none")
ggp3
ggsave(plot = ggp3, filename = paste0("../results/_figures/Trajectories/AllSamples_women_post.png"), 
       height = 7, width = 14)

#' ## Pre-menopausal women
ggp4  = ggplot(myTab6[flag_pre==T & age<=60,], aes(x=age, y=TC, group=BSU_ID, col=sex2, 
                                          fill=sex2, shape = as.factor(lipLowMed))) +
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point(color="black")+
  labs(x="Age (years)", y="TC (mmol/L)",title = "Pre-menopausal women") +
  scale_shape_manual(values = c(21,24),
                     labels = c("no treatment","statin treatment"))+
  scale_fill_manual(values = c("orange"),
                    labels = c("female"))+
  scale_colour_manual(values = c("orange"),
                      labels = c("female"))+
  theme_classic() +
  theme(legend.position = "none")
ggp4
ggsave(plot = ggp4, filename = paste0("../results/_figures/Trajectories/AllSamples_women_pre.png"), 
       height = 7, width = 14)

#' # Boxplots ####
#' ***
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
                    labels = c( "male","post-menopausal women","pre-menopausal women"))+
  labs(x="Age (grouped)",
       y="Total cholesterol (in mmol/l)", 
       fill="Sex",
       title = paste0("Boxplots of TC levels per age group (all samples)")) +
  theme_classic() +
  theme(legend.position = "none")
ggp10
ggsave(plot = ggp10, filename = paste0("../results/_figures/Boxplots/AllSamples_allByGroup.png"), 
       height = 7, width = 14)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

