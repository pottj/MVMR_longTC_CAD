#' ---
#' title: "Figures - trajectories"
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
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()

#' # Full data set  ####
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
ggsave(plot = ggp0, filename = paste0("../results/_figures/Trajectories/Subset_all.png"), 
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
ggsave(plot = ggp1, filename = paste0("../results/_figures/Trajectories/Subset_men.png"), 
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
ggsave(plot = ggp2, filename = paste0("../results/_figures/Trajectories/Subset_women.png"), 
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
ggsave(plot = ggp3, filename = paste0("../results/_figures/Trajectories/Subset_women_post.png"), 
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
ggsave(plot = ggp4, filename = paste0("../results/_figures/Trajectories/Subset_women_pre.png"), 
       height = 7, width = 14)


#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

