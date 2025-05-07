#' ---
#' title: "Get exposure data - plots CAD"
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
load(paste0(data_QC,"/01_Prep_01_UKB_filtered_CAD.RData"))
names(myTab)

myTab[sex==0,sex:=2]
myTab[sex==2,sex2 := "female"]
myTab[sex==1,sex2 := "male"]

myTab[,min(age)]
myTab[,max(age)]
myTab[age<=45, age3 := "38-45"]
myTab[age<=50 & age>45, age3 := "46-50"]
myTab[age<=55 & age>50, age3 := "51-55"]
myTab[age<=60 & age>55, age3 := "56-60"]
myTab[age<=65 & age>60, age3 := "61-65"]
myTab[age>65, age3 := "66-72"]

myTab[,table(sex2,age3)]
myTab[,table(group,age3)]

myTab = myTab[!is.na(group) & !is.na(CAD),]

#' # Prevalence per age group and sex ####
#' ***
dumTab1 = myTab[CAD==1,.N,by=c("age3","sex")]
dumTab = myTab[,.N,by=c("age3","sex")]
setorder(dumTab1,sex,age3)
setorder(dumTab,sex,age3)
stopifnot(dumTab$age3 == dumTab1$age3)
stopifnot(dumTab$sex == dumTab1$sex)
dumTab[,CAD := dumTab1$N]
dumTab[,prevalence := CAD/N]
dumTab

dumTab[sex==2,sex2 :="women"]
dumTab[sex==1,sex2 :="men"]
dumTab[,prev2 := round(prevalence,3)*100]

#' # Barplot by age group and sex ####
#' ***
ggp0  = ggplot(data=dumTab, aes(x=age3, y=prev2, fill=sex2)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=prev2), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values = c("steelblue","darkred"),
                    labels = c("male", "female"))+
  labs(x="Age (in years)",
       y="CAD prevalence (in %)", 
       fill="Sex",
       title = paste0("CAD prevalence per age group")) +
  theme_minimal()

ggp0
ggsave(plot = ggp0, filename = paste0("../results/_figures/Barplots/CAD_AgeBySex.png"),
       height = 7, width = 10)

#' # Prevalence per age group and group ####
#' ***
dumTab1 = myTab[CAD==1,.N,by=c("age3","group")]
dumTab = myTab[,.N,by=c("age3","group")]
setorder(dumTab1,group,age3)
setorder(dumTab,group,age3)
stopifnot(dumTab$age3 == dumTab1$age3)
stopifnot(dumTab$group == dumTab1$group)
dumTab[,CAD := dumTab1$N]
dumTab[,prevalence := CAD/N]
dumTab

dumTab[,prev2 := round(prevalence,3)*100]

#' # Barplot by age group and group ####
#' ***
ggp1  = ggplot(data=dumTab, aes(x=age3, y=prev2, fill=group)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=prev2), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values = c("steelblue","darkgreen","orange"),
                    labels = c( "male","female \n(post-menopausal)","female \n(pre-menopausal)"))+
  labs(x="Age (in years)",
       y="CAD prevalence (in %)", 
       fill="Sex",
       title = paste0("CAD prevalence per age group")) +
  theme_minimal()

ggp1
ggsave(plot = ggp1, filename = paste0("../results/_figures/Barplots/CAD_AgeByGroup.png"),
       height = 7, width = 10)

#' # Prevalence per age and sex ####
#' ***
myTab[,age4 := round(age,0)]
dumTab1 = myTab[CAD==1,.N,by=c("age4","sex")]
dumTab = myTab[,.N,by=c("age4","sex")]
setorder(dumTab1,sex,age4)
setorder(dumTab,sex,age4)
dumTab1[,dumID := paste(sex,age4,sep=" -")]
dumTab[,dumID := paste(sex,age4,sep=" -")]
matched = match(dumTab$dumID,dumTab1$dumID)
table(dumTab$age4 == dumTab1$age4[matched])
table(dumTab$sex == dumTab1$sex[matched])

dumTab[,CAD := dumTab1[matched,N]]
dumTab[,prevalence := CAD/N]

dumTab[sex==2,sex2 :="women"]
dumTab[sex==1,sex2 :="men"]
dumTab[,prev2 := round(prevalence,3)*100]
dumTab = dumTab[!is.na(prev2),]
dumTab = dumTab[prev2 !=100,]

#' # Barplot by age group and sex ####
#' ***
ggp2  = ggplot(data=dumTab, aes(x=age4, y=prev2, fill=sex2)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=prev2), vjust=1.6, color="white",
            position = position_dodge(0.9), size=2)+
  scale_fill_manual(values = c("steelblue","darkred"),
                    labels = c("male", "female"))+
  labs(x="Age (in years)",
       y="CAD prevalence (in %)", 
       fill="Sex",
       title = paste0("CAD prevalence per age")) +
  theme_minimal()

ggp2
ggsave(plot = ggp2, filename = paste0("../results/_figures/Barplots/CAD_AgeBySex2.png"),
       height = 7, width = 15)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

