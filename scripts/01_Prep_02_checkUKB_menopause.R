#' ---
#' title: "Get exposure data - strata specific"
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
#' In the previous script, I extract the potential samples. This should be the same sample size as in my previous project. Now I want to generate the relevant subsets: 
#' 
#' - men 
#' - women
#' - premenopausal women (EHR - UKB baseline, including women answering the question "had menopause" with "no" at baseline & age <= 60 years)
#' - postmenopausal women (UKB baseline - EHR, including women answering the question "had menopause" with "yes" at baseline & age >= 51 years )
#' 
#' For all 4 subsets, I want two settings: 
#' 
#' - all samples
#' - age within 40-70, no lipid-lowering medication, first TC value per year 
#'  
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()

#' # Load UKB data ####
#' ***
load(paste0(data_QC,"/01_Prep_01_UKB_GP_TC_GLGC.RData"))
myTab7[,table(sex)]
myTab7[,table(sex,group)]

#' # Stratify women ####
#' ***
#' ## Postmenopausal women ####
myTab7_post = copy(myTab7)
myTab7_post = myTab7_post[group=="postmenopausal" & age>=51,]

myTab6_post = copy(myTab6)
myTab6_post = myTab6_post[BSU_ID %in% myTab7_post$ID,]

counter = seq(1,dim(myTab7_post)[1],1000)
for(i in 1:dim(myTab7_post)[1]){
  if(i %in% counter)message("Working on i=",i)
  myID = myTab7_post[i,ID]
  myAge = myTab7_post[i,age]
  myTab6_post = myTab6_post[!(BSU_ID==myID & exposure_age<myAge),]
}

#' Now I have to check if all individuals still have 3 or more observations
#' 
observations_post = myTab6_post[,.N,by=BSU_ID]
hist(observations_post$N)
table(observations_post$N>2)
myTab6_post = myTab6_post[BSU_ID %in% observations_post[N>2,BSU_ID],]
myTab7_post = myTab7_post[ID %in% observations_post[N>2,BSU_ID],]
observations_post = myTab6_post[,.N,BSU_ID]
summary(observations_post$N)

#' ## Premenopausal women ####
myTab7_pre = copy(myTab7)
myTab7_pre = myTab7_pre[group=="premenopausal" & age<=60,]

myTab6_pre = copy(myTab6)
myTab6_pre = myTab6_pre[BSU_ID %in% myTab7_pre$ID,]

counter = seq(1,dim(myTab7_pre)[1],500)
for(i in 1:dim(myTab7_pre)[1]){
  if(i %in% counter)message("Working on i=",i)
  myID = myTab7_pre[i,ID]
  myAge = myTab7_pre[i,age]
  myTab6_pre = myTab6_pre[!(BSU_ID==myID & exposure_age>myAge),]
}
#' Now I have to check if all individuals still have 2 or more observations (more relaxed, because it already is a small data set)
#' 
observations_pre = myTab6_pre[,.N,by=BSU_ID]
hist(observations_pre$N)
table(observations_pre$N>1)
myTab6_pre = myTab6_pre[BSU_ID %in% observations_pre[N>1,BSU_ID],]
myTab7_pre = myTab7_pre[ID %in% observations_pre[N>1,BSU_ID],]
observations_pre = myTab6_pre[,.N,BSU_ID]
summary(observations_pre$N)

#' ## Add flags to main data sets ####
myTab6[,table(duplicated(dumID))]

myTab6[,flag_post := F]
myTab6[dumID %in% myTab6_post$dumID,flag_post := T]

myTab6[,flag_pre := F]
myTab6[dumID %in% myTab6_pre$dumID,flag_pre := T]

myTab7[,flag_post := F]
myTab7[ID %in% myTab7_post$ID,flag_post := T]

myTab7[,flag_pre := F]
myTab7[ID %in% myTab7_pre$ID,flag_pre := T]

#' # Save all samples approach ####
#' ***
#' Sample size per strata
#' 
#' - men: 35,726
#' - women: 38,052
#' - post-menopausal women: 21,356
#' - pre-menopausal women: 3,642
#' 
save(myTab6, myTab7, file = paste0(data_QC,"/01_Prep_02_UKB_GP_allSamples.RData"))

#' # Get subset ####
#' ***
#' age within 40-70, no lipid-lowering medication ever, 1 observation per year 
#' 
#' ## No lipid-lowering medication
myTab6 = myTab6[lipLowMed==0,]

#' ## 1 observation per year
myTab6[,age2:= ceiling(exposure_age)]
myTab6[,dumID2 := paste0(BSU_ID,"_",age2)]
myTab6 = myTab6[!duplicated(dumID2),]

#' ## Age group
myTab6[,min(age2)]
myTab6[,max(age2)]
myTab6 = myTab6[age2>=40 & age2<=70,]

#' ## Check number of observations
observations = myTab6[,.N,by=BSU_ID]
hist(observations$N)
table(observations$N>2)
myTab6 = myTab6[BSU_ID %in% observations[N>2,BSU_ID],]
myTab7 = myTab7[ID %in% observations[N>2,BSU_ID],]
observations = myTab6[,.N,BSU_ID]
summary(observations$N)

#' ## Stratify again ####
#' 
#' ### Postmenopausal women ####
myTab7_post = copy(myTab7)
myTab7_post = myTab7_post[group=="postmenopausal" & age>=51,]

myTab6_post = copy(myTab6)
myTab6_post = myTab6_post[BSU_ID %in% myTab7_post$ID,]

counter = seq(1,dim(myTab7_post)[1],1000)
for(i in 1:dim(myTab7_post)[1]){
  if(i %in% counter)message("Working on i=",i)
  myID = myTab7_post[i,ID]
  myAge = myTab7_post[i,age]
  myTab6_post = myTab6_post[!(BSU_ID==myID & exposure_age<myAge),]
}

#' Now I have to check if all individuals still have 3 or more observations
#' 
observations_post = myTab6_post[,.N,by=BSU_ID]
hist(observations_post$N)
table(observations_post$N>2)
myTab6_post = myTab6_post[BSU_ID %in% observations_post[N>2,BSU_ID],]
myTab7_post = myTab7_post[ID %in% observations_post[N>2,BSU_ID],]
observations_post = myTab6_post[,.N,BSU_ID]
summary(observations_post$N)

#' ### Premenopausal women ####
myTab7_pre = copy(myTab7)
myTab7_pre = myTab7_pre[group=="premenopausal" & age<=60,]

myTab6_pre = copy(myTab6)
myTab6_pre = myTab6_pre[BSU_ID %in% myTab7_pre$ID,]

counter = seq(1,dim(myTab7_pre)[1],500)
for(i in 1:dim(myTab7_pre)[1]){
  if(i %in% counter)message("Working on i=",i)
  myID = myTab7_pre[i,ID]
  myAge = myTab7_pre[i,age]
  myTab6_pre = myTab6_pre[!(BSU_ID==myID & exposure_age>myAge),]
}
#' Now I have to check if all individuals still have 2 or more observations (more relaxed, because it already is a small data set)
#' 
observations_pre = myTab6_pre[,.N,by=BSU_ID]
hist(observations_pre$N)
table(observations_pre$N>1)
myTab6_pre = myTab6_pre[BSU_ID %in% observations_pre[N>1,BSU_ID],]
myTab7_pre = myTab7_pre[ID %in% observations_pre[N>1,BSU_ID],]
observations_pre = myTab6_pre[,.N,BSU_ID]
summary(observations_pre$N)

#' # Add flags to main data sets ####
#' ***
myTab6[,table(duplicated(dumID))]

myTab6[,flag_post := F]
myTab6[dumID %in% myTab6_post$dumID,flag_post := T]

myTab6[,flag_pre := F]
myTab6[dumID %in% myTab6_pre$dumID,flag_pre := T]

myTab7[,flag_post := F]
myTab7[ID %in% myTab7_post$ID,flag_post := T]

myTab7[,flag_pre := F]
myTab7[ID %in% myTab7_pre$ID,flag_pre := T]

#' # Save sub-set approach ####
#' ***
#' Sample size per strata
#' 
#' - men:   35,726    --> 23,524      
#' - women: 38,052    --> 27,666
#' - post-menopausal women: 21,987    --> 10,154
#' - pre-menopausal women:   3,720    -->  2,365
#' 
save(myTab6, myTab7, file = paste0(data_QC,"/01_Prep_02_UKB_GP_noLipidLoweringSubSet.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

