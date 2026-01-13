#' ---
#' title: "Get UKB data"
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
#' I need to rerun this in an attempt to harmonize my group criteria for young and old men (should be similar to what I said on the poster). 
#' 
#' Step 1: Load UKB data and reduce to 
#'    - White British ancestry
#'    - No pairwise kinship
#'    - LDL-C data available at baseline
#'    - consent still given
#'    
#' Step 2: define groups (men overlap will be randomly split)
#'    - post-menopausal women at baseline: "had menopause?" == "yes" (coded as 1) & age >50
#'    - pre-menopausal women at baseline: "had menopause?" == "no (coded as 0) & age<61
#'    - elderly men at baseline: age >50 
#'    - younger men at baseline: age <61
#'    
#' Step 3: Load GP data and reduce to ID overlap
#'    
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # UKB data ####
#' ***
#' ## Load data 
myTab_20 = fread(UKB_phenotypes, header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:18506]

exposure = c("f.30690.0.0","f.30690.1.0")
table(is.element(exposure,myAnnot$colNm))

covars = c("f.eid", "f.31.0.0","f.53.0.0","f.53.1.0","f.2724.0.0","f.2724.1.0",
           paste0("f.20003.0.",c(0:47)),paste0("f.20003.1.",c(0:47)),
           "f.21000.0.0","f.21022.0.0","f.22000.0.0","f.22001.0.0",
           paste0("f.22009.0.",1:10), "f.22021.0.0")
table(is.element(covars,myAnnot$colNm))

myAnnot = myAnnot[colNm %in% c(exposure,covars),]

x = myAnnot[,colNR]
myTab_cross <- fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
names(myTab_cross)
names(myTab_cross) = c("ID","sex","date_init","date_FU","menopause_init","menopause_FU",
                       paste("medication_init",1:48,sep="_"),
                       paste("medication_FU",1:48,sep="_"),
                       "ancestry","age","GenotypeBatch","GeneticSex",
                       paste("PC",1:10,sep="_"),
                       "kinship","TC_init","TC_FU")

save(myTab_cross, file = paste0(data_QC,"UKB/01_Prep_TrajGWAS_unfiltered.RData"))
#load(paste0(data_QC,"UKB/01_Prep_TrajGWAS_unfiltered.RData"))

#' ## Filter data
#' Ancestry, kinship, baseline TC, consent, sex in database is equal to genetic sex
myTab_UKB = copy(myTab_cross) 
myTab_UKB = myTab_UKB[ancestry == 1001,]
myTab_UKB = myTab_UKB[kinship == 0,]
myTab_UKB = myTab_UKB[!is.na(TC_init),]
myTab_UKB[,table(is.na(TC_FU))]
ToExclude = fread(gsub("ukb672224.tab","withdraw98032_20250818.csv",UKB_phenotypes))
myTab_UKB = myTab_UKB[!is.element(ID,ToExclude$V1),]
myTab_UKB = myTab_UKB[sex == GeneticSex,]

#' ## Get medication 
#' C10 = lipid modifying agents
#' 
codingTable = data.table(read_xlsx(UKB_MedCoding,sheet=1))
lipidLowering = codingTable[grepl("C10",ATC_Code),Coding]

myMeds = names(myTab_UKB)[grep("medication_init",names(myTab_UKB))]
myTab_UKB[,lipidLow_init := 0]
for(i in 1:length(myMeds)){
  #i=1
  myTab_UKB[get(myMeds[i]) %in% lipidLowering,lipidLow_init := 1]
}
myTab_UKB[,get("myMeds"):=NULL]

myMeds = names(myTab_UKB)[grep("medication_FU",names(myTab_UKB))]
myTab_UKB[,lipidLow_FU := 0]
for(i in 1:length(myMeds)){
  #i=1
  myTab_UKB[get(myMeds[i]) %in% lipidLowering,lipidLow_FU := 1]
}
myTab_UKB[,get("myMeds"):=NULL]
dim(myTab_UKB)
myTab_UKB[,table(lipidLow_init,lipidLow_FU)]

#' ## Check genotyping batch
myTab_UKB[,table(GenotypeBatch)]
myTab_UKB[,Array := "Axiom"]
myTab_UKB[GenotypeBatch<0,Array := "BiLEVE"]

#' # Get groups ####
#' ***
#' ## Women
#' - post-menopausal women at baseline: "had menopause?" == "yes" (coded as 1) & age >50
#' - pre-menopausal women at baseline: "had menopause?" == "no (coded as 0) & age<61
#' 
myTab_UKB[,table(sex,menopause_init)]
myTab_UKB[,table(sex,menopause_FU)]

myTab_UKB[sex==0 & menopause_init==0 & age <=60, group := "premenopausal"]
myTab_UKB[sex==0 & menopause_init==1 & age >=51, group := "postmenopausal"]
myTab_UKB[group=="postmenopausal" & menopause_FU!=1, group := NA]

myTab_UKB[sex==0,table(is.na(group))]
myTab_UKB[sex==0,table((group))]

#' ## Men
#' - elderly men at baseline: age >50 
#' - younger men at baseline: age <61
#' 
myTab_UKB[sex==1 & age <50, group := "younger"]
myTab_UKB[sex==1 & age >60, group := "elderly"]

myTab_UKB[sex==1,table(is.na(group))]
myTab_UKB[sex==1,table((group))]

#' Okay, to split the overlapping IDs so that we have a similar young/old ratio as in women (28% young, 72% old)
myIDs = myTab_UKB[is.na(group) & sex==1,ID]

set.seed(2025)
mySample = sample(1:length(myIDs),size = 8000,replace = F)

myTab_UKB[ID %in% myIDs[mySample], group := "younger"]
myTab_UKB[sex==1 & is.na(group), group := "elderly"]

myTab_UKB[sex==1,table(is.na(group))]
myTab_UKB[sex==1,table((group))]

#' ## Reduce to relevant columns
names(myTab_UKB)
myTab_UKB[,ancestry := NULL]
myTab_UKB[,kinship := NULL]
myTab_UKB[,menopause_init := NULL]
myTab_UKB[,menopause_FU := NULL]
myTab_UKB[,GenotypeBatch := NULL]
myTab_UKB[,GeneticSex := NULL]

#' # GP data  ####
#' ***
#' total cholesterol: cholesterol codes from Denaxas supplement table
#' Lipid-Regulating Drugs: BNF code starting with "02.12"
#' 
myTab_GP_TC = fread(UKB_GPdata)
myTab_GP_TC = myTab_GP_TC[eid %in% myTab_UKB$ID,]
myTab_GP_TC[,length(unique(eid))]

#' Load coding table for primary care data (Denaxas et al, Supplementary Table 2)
codingTable = fread(UKB_GPcoding)
codingTable[,table(duplicated(readcode))]

#' Filter for reads from the TC coding
myTab_GP_TC = myTab_GP_TC[read_3 %in% codingTable$readcode | read_2 %in% codingTable$readcode,]
myTab_GP_TC[,table(read_3)]
myTab_GP_TC[,table(read_2)]

#' Check that there is information on ID, date and TC value
names(myTab_GP_TC)
myTab_GP_TC[,table(is.na(eid))]

myTab_GP_TC[,table(is.na(event_dt))]
myTab_GP_TC[,table(event_dt=="")]
myTab_GP_TC[event_dt==""]
myTab_GP_TC = myTab_GP_TC[event_dt!=""]

myTab_GP_TC[,table(is.na(value1))]
myTab_GP_TC[,table(value1=="")]
myTab_GP_TC[,table(is.na(value2))]
myTab_GP_TC[,table(value2=="")]
myTab_GP_TC[,table(is.na(value3))]
myTab_GP_TC[,table(value3=="")]

myTab_GP_TC[value2 != "",table(value1)]
myTab_GP_TC[value3 != "",table(value1)]
myTab_GP_TC[,table(value3=="",value2=="")]
myTab_GP_TC[(value1 == "" | value1=="OPR003") & value2 != "",value1 := value2]
myTab_GP_TC[(value1 == "" | value1=="OPR003") & value3 != "",value1 := value3]

myTab_GP_TC = myTab_GP_TC[value1 != ""]
myTab_GP_TC[,value1 := as.numeric(value1)]
myTab_GP_TC = myTab_GP_TC[!is.na(value1),]
myTab_GP_TC[,table(value1==0)]
myTab_GP_TC = myTab_GP_TC[value1!=0]
myTab_GP_TC[,table(value1<0)]

#' Check distribution of TC values - similar to Ko et al., I filter everything above 30
hist(myTab_GP_TC$value1)
hist(myTab_GP_TC[value1<30,value1])
myTab_GP_TC = myTab_GP_TC[value1<30,]

#' Now reduce to relevant columns: ID, date, TC value
myTab_GP_TC = myTab_GP_TC[,c(1,3,6)]
names(myTab_GP_TC) = c("ID","date","TC")

#' Now I can recreate the time since UKB baseline
myTab_GP_TC[,day := substr(date,1,2)]
myTab_GP_TC[,month := substr(date,4,5)]
myTab_GP_TC[,year := substr(date,7,10)]
myTab_GP_TC[,date2 := paste(year,month,day,sep = "-")]
myTab_GP_TC[,date2 := as.Date(date2)]
matched = match(myTab_GP_TC$ID,myTab_UKB$ID)
table(is.na(matched))
myTab_GP_TC[,date_UKBBaseline := myTab_UKB[matched,date_init]]
myTab_GP_TC[,age_UKBBaseline := myTab_UKB[matched,age]]
myTab_GP_TC[,date_UKBBaseline := as.Date(date_UKBBaseline)]
myTab_GP_TC[,dateDif := date2 - date_UKBBaseline]
myTab_GP_TC[,dateDifYears := as.numeric(dateDif)/365]
myTab_GP_TC[,age := age_UKBBaseline + dateDifYears]
myTab_GP_TC[dateDifYears>0,flag := "after baseline"]
myTab_GP_TC[dateDifYears<0,flag := "before baseline"]
myTab_GP_TC[dateDifYears==0,]

#' There are some entries of the exact same date as of the UKB baseline. In these cases I keep the UKB data
myTab_GP_TC = myTab_GP_TC[dateDifYears!=0,]

names(myTab_GP_TC)
myTab_GP_TC = myTab_GP_TC[,c(1,7,3,12,13)]
myTab_GP_TC[,table(is.na(flag))]
myTab_GP_TC[,table(is.na(age))]
myTab_GP_TC[,table(age<0)]
myTab_GP_TC[age<0,table(date2)]
myTab_GP_TC = myTab_GP_TC[age>0,]
myTab_GP_TC[,min(age)]
myTab_GP_TC[,max(age)]
myTab_GP_TC = myTab_GP_TC[age>20,]

#' Check duplicates
myTab_GP_TC[,dumID := paste(ID,date2,sep="__")]
myTab_GP_TC[,table(duplicated(dumID))]
dupIDs = myTab_GP_TC[duplicated(dumID),dumID]
myTab_GP_TC[dumID %in% dupIDs,]
myTab_GP_TC = myTab_GP_TC[!duplicated(dumID),]
myTab_GP_TC[,dumID := NULL]

#' Get lipid-lowering medication
myTab_GP_scripts = fread(UKB_GPscripts)
myTab_GP_scripts = myTab_GP_scripts[eid %in% myTab_UKB$ID,]
myTab_GP_scripts[,length(unique(eid))]
myTab_GP_scripts = myTab_GP_scripts[substr(bnf_code,1,5) == "02.12",]

#' Now I have to add the lipid lowering treatment information
myTab_GP_scripts[,table(is.na(issue_date))]
myTab_GP_scripts[,table(issue_date=="")]
myTab_GP_scripts[issue_date==""]
myTab_GP_scripts = myTab_GP_scripts[issue_date!=""]

myTab_GP_scripts[,day := substr(issue_date,1,2)]
myTab_GP_scripts[,month := substr(issue_date,4,5)]
myTab_GP_scripts[,year := substr(issue_date,7,10)]
myTab_GP_scripts[,date2 := paste(year,month,day,sep = "-")]
myTab_GP_scripts[,date2 := as.Date(date2)]
myTab_GP_scripts[,table(is.na(date2))]
myTab_GP_scripts[,table(date2=="")]
myTab_GP_scripts[,table(grepl("1901",date2))]
myTab_GP_scripts[,table(grepl("1902",date2))]
myTab_GP_scripts[,table(grepl("1903",date2))]
myTab_GP_scripts = myTab_GP_scripts[!grepl("1902",date2)]
myTab_GP_scripts = myTab_GP_scripts[!grepl("2037",date2)]

#' I will assume that after the first prescription they are always on treatment
FirstScript = myTab_GP_scripts[,min(date2),by = eid]
FirstScript = FirstScript[eid %in% myTab_GP_TC$ID,]
myTab_GP_TC[!is.element(ID,FirstScript$eid),statin := 0]

for(i in 1:dim(FirstScript)[1]){
  #i=1
  myTab_GP_TC[ID == FirstScript[i,eid] & date2<=FirstScript[i,V1],statin := 0]
  myTab_GP_TC[ID == FirstScript[i,eid] & date2>FirstScript[i,V1],statin := 1]
  
}
myTab_GP_TC = myTab_GP_TC[,c(1,2,3,4,6,5)]
names(myTab_GP_TC)
setnames(myTab_GP_TC,"date2","date")
myTab_GP_TC[,table(is.na(statin))]
myTab_GP_TC[,table(statin)]

#' # Merge UKB and GP data ####
#' ***
myTab_Baseline = copy(myTab_UKB)
myTab_Baseline = myTab_Baseline[ID %in% myTab_GP_TC$ID,]
myTab_Baseline = myTab_Baseline[,c(1,3,16,5,18)]
myTab_Baseline[,flag:="UKB baseline"]
names(myTab_Baseline) = names(myTab_GP_TC)
myTab_Baseline[,date := as.Date(date)]

myTab_FollowUp = copy(myTab_UKB)
myTab_FollowUp = myTab_FollowUp[ID %in% myTab_GP_TC$ID & !is.na(TC_FU),]
myTab_FollowUp[,dateDif := date_FU - date_init]
myTab_FollowUp[,dateDifYears := as.numeric(dateDif)/365]
myTab_FollowUp[,age_FU := age + dateDifYears]
myTab_FollowUp = myTab_FollowUp[,c(1,4,17,24,19)]
myTab_FollowUp[,flag:="UKB follow-up"]
names(myTab_FollowUp) = names(myTab_GP_TC)
myTab_FollowUp[,date := as.Date(date)]

myTab_long = rbind(myTab_Baseline,myTab_FollowUp,myTab_GP_TC)
setorder(myTab_long,ID,date)

#' Match fixed covariables (sex, group, array, genetic PCs)
matched = match(myTab_long$ID,myTab_UKB$ID)
table(is.na(matched))
stopifnot(myTab_long$ID == myTab_UKB[matched,ID])
myTab_long[,sex := myTab_UKB[matched,sex]]
myTab_long[sex==0,sex := 2]
myTab_long[,group := myTab_UKB[matched,group]]
myTab_long[group=="postmenopausal",group := "older"]
myTab_long[group=="elderly",group := "older"]
myTab_long[group=="premenopausal",group := "younger"]
myTab_long[,table(sex,group)]
myTab_long[,table(sex,is.na(group))]
myTab_long[,array := myTab_UKB[matched,Array]]
myTab_long = cbind(myTab_long,myTab_UKB[matched,6:15])

#' # Checks ####
#' *** 
#' Check 1: all IDs must have UKB baseline
test = myTab_long[,.N,by=c("ID","flag")]
test[,table(flag)]
myTab_long[,length(unique(ID))]

#' Check 2: check range of TC values
hist(myTab_long$TC)
summary(myTab_long$TC)

#' I will restrict everything to mean + 8*SD from the UKB baseline data (whole UKB, see https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30690)
upperLimit = 5.69263 + 8*1.14657
myTab_long[,table(TC>upperLimit)]
myTab_long = myTab_long[TC<=upperLimit,]
hist(myTab_long$TC)
summary(myTab_long$TC)

#' Check 3: statin always 1 after first mentioned
test1 = myTab_long[statin==0,max(date),by=ID]
test2 = myTab_long[statin==1,min(date),by=ID]
matched = match(test2$ID,test1$ID)
test2[,date_noStatin := test1[matched,V1]]
test2 = test2[!is.na(date_noStatin),]
test2 = test2[date_noStatin>V1,]

for(i in 1:dim(test2)[1]){
  #i=1
  myTab_long[ID == test2[i,ID] & date<test2[i,V1],statin := 0]
  myTab_long[ID == test2[i,ID] & date>=test2[i,V1],statin := 1]
  
}

test1 = myTab_long[statin==0,max(date),by=ID]
test2 = myTab_long[statin==1,min(date),by=ID]
matched = match(test2$ID,test1$ID)
test2[,date_noStatin := test1[matched,V1]]
test2 = test2[!is.na(date_noStatin),]
test2 = test2[date_noStatin>V1,]
test2

#' Check 4: number of time points per individual and setting
TP_overall = myTab_long[,.N,by="ID"]
hist(TP_overall$N)
summary(TP_overall$N)

#' Okay, because of the extreme value filtering there are some with just one value. I will remove them
TP_overall[,table(N==1)]
myTab_long = myTab_long[!is.element(ID, TP_overall[N==1,ID])]
TP_overall = myTab_long[,.N,by=c("ID","sex")]
TP_overall[,table(N==1,sex)]
TP_overall[,table(N>=3,sex)]

#' Now check for young ones
TP_younger = myTab_long[group=="younger" & flag %in% c("before baseline","UKB baseline"),.N,by=c("ID","sex")]
hist(TP_younger$N)
summary(TP_younger$N)
TP_younger[,table(N==1,sex)]
TP_younger[,table(N>=3,sex)]
TP_younger[,table(N>=2,sex)]

#' Now check for old ones
TP_older = myTab_long[group=="older" & flag %in% c("after baseline","UKB baseline","UKB follow-up"),.N,by=c("ID","sex")]
hist(TP_older$N)
summary(TP_older$N)
TP_older[,table(N==1,sex)]
TP_older[,table(N>=3,sex)]
TP_older[,table(N>=2,sex)]

#' # Save data ####
#' ***
myTab_long[,age2 := age^2]
myTab_long[,plot(age,age2)]

save(myTab_long, file = paste0(data_QC,"/UKB/01_Prep_TrajGWAS.RData"))
write.table(myTab_long[,unique(ID)],file = paste0(data_QC,"/UKB/01_Prep_TrajGWAS_SampleList.txt"), 
            col.names = F, row.names = F, quote = F)

#' # Create input for trajGWAS ####
#' ***
#' I need to create the csv for the trajGWAS input and the sample list for the PLINK extraction. In addition, I want some basic description per sub-population as well (sample size, mean TC, sd TC, number of time points, average age when statin was first described, average age at first time point, average age at last time point)
#' 

groups = c("men","women","older men","postmenopausal","premenopausal","young men")
groups2 = c("men","women","old","post","pre","young")

dumTab = foreach(i = 1:length(groups))%do%{
  #i=1
  data = copy(myTab_long)
  if(i %in% c(1,3,6)) data = data[sex==1,]
  if(i %in% c(2,4,5)) data = data[sex==2,]
  if(i %in% c(3,4)) data = data[group=="older" & flag %in% c("after baseline","UKB baseline","UKB follow-up"),]
  if(i %in% c(5,6)) data = data[group=="younger" & flag %in% c("before baseline","UKB baseline"),]
  test = data[,.N,by=c("ID")]
  data = data[ID %in% test[N>1,ID]]
  
  # Save group specific data
  write.csv(data,
            file = paste0(data_QC,"/UKB/01_Prep_TrajGWAS_",groups2[i],".csv"), 
            row.names = FALSE)
  write.table(data[,unique(ID)],
              file = paste0(data_QC,"/UKB/01_Prep_TrajGWAS_SampleList_",groups2[i],".txt"), 
              col.names = F, row.names = F, quote = F)
  dum1 = data[,min(age),by=ID]
  dum2 = data[,max(age),by=ID]
  dum3 = data[,.N,by=ID]
  setnames(dum1,"V1","minAge")
  dum1[,maxAge := dum2$V1]
  dum1[,timeDiff := maxAge - minAge]
  dum1[,timePoints := dum3$N]
  dum4 = data[statin==1,min(age),by = ID]
  
  res1 = data.table(setting = groups2[i],
                    sampleSize = dim(dum1)[1],
                    meanTC = data[,mean(TC)],
                    sdTC = data[,sd(TC)],
                    observations = dim(data)[1],
                    timepoints_median = dum1[,as.numeric(summary(timePoints)[3])],
                    timepoints_1st = dum1[,as.numeric(summary(timePoints)[2])],
                    timepoints_3rd = dum1[,as.numeric(summary(timePoints)[5])],
                    age_1st_TP = dum1[,mean(minAge)],
                    age_last_TP = dum1[,mean(maxAge)],
                    statin_treated = dim(dum4)[1],
                    age_1st_statin = dum4[,mean(V1)])
  
  res1
  
}

res = rbindlist(dumTab)
res
save(res,file="../results/UKB/01_descriptiveTable.RData")

#' Simple plots
#' 
myTab_long[,group2 := paste(sex,group,sep=" - ")]
myTab_long[,group2 := gsub("1 - ","males - ",group2)]
myTab_long[,group2 := gsub("2 - ","females - ",group2)]

ggp0  = ggplot(myTab_long[!is.na(group),], aes(x=age, y=TC, group=ID, col=group2, fill=group2, shape=as.factor(statin))) +
  facet_wrap(~group2,scales="free") +
  geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point()+
  labs(x="Age (years)", y="TC (mmol/L)", title = "All samples") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(filename = "../results/UKB/01_Trajectories_age_all.png",plot = ggp0)

ggp1  = ggplot(myTab_long[!is.na(group) & statin==0,], aes(x=age, y=TC, group=ID, col=group2, fill=group2)) +
  facet_wrap(~group2,scales="free") +
  geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point()+
  labs(x="Age (years)", y="TC (mmol/L)", title = "Statin-free samples") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(filename = "../results/UKB/01_Trajectories_age_noStatin.png",plot = ggp1)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
