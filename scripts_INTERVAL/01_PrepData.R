#' ---
#' title: "Get exposure data in INTERVAL study"
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
#' Step 1: Load INTERVAL data and reduce to 
#'    - White British ancestry
#'    - LDL-C data available at baseline
#'    
#' Step 2: define groups (men overlap will be randomly split)
#'    - post-menopausal women at baseline: "had menopause?" == "yes" (coded as 1) & age >50
#'    - pre-menopausal women at baseline: "had menopause?" == "no (coded as 0) & age<61
#'    - elderly men at baseline: age >50 
#'    - younger men at baseline: age <61
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")

#' # INTERVAL phenotypes ####
#' ***
#' ## Load data 
data = fread(INTERVAL_phenoptype)
head(data)

map = fread(INTERVAL_pheno_geno_map)
head(map)
map[,table(is.na(Affymetrix_gwasQC_bl),is.na(Affymetrix_gwasQC_24m))]
map[!is.na(Affymetrix_gwasQC_bl) & !is.na(Affymetrix_gwasQC_24m)]
map[, ID2 := Affymetrix_gwasQC_bl]
map[is.na(Affymetrix_gwasQC_bl),ID2 := Affymetrix_gwasQC_24m]
map[is.na(ID2)]
map = map[!is.na(ID2)]

matched = match(data$identifier,map$identifier)
data[,genID := map[matched,ID2]]

genPCs = fread(INTERVAL_SNP_PCs) 
table(data$genID %in% genPCs$ID)
genRelation = fread(INTERVAL_SNP_relation)
table(data$genID %in% genRelation$V1)

#' Filter data for genetic data available & white ancestry & no relationship 
data = data[!is.na(genID),]

data[,table(ethnicPulse)]
data[,table(ethnicPulse,ethnic_bl)]
data = data[ethnicPulse %in% c("Eng/W/Scot/NI/Brit","Other White backgrnd","White Irish"),]

data[,table(genID %in% genRelation$V1)]
data = data[!is.element(genID,genRelation$V1),]

#' Quick age check
hist(data$agePulse)

#' Quick sex check
data[,table(sexPulse)]
data[,table(sexPulse,agePulse<50)]

#' # Menopause and lipid medication ####
#' ***
#' "Have you had your menopause?"
#' 
#' - 1 = yes, 
#' - 2 = no, 
#' - 3 = not sure – had hysterectomy, 
#' - 4 = not sure – other reason, 
#' - 999 = don’t know / prefer not to answer
#' 
#' "Do you currently take any medications that control your cholesterol?"
#' 
#' - 1 = yes, 
#' - 2 = no, 
#' - 999 = don’t know
#' 
#' Check menopause
data[,table(sexPulse,menopause_bl)]

#' Remove men who might be trans
filt = data$sexPulse == 1 & !is.na(data$menopause_bl)
table(filt)
data = data[!filt,]

#' Simplify to 2 categories: yes/no, NA if not these two
data[is.element(menopause_bl,c(3,4,999)),menopause_bl :=NA]
data[,table(sexPulse,menopause_bl)]
data[,table(sexPulse,is.na(menopause_bl))]

#' Check lipid medication
data[,table(med_lipid_24m)]
data[,table(med_lipid_48m)]

#' Simplify to 2 categories: yes/no, NA if not these two
data[is.element(med_lipid_24m,c(3,4,999)),med_lipid_24m :=NA]
data[is.element(med_lipid_48m,c(3,4,999)),med_lipid_48m :=NA]

#' create strata
data[sexPulse==2 & menopause_bl==2 & agePulse <=45, strata := "premenopausal"]
data[sexPulse==1 & agePulse <=45, strata := "young men"]
data[sexPulse==2 & menopause_bl==1 & agePulse >50, strata := "postmenopausal"]
data[sexPulse==1 & agePulse >50, strata := "older men"]
data[sexPulse==2, strata2 := "women"]
data[sexPulse==1, strata2 := "men"]

data[,table(is.na(strata))]
data[,table(is.na(strata2))]

#' # Reformat data ####
#' ***
#' I want the data in the long format, with the following columns: 
#' 
#' - fixed over time: ID, sex, ethnic, menopausse status at baseline
#' - changing over time: age, attendance date, medication, TC value, visit (bl, 24m or 48m)
#' - helper columns (to be deleted at the end): age at baseline, date at baseline
#' 
myNames1 = names(data)[c(1,22,2,4,10,23,    5,3,13,12)]
myNames2 = names(data)[c(1,22,2,4,10,23,    16,3,15,12,   5)]
myNames3 = names(data)[c(1,22,2,4,10,23,    19,3,21,18,   5)]

colsOut<-setdiff(colnames(data),myNames1)
data1 = copy(data)
data1[,get("colsOut"):=NULL]
setcolorder(data1,myNames1)
names(data1) = c("ID","genID","sex","ethnic","meno","strata","date","age","CHOL","med_lipid")
data1[,table(is.na(CHOL))]
data1 = data1[!is.na(CHOL),]

data1[,day := as.numeric(substr(date,1,2))]
data1[,month := substr(date,3,5)]
monthTable = data.table(months = unique(data1$month))
setorder(monthTable,months)
monthTable[,NR := c(4,8,12,2,1,7,6,3,5,11,10,9)]
setorder(monthTable,NR)
matched = match(data1$month,monthTable$months)
data1[,month2 := monthTable[matched,NR]]
data1[,year := as.numeric(substr(date,6,9))]
data1[,date := paste(year,month2,day,sep = "-")]
data1[,date := as.Date(date)]

colsOut<-setdiff(colnames(data),myNames2)
data2 = copy(data)
data2[,get("colsOut"):=NULL]
setcolorder(data2,myNames2)
names(data2) = c("ID","genID","sex","ethnic","meno","strata","date","age","CHOL","med_lipid","date_bl")
data2[,table(is.na(CHOL))]
data2 = data2[!is.na(CHOL),]

data2[,day := as.numeric(substr(date_bl,1,2))]
data2[,month := substr(date_bl,3,5)]
matched = match(data2$month,monthTable$months)
data2[,month2 := monthTable[matched,NR]]
data2[,year := as.numeric(substr(date_bl,6,9))]
data2[,date_bl := paste(year,month2,day,sep = "-")]
data2[,date_bl := as.Date(date_bl)]

data2[,day := as.numeric(substr(date,1,2))]
data2[,month := substr(date,3,5)]
matched = match(data2$month,monthTable$months)
data2[,month2 := monthTable[matched,NR]]
data2[,year := as.numeric(substr(date,6,9))]
data2[,date := paste(year,month2,day,sep = "-")]
data2[,date := as.Date(date)]

data2[,diff := as.numeric(date -date_bl)/365.25]
data2[,age := age + diff ]

colsOut<-setdiff(colnames(data),myNames3)
data3 = copy(data)
data3[,get("colsOut"):=NULL]
setcolorder(data3,myNames3)
names(data3) = c("ID","genID","sex","ethnic","meno","strata","date","age","CHOL","med_lipid","date_bl")
data3[,table(is.na(CHOL))]
data3 = data3[!is.na(CHOL),]

data3[,day := as.numeric(substr(date_bl,1,2))]
data3[,month := substr(date_bl,3,5)]
matched = match(data3$month,monthTable$months)
data3[,month2 := monthTable[matched,NR]]
data3[,year := as.numeric(substr(date_bl,6,9))]
data3[,date_bl := paste(year,month2,day,sep = "-")]
data3[,date_bl := as.Date(date_bl)]

data3[,day := as.numeric(substr(date,1,2))]
data3[,month := substr(date,3,5)]
matched = match(data3$month,monthTable$months)
data3[,month2 := monthTable[matched,NR]]
data3[,year := as.numeric(substr(date,6,9))]
data3[,date := paste(year,month2,day,sep = "-")]
data3[,date := as.Date(date)]

data3[,diff := as.numeric(date -date_bl)/365.25]
data3[,age := age + diff ]

data1[,scan := "0m"]
data2[,scan := "24m"]
data3[,scan := "48m"]

myTab = rbind(data1[,c(1:10,15)],data2[,c(1:10,17)],data3[,c(1:10,17)])
myTab[,table(strata)]
myTab[,table(is.na(strata))]
test = myTab[,.N,by=c("ID","strata")]
test[,table(strata,N)]


ggp0  = ggplot(myTab[!is.na(strata),], aes(x=scan, y=CHOL, group=ID, col=strata, fill=strata)) +
  facet_wrap(~strata,scales="free") +
  geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point()+
  labs(x="Scans", y="TC (mmol/L)", title = "All samples") +
  theme_classic() +
  theme(legend.position = "none")
ggp0

ggsave(filename = "../results/INTERVAL/01_Trajectories_scan.png",plot = ggp0)

ggp1  = ggplot(myTab[!is.na(strata),], aes(x=age, y=CHOL, group=ID, col=strata, fill=strata)) +
  facet_wrap(~strata,scales="free") +
  geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
  geom_line(aes(alpha=0.01)) + 
  geom_point()+
  labs(x="Age (years)", y="TC (mmol/L)", title = "All samples") +
  theme_classic() +
  theme(legend.position = "none")
ggp1

ggsave(filename = "../results/INTERVAL/01_Trajectories_age.png",plot = ggp1)

#' # Checks ####
#' *** 
#' Check sample ID and genetic ID
#' 
psam = fread(paste0(INTERVAL_SNP_data,"impute_dedup_1_interval.psam"))
head(psam)
head(myTab)

table(psam$IID %in% myTab$genID)
table(myTab$genID %in% psam$IID)

#' Check 2: check range of TC values
hist(myTab$CHOL)
summary(myTab$CHOL)

#' Check 3: lipid lowering medication coding
#' 
myTab[,table(scan,is.na(med_lipid))]
myTab[,table(scan,(med_lipid))]

#' I want 0 to code for no treatment or no information, and 1 as treated with lipid-lowering medication. 
myTab[is.na(med_lipid),med_lipid := 0]
myTab[med_lipid==2,med_lipid := 0]

#' Check 4: sample sizes
test = myTab[,.N,by=c("ID","strata")]
test[,table(N,strata)]
test[,table(strata)]
test[N>1,table(strata)]

#' 
#' - older men:               7,798 (with 2 or more time points: 5,894)
#' - postmenopausal women:    4,362 (with 2 or more time points: 3,443)
#' - premenopausal women:    10,215 (with 2 or more time points: 5,523)
#' - younger men:             8,106 (with 2 or more time points: 4,799)
#'  
test = myTab[,.N,by=c("ID","sex")]
test[,table(N,sex)]
test[,table(sex)]
test[N>1,table(sex)]

#'
#' - men:               18,159 (with 2 or more time points: 12,359)
#' - women:             18,534 (with 2 or more time points: 11,737)
#'
#' # Save data ####
#' ***
myTab[,age2 := age^2]
myTab[,plot(age,age2)]

matched = match(myTab$genID,genPCs$ID)
table(is.na(matched))
myTab = cbind(myTab,genPCs[matched,c(2:12)])

save(myTab, file = paste0(data_QC,"/INTERVAL/01_Prep_TrajGWAS.RData"))
write.table(myTab[,unique(genID)],file = paste0(data_QC,"/INTERVAL/01_Prep_TrajGWAS_SampleList.txt"), 
            col.names = F, row.names = F, quote = F)

#' # Create input for trajGWAS ####
#' ***
#' I need to create the csv for the trajGWAS input and the sample list for the PLINK extraction. In addition, I want some basic description per sub-population as well (sample size, mean TC, sd TC, number of time points, average age when statin was first described, average age at first time point, average age at last time point)
#' 
groups = c("men","women","older men","postmenopausal","premenopausal","young men")
groups2 = c("men","women","old","post","pre","young")

dumTab = foreach(i = 1:length(groups))%do%{
  #i=1
  data = copy(myTab)
  if(i<3) data = data[sex==i,]
  if(i>2) data = data[strata == groups[i]]
  test = data[,.N,by=c("ID")]
  data = data[ID %in% test[N>1,ID]]
  
  # Save group specific data
  write.csv(data,
            file = paste0(data_QC,"/INTERVAL/01_Prep_TrajGWAS_",groups2[i],".csv"), 
            row.names = FALSE)
  write.table(data[,unique(genID)],
              file = paste0(data_QC,"/INTERVAL/01_Prep_TrajGWAS_SampleList_",groups2[i],".txt"), 
              col.names = F, row.names = F, quote = F)
  dum1 = data[,min(age),by=ID]
  dum2 = data[,max(age),by=ID]
  dum3 = data[,.N,by=ID]
  setnames(dum1,"V1","minAge")
  dum1[,maxAge := dum2$V1]
  dum1[,timeDiff := maxAge - minAge]
  dum1[,timePoints := dum3$N]
  dum4 = data[med_lipid==1,min(age),by = ID]
  
  res1 = data.table(setting = groups2[i],
                    sampleSize = dim(dum1)[1],
                    meanTC = data[,mean(CHOL)],
                    sdTC = data[,sd(CHOL)],
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
save(res,file="../results/INTERVAL/01_descriptiveTable.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

