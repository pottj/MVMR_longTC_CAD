#' ---
#' title: "Get exposure data"
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

#' # UKB Phenotype Data ####
#' ***
#' I want TC data from the BSU UKB data and from the GP records. 
#' 
#' ## BSU data set ####
#' 
#' - totchol          --> 30690
#' 
#' In addition, I want the following columns/covariables: 
#' 
#' - ID                 --> eid
#' - Biological sex     --> 31
#' - Assessment date    --> 53
#' - Menopausal status  --> 2724
#' - Medication         --> 20003
#' - Ancestry           --> 21000
#' - Age                --> 21022
#' - Genotype bach      --> 22000
#' - Genetic PCs        --> 22009
#' - Kinship            --> 22021

myTab_20 = fread(UKB_phenotypes, header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:18506]

#' get relevant parameters: 
#'  
exposure = c("f.30690.0.0","f.30690.1.0")
table(is.element(exposure,myAnnot$colNm))

covars = c("f.eid", "f.31.0.0","f.53.0.0","f.53.1.0","f.2724.0.0","f.2724.1.0",
           paste0("f.20003.1.",c(0:47)),paste0("f.20003.0.",c(0:47)),
           "f.21000.0.0","f.21022.0.0","f.22000.0.0",paste0("f.22009.0.",1:10), "f.22021.0.0")
table(is.element(covars,myAnnot$colNm))

myAnnot = myAnnot[colNm %in% c(exposure,covars),]

x = myAnnot[,colNR]
myTab_cross <- fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
names(myTab_cross)
names(myTab_cross) = c("ID","sex","date_init","date_FU","menopause_init","menopause_FU",
                       paste("medication_init",1:48,sep="_"),
                       paste("medication_FU",1:48,sep="_"),
                       "ancestry","age","GenotypeBatch",
                       paste("PC",1:10,sep="_"),
                       "kinship","TC_init","TC_FU")

save(myTab_cross, file = paste0(data_QC,"01_Prep_01_UKB_unfiltered.RData"))
load(paste0(data_QC,"01_Prep_01_UKB_unfiltered.RData"))

#' Filter for White British People without any kinship
myTab = copy(myTab_cross) 
myTab = myTab[ancestry == 1001,]
myTab = myTab[kinship == 0,]

save(myTab, file = paste0(data_QC,"01_Prep_01_UKB_filtered_Ancestry_Kinship.RData"))
load(paste0(data_QC,"01_Prep_01_UKB_filtered_Ancestry_Kinship.RData"))

#' Get medication: 
#' 
#' For lipid lowering medication I use C10 - lipid modifying agents
#' 
codingTable = data.table(read_xlsx(MedicationCoding,sheet=1))
lipidLowering = codingTable[grepl("C10",ATC_Code),Coding]

myMeds = names(myTab)[grep("medication_init",names(myTab))]
myTab[,lipidLow_init := 0]
for(i in 1:length(myMeds)){
  #i=1
  myTab[get(myMeds[i]) %in% lipidLowering,lipidLow_init := 1]
}
myTab[,get("myMeds"):=NULL]

myMeds = names(myTab)[grep("medication_FU",names(myTab))]
myTab[,lipidLow_FU := 0]
for(i in 1:length(myMeds)){
  #i=1
  myTab[get(myMeds[i]) %in% lipidLowering,lipidLow_FU := 1]
}
myTab[,get("myMeds"):=NULL]
dim(myTab)

#' Check for consent
ToExclude = fread(gsub("ukb672224.tab","withdraw98032_19.txt",UKB_phenotypes))
table(is.element(myTab$ID,ToExclude$V1))

#' Filter for no NA in TC baseline data
myTab = myTab[!is.na(TC_init),]
myTab[,table(is.na(TC_FU))]

#' Check menopausal status - restrict to women who are sure (0 - no, 1 - yes)
myTab[,table(sex,menopause_init)]
# myTab = myTab[sex==1 | (sex==0 & menopause_init %in% c(0,1)),]
myTab[,table(sex,menopause_init)]
myTab[,table(menopause_init,menopause_FU)]
myTab[sex==1, group := "men"]
myTab[sex==0 & menopause_init==0, group := "premenopausal"]
myTab[sex==0 & menopause_init==1, group := "postmenopausal"]
myTab[group=="postmenopausal" & menopause_FU!=1, group := NA]

#' Check genotyping batch
myTab[,table(GenotypeBatch)]
myTab[,Array := "Axiom"]
myTab[GenotypeBatch<0,Array := "BiLEVE"]
myTab[,table(Array,sex)]

# save filtered data
save(myTab, file = paste0(data_QC,"01_Prep_01_UKB_filtered_Ancestry_Kinship_Meds.RData"))

#' ## GP data set ####
loaded1 = load(UKB_EHR)
myTab2 = get(loaded1)
setDT(myTab2)
myKey = fread(UKB_EHR_IDs)

myTab2 = myTab2[idno %in% myKey$Adiposity_sample_ID,]
matched = match(myTab2$idno,myKey$Adiposity_sample_ID)
myTab2[,BSU_ID := myKey[matched, BSUid]]
myTab2[,.N,by=exposure_type]
myTab2 = myTab2[exposure_type %in% c("tchol","statins"),]

#' Okay, I really need the tchol data, but for statins, I just want to know yes/no
TC_check = myTab2[exposure_type == "tchol",.N,by = c("BSU_ID")]
Statin_check = myTab2[exposure_type == "statins",.N,by = c("BSU_ID")]

#' Check overlap
table(is.element(myTab$ID,TC_check$BSU_ID))
table(is.element(TC_check$BSU_ID,myTab$ID))

myTab = myTab[ID %in% TC_check$BSU_ID, ]
myTab2 = myTab2[BSU_ID %in% myTab$ID, ]

myTab2_statins = copy(myTab2)[exposure_type == "statins",]
myTab2 = myTab2[exposure_type == "tchol",]

#' # Merge data ####
#' ***
#' How do I want my data to look like? 
#' 
#' Longitudinal data: ID - date - TC - age - medication - source - dumID
#' 
myTab2[,dumID := paste(BSU_ID,exposure_date, sep="_")]
myTab2[,table(duplicated(dumID))]
myDups = myTab2[duplicated(dumID),dumID]
myTab3 = copy(myTab2)
myTab3 = myTab3[dumID %in% myDups,]
Dup_check = myTab3[,.N, by = c("dumID","exposure_value")]
table(duplicated(Dup_check$dumID))

bad_dups = Dup_check[duplicated(dumID),dumID]
myTab2 = myTab2[!is.element(dumID,bad_dups)]
myTab2 = myTab2[!duplicated(dumID)]

#' Add medication
myTab2_statins[,dumID := paste(BSU_ID,exposure_date, sep="_")]
Statins_minAge = myTab2_statins[,min(exposure_age),by=BSU_ID]
myTab2[, lipLowMed := 0]

counter = seq(1,30000,1000)
for(i in 1:dim(Statins_minAge)[1]){
  #i=1
  if(i %in% counter)message("Working on i=",i)
  myRow = Statins_minAge[i,]
  myTab2[BSU_ID == myRow$BSU_ID & exposure_age>=myRow$V1, lipLowMed := 1]
}

myTab2[,table(lipLowMed)]
myTab2[,source := "GP"]

save(myTab2, file=paste0(data_QC,"01_Prep_01_UKB_GP.RData"))

#' Okay, now I add my UKB data set

myTab[!is.na(date_FU),time_diff := date_FU - date_init]
myTab[!is.na(date_FU),time_diff_years := time_diff/365]
myTab[!is.na(date_FU),age_FU := age + time_diff_years]
myTab[,dumID_init := paste(ID,date_init,sep="_")]
myTab[,dumID_FU := paste(ID,date_FU,sep="_")]
myTab[,BSU_ID := ID]

myTab4 = copy(myTab)
myTab4[,source := "UKB_Baseline"]
myNames = c("ID","date_init","TC_init","age","BSU_ID","dumID_init","lipidLow_init","source")
colsOut = setdiff(colnames(myTab4),myNames)
myTab4[,get("colsOut"):=NULL]
setcolorder(myTab4,myNames)
dim(myTab4)
table(is.element(myTab4$dumID_init,myTab2$dumID))
myTab2 = myTab2[!is.element(dumID, myTab4$dumID_init)]

myTab5 = copy(myTab)
myTab5[,source := "UKB_FollowUp"]
myNames = c("ID","date_FU","TC_FU","age_FU","BSU_ID","dumID_FU","lipidLow_FU","source")
colsOut = setdiff(colnames(myTab5),myNames)
myTab5[,get("colsOut"):=NULL]
setcolorder(myTab5,myNames)
myTab5 = myTab5[!is.na(TC_FU),]
table(is.element(myTab5$dumID_FU,myTab2$dumID))
myTab2 = myTab2[!is.element(dumID, myTab5$dumID_init)]

names(myTab4) = c("idno","exposure_date","exposure_value","exposure_age","BSU_ID","dumID","lipLowMed","source")
names(myTab5) = c("idno","exposure_date","exposure_value","exposure_age","BSU_ID","dumID","lipLowMed","source")

#' Harmonize the class of the date columns
myTab2[,date := gsub(".*_","",dumID)]
myTab4[,date := gsub(".*_","",dumID)]
myTab5[,date := gsub(".*_","",dumID)]

myTab2[,exposure_date:=NULL]
myTab4[,exposure_date:=NULL]
myTab5[,exposure_date:=NULL]

myTab6 = rbind(myTab2, myTab4, myTab5, fill=T)
setorder(myTab6,dumID)

#' I have the following problem: in the UKB, I only have age rounded to years, e.g. 48. In the GP data, age is rounded to 3 digits, e.g. 48.784. When I sort by chronological date, the age is not always increasing. I will simply exclude these ages, as I have no clever way to correct for this problem (how to calculate the birth date from age in R).
#' 
myTab6[,flag:=T]
counter = seq(1,700000,25000)
for(i in 2:dim(myTab6)[1]){
  # i=2
  if(i %in% counter)message("Working on i=",i)
  if(myTab6[i,BSU_ID] == myTab6[i-1,BSU_ID]){
    check = myTab6[i,exposure_age] - myTab6[i-1,exposure_age]
    if(sign(check) == 1){ 
      myTab6[i,flag:=T]
    }else{
      myTab6[i,flag:=F]
    }
  }else{
    myTab6[i,flag:=T]
  }
}

myTab6[,table(flag)]
myTab6 = myTab6[flag==T,]

# okay, now recheck that the lipid lowering medication is still homogenous (1 after first mentioned in GP file / UKB)
Statins_minAge2 = myTab6[lipLowMed==1,min(exposure_age),by=BSU_ID]

counter = seq(1,31000,1000)
for(i in 1:dim(Statins_minAge2)[1]){
  #i=1
  if(i %in% counter)message("Working on i=",i)
  myRow = Statins_minAge2[i,]
  myTab6[BSU_ID == myRow$BSU_ID & exposure_age>=myRow$V1, lipLowMed := 1]
}

myTab6[,table(lipLowMed)]

#' Filter for at least 3 time points
TimePoint_Check = myTab6[,.N,BSU_ID]
myTab6 = myTab6[BSU_ID %in% TimePoint_Check[N>2,BSU_ID],]
TimePoint_Check2 = myTab6[,.N,BSU_ID]
summary(TimePoint_Check2$N)

#' Filter for max. median + 6*median
myTab6 = myTab6[BSU_ID %in% TimePoint_Check2[N<43,BSU_ID],]
TimePoint_Check2 = myTab6[,.N,BSU_ID]

#' Get all the covariates from the UKB baseline table
myTab7 = copy(myTab)
myTab7 = myTab7[ID %in% myTab6$BSU_ID,]
myTab7[,table(sex)]
myTab7[,table(group)]

#' When I read the tables correctly: 
#' 
#' - there are 67,285 samples: 
#'    - 35,726 men
#'    - 38,052 women
#'    - 25,192 women postmenopausal (not yet filtered for age >= baseline)
#'    -  6,341 women premenopausal  (not yet filtered for age <= baseline)
#'    
#' # Some data QC ####
#' ***
#' I want to check the HDLC and TC data and filter outliers if necessary
#' 
#' In UKB data: mmol/l, values between 0.601 and 15.46
#' 
hist(myTab6$exposure_value)
mean(myTab6$exposure_value,na.rm=T)
sd(myTab6$exposure_value,na.rm=T)
summary(myTab6$exposure_value,na.rm=T)
mean(myTab6$exposure_value,na.rm=T) + 8 * sd(myTab6$exposure_value,na.rm=T)
table(myTab6$exposure_value>mean(myTab6$exposure_value,na.rm=T) + 8 * sd(myTab6$exposure_value,na.rm=T))

#' # Save data ####
#' ***
save(myTab6, myTab7, file = paste0(data_QC,"/01_Prep_01_UKB_GP_TC_GLGC.RData"))
write.table(myTab7$ID,file = paste0(data_QC,"/01_Prep_01_SampleList_TC_GLGC.txt"), 
            col.names = F, row.names = F, quote = F)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

