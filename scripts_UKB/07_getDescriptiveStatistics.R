#' ---
#' title: "Get descriptive data"
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

#' # Load UKB Data ####
#' ***
#' 
load("../results/UKB/01_descriptiveTable.RData") 
load(paste0(data_QC,"/UKB/01_Prep_TrajGWAS.RData"))

myFiles = list.files(path=paste0(data_QC,"/UKB/"),pattern = ".csv")

dumTab = foreach(i = 1:length(myFiles))%do%{
  #i=1
  data = fread(paste0(data_QC,"/UKB/",myFiles[i]))
  trait = gsub("01_Prep_TrajGWAS_","",myFiles[i])
  trait = gsub(".csv","",trait)
  data[,group2 := trait]
  data
}

myTab = rbindlist(dumTab)
myTab[,.N,by = group2]

#' # Save tables ####
#' ***
tosave4 = data.table(data = c("res"), 
                     SheetNames = c("Merged"))
excel_fn = paste0("../results/UKB/07_descriptiveTables.xlsx")
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

write.table(res, file = paste0("../results/UKB/07_descriptiveTables.txt"),
            col.names = T,row.names = F,quote = F,sep="\t")

#' # Trajectories ####
#' ***
plotData = copy(myTab)[group2 %in% c("old","young","post","pre")]
plotData[,table(sex,group)]
plotData[,group3 := paste(sex,group,sep=" - ")]
plotData[,group3 := gsub("1 - ","males - ",group3)]
plotData[,group3 := gsub("2 - ","females - ",group3)]
plotData[,table(group2,group3)]
plotData = plotData[age<85,]
plotData = plotData[statin==0,]

#' Simple plots per baseline group
p1 = ggplot(data = plotData, aes(x = age, y = TC, group = ID, col=group3, fill=group3))
p1 + geom_line() + 
  facet_wrap(. ~ group3,scales = "free_x") +
  stat_smooth(aes(group = 1)) + 
  #stat_summary(aes(group = 1),geom = "point", fun.y = mean, shape = 17, size = 3) + 
  labs(x="Age (years)", y="TC (mmol/L)") +
  theme_classic()+
  guides(col="none",fill="none",shape="none")

ggsave(filename = paste0("../results/UKB/07_trajectories/Plot1a_allLines.png"),width = 12,height = 8)

p1 + #geom_line() + 
  facet_wrap(. ~ group3, scales = "free_x") +
  stat_smooth(aes(group = 1)) + 
  #stat_summary(aes(group = 1),geom = "point", fun.y = mean, shape = 17, size = 3) + 
  labs(x="Age (years)", y="TC (mmol/L)") +
  theme_classic() +
  guides(col="none",fill="none",shape="none")

ggsave(filename = paste0("../results/UKB/07_trajectories/Plot1b_summaryPerGroup.png"),width = 12,height = 8)

plotData = copy(myTab)[group2 %in% c("men","women")]
plotData[,table(sex,group)]
plotData[,group3 := sex]
plotData[,group3 := gsub("1","males",group3)]
plotData[,group3 := gsub("2","females",group3)]
plotData[,table(group2,group3)]
plotData = plotData[age<85,]
plotData = plotData[statin==0,]

#' Simple plots per baseline group
p2 = ggplot(data = plotData, aes(x = age, y = TC, group = ID, col=group3, fill=group3))
p2 + geom_line() + 
  facet_wrap(. ~ group3,scales = "free_x") +
  stat_smooth(aes(group = 1)) + 
  #stat_summary(aes(group = 1),geom = "point", fun.y = mean, shape = 17, size = 3) + 
  labs(x="Age (years)", y="TC (mmol/L)") +
  theme_classic()+
  guides(col="none",fill="none",shape="none")

ggsave(filename = paste0("../results/UKB/07_trajectories/Plot2a_allLines.png"),width = 12,height = 8)

p2 + #geom_line() + 
  facet_wrap(. ~ group3, scales = "free_x") +
  stat_smooth(aes(group = 1)) + 
  #stat_summary(aes(group = 1),geom = "point", fun.y = mean, shape = 17, size = 3) + 
  labs(x="Age (years)", y="TC (mmol/L)") +
  theme_classic() +
  guides(col="none",fill="none",shape="none")

ggsave(filename = paste0("../results/UKB/07_trajectories/Plot2b_summaryPerSex.png"),width = 12,height = 8)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

