#' ---
#' title: "Test non-genetic model (subset)"
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
#' Here I want to test the non-genetic model. I will adjust for age and lipid lowering medication, and will include random effects for intercept, slope, and variability. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile.R")
.libPaths()

#' # All samples ####
#' ***
load(paste0(data_QC,"/01_Prep_02_UKB_GP_noLipidLoweringSubSet.RData"))
names(myTab6)
setnames(myTab6,"exposure_value","TC")
setnames(myTab6,"exposure_age","age")

matched = match(myTab6$BSU_ID,myTab7$ID)
myTab6[,sex := myTab7[matched,sex]]
myTab6[sex==0,sex:=2]

#' ## Men ####
data1 = copy(myTab6)
data1 = data1[sex==1,]
data1[,ID := as.factor(BSU_ID)]

time1<-Sys.time()
resulta = tryCatch({
  mod = gamlss(TC ~ age + re(random=~1+age|ID),   
               sigma.formula = ~ age + re(random=~1|ID), 
               data = na.omit(data1), family = "NO")
  dummy = summary(mod)
  
  res = data.table(model = rep("R3",each=2),
                   estimate = c("mu_age","sigma_age"),
                   beta = dummy[c(2,4),1],
                   SE = dummy[c(2,4),2],
                   tval = dummy[c(2,4),3],
                   pval = dummy[c(2,4),4], 
                   comment = rep("OK",2))
  
}, error = function(e) {
  res = data.table(model = rep("R3",each=2),
                   estimate = c("mu_age","sigma_age"),
                   comment = rep("GAMLSS failed to converge",2))
  return(res)
})

message("Time for model with 3 random effects: " ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

time2<-Sys.time()
resultb = tryCatch({
  mod = gamlss(TC ~ age + re(random=~1|ID),   
               sigma.formula = ~ age + re(random=~1|ID), 
               data = na.omit(data1), family = "NO")
  dummy = summary(mod)
  
  res = data.table(model = rep("RIRI",each=2),
                   estimate = c("mu_age","sigma_age"),
                   beta = dummy[c(2,4),1],
                   SE = dummy[c(2,4),2],
                   tval = dummy[c(2,4),3],
                   pval = dummy[c(2,4),4], 
                   comment = rep("OK",2))
  
}, error = function(e) {
  res = data.table(model = rep("RIRI",each=2),
                   estimate = c("mu_age","sigma_age"),
                   comment = rep("GAMLSS failed to converge",2))
  return(res)
})

message("Time for model with 2 random intercepts: " ,round(difftime(Sys.time(),time2,units = "mins"),3)," minutes")

time3<-Sys.time()
resultc = tryCatch({
  mod = gamlss(TC ~ age + re(random=~1+age|ID),   
               sigma.formula = ~ age, 
               data = na.omit(data1), family = "NO")
  dummy = summary(mod)
  
  res = data.table(model = rep("RIRS",each=2),
                   estimate = c("mu_age","sigma_age"),
                   beta = dummy[c(2,4),1],
                   SE = dummy[c(2,4),2],
                   tval = dummy[c(2,4),3],
                   pval = dummy[c(2,4),4], 
                   comment = rep("OK",2))
  
}, error = function(e) {
  res = data.table(model = rep("RIRS",each=2),
                   estimate = c("mu_age","sigma_age"),
                   comment = rep("GAMLSS failed to converge",2))
  return(res)
})

message("Time for model with random intercept and random slope: " ,round(difftime(Sys.time(),time3,units = "mins"),3)," minutes")

time4<-Sys.time()
resultd = tryCatch({
  mod = gamlss(TC ~ age + re(random=~1|ID),   
               sigma.formula = ~ age, 
               data = na.omit(data1), family = "NO")
  dummy = summary(mod)
  
  res = data.table(model = rep("RI",each=2),
                   estimate = c("mu_age","sigma_age"),
                   beta = dummy[c(2,4),1],
                   SE = dummy[c(2,4),2],
                   tval = dummy[c(2,4),3],
                   pval = dummy[c(2,4),4], 
                   comment = rep("OK",2))
  
}, error = function(e) {
  res = data.table(model = rep("RI",each=2),
                   estimate = c("mu_age","sigma_age"),
                   comment = rep("GAMLSS failed to converge",2))
  return(res)
})

message("Time for model with random intercept: " ,round(difftime(Sys.time(),time4,units = "mins"),3)," minutes")

#' # Combine and save ####
#' ***
res = rbind(resulta,resultb,resultc,resultd, use.names=T, fill=T)
res[estimate=="mu_age",]
res[estimate=="mu_lipLowMed",]
res[estimate=="sigma_age",]
res[estimate=="sigma_lipLowMed",]

res[,sampleSize := dim(data1)[1]]
res[,group := "men"]
res[,setting := "subset"]

save(res, file = "../results/01_Prep_07_NonGeneticAssociations_men_subset.RData")

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
