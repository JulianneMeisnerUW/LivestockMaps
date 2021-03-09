rm(list = ls())
set.seed(1154)
library(rgdal)
library(sp)
library(splines)
library(raster)
library(rgeos)
library(dplyr)
library(stringr)
library(robustHD)
library(parallel)

source("gformula/gformula_functions.R")

folderData<-paste0("gform_results/Malawi/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

n.cores<-8

type<-"rhod"
species<-"cattle"
mediator<-"NDVI"

HAT<-readOGR(paste0("Outcome_data/Malawi/HAT_final"), layer="HAT_finalPO")
colnames(HAT@data)<-c("cluster", "year","focus","cases","density_med_c","density_sd_c", "density_q025_c","density_q975_c",   
                              "density_med_p","density_sd_p","density_q025_p","density_q975_p","source","prox1km",  "intercept",       
                              "WorldPop","LandScan","urban","nightlights","wealth","wealth_sd","wealth_q025","wealth_975",
                                "travel_hrs", "protected", "district", "disaster",         
                              "dis_num","dis_deaths","NDVI","LST","elevation", "st.cluster","dist_clust", "pigs_POsum", "pigs_POmean", "cattle_POsum", "cattle_POmean")

'%ni%' <- Negate('%in%')

HAT@data$st.cluster<-1:nrow(HAT@data)

HAT@data$disaster<-ifelse(HAT@data$disaster=="None", 0,1)

#Set up predictors
Y<-c("wealth", "NDVI",  "cases", "livestock")
X<-list(c("livestock_lag1", "wealth_lag"),
          c("livestock_lag1"),
          c("NDVI_lag",  "PO_lag2", "livestock_lag2", "offset(log(WorldPop))"),
          c("wealth", "NDVI", "livestock_lag1"))
family<-c("gaussian", "gaussian", "poisson", "gaussian")

lag1vars<-c("NDVI", "LST", "wealth", "livestock", "disaster")
lag1names<-paste0(lag1vars, "_lag")
lag1names[which(lag1names=="livestock_lag")]<-"livestock_lag1"

lag2vars<-c("livestock", "PO")
lag2names<-paste0(lag2vars, "_lag2")

cat.vars<-c("disaster", "disaster_lag", "protected")

if(length(Y)!=length(X)){
  warning("Mismatch between predictors and time-varying confounders")
}

HAT@data$cluster<-rep(NA, nrow(HAT@data))

years<-unique(HAT@data$year)
 
for (year in years){
    sub<-HAT[which(HAT@data$year==year),]
    idx<-sub@data$st.cluster
    cluster<-1:nrow(sub@data)
    HAT@data$cluster[which(HAT@data$st.cluster%in%idx)]<-cluster
}

#Remove clusters with 0 pop for WorldPop in any year; only keep clusters with at least 1 cases reported over the study period
zero_list<-c()
cases_list<-c()

for(year in years){
    sub<-HAT[which(HAT@data$year==year),]
    zero<-which(sub@data$WorldPop==0)
    zero_clus<-sub@data$cluster[zero]
    zero_list<-c(zero_list, zero_clus)
}
  
zero<-unique(zero_list) 

HAT<-HAT[which(HAT@data$cluster%ni%zero),]

#Re-set cluster
for (year in years){
  sub<-HAT[which(HAT@data$year==year),]
  idx<-sub@data$st.cluster
  cluster<-1:nrow(sub@data)
  HAT@data$cluster[which(HAT@data$st.cluster%in%idx)]<-cluster
}

HAT<-HAT@data
  
#Collapse on 5 year intervals: 2000-2004, 2005-2009, 2010-2014, 2015-2018
HAT$interval<-ifelse(HAT$year<2005,1, 
                     ifelse(HAT$year>=2005&HAT$year<2010,2,
                            ifelse(HAT$year>=2010&HAT$year<2015, 3, 4)))

tm.vec<-c("a1m0", "a0m0", "a1m1")

results.mean=list()
results.total=list()
results.inc=list()

nSamp=100

for(i in 1:nSamp){
  res<-gformMed(tm.vec, HAT, years, type, species, mediator, 
              Y.vars=Y, X.vars=X, family, lag1vars, 
              lag2vars, lag1names, 
              lag2names, cat.vars)
  
  results.mean=c(results.mean, list(res$mean))
  results.total=c(results.total, list(res$total))
  results.inc=c(results.inc, list(res$inc))
  print(i)
}
    
if(species=="cattle"){
  save(file=paste0("gformula/Malawi/results/gform_cattle_mean_med_", mediator, "_", type, ".Rdata"), results.mean)
  save(file=paste0("gformula/Malawi/results/gform_cattle_total_med_",mediator, "_",type, ".Rdata"), results.total)
  save(file=paste0("gformula/Malawi/results/gform_cattle_inc_med",mediator, "_",type, ".Rdata"), results.inc)
}else{
  save(file=paste0("gformula/Malawi/results/gform_pigs_mean_med_",mediator, "_",type, ".Rdata"), results.mean)
  save(file=paste0("gformula/Malawi/results/gform_pigs_total_med_",mediator, "_",type, ".Rdata"), results.total)
  save(file=paste0("gformula/Malawi/results/gform_pigs_inc_med_",mediator, "_",type, ".Rdata"), results.inc)
}

