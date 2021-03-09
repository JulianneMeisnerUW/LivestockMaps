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
source("codes/gformula_functions.R")

folderData<-paste0("gform_results/DRC/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

type<-"rhod"
species<-"cattle"
test<-"yes"

HAT<-readOGR("Data/Outcome_data/DRC/HAT_final", layer="HAT_finalPO", p4s=proj4string(humdata))
colnames(HAT@data)<-c("cluster", "interaction", "year","focus","cases","density_med_c", "density_sd_c",  
                      "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                      "density_q025_p", "density_q975_p", "source", "intercept",
                      "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                      "wealth_q975", "travel_hrs","protected","adm2", "confilct", "district", "disaster", "dis_num", "dis_deaths", "NDVI", "LST", "elevation", 
                      "st_cluster", "dist_cluster", "cattle_PO_mean", "cattle_PO_sum", "pigs_PO_mean", "pigs_PO_sum", "prox1km")

'%ni%' <- Negate('%in%')

HAT@data$st.cluster<-1:nrow(HAT@data)

HAT@data$disaster<-ifelse(HAT@data$disaster=="None", 0,1)

#Set up predictors
Y<-c("wealth", "NDVI", "LST", "cases", "livestock")
if(type=="rhod"){
  X<-list(c("livestock_lag1", "conflict", "disaster", "wealth_lag"),
          c("livestock_lag1", "protected", "elevation", "conflict", "disaster"),
          c("livestock_lag1", "protected", "elevation", "conflict", "disaster"),
          c("NDVI_lag", "LST_lag", "conflict_lag", "disaster_lag", "wealth_lag", "elevation", "protected", "PO_lag2", "livestock_lag2", "offset(log(WorldPop))"),
          c("wealth", "protected", "conflict", "disaster", "elevation", "LST", "NDVI", "livestock_lag1"))
}else{
  X<-list(c("livestock_lag1", "conflict", "disaster", "wealth_lag"),
          c("livestock_lag1", "elevation", "conflict", "disaster"),
          c("livestock_lag1", "elevation", "conflict", "disaster"),
          c("NDVI_lag", "LST_lag", "conflict_lag", "disaster_lag", "wealth_lag", "elevation", "PO_lag2", "livestock_lag2", "offset(log(WorldPop))"),
          c("wealth", "conflict", "disaster", "elevation", "LST", "NDVI", "livestock_lag1"))
}
family<-c("gaussian", "gaussian", "gaussian", "poisson", "gaussian")

lag1vars<-c("NDVI", "LST", "wealth", "livestock", "conflict", "disaster")
lag1names<-paste0(lag1vars, "_lag")
lag1names[which(lag1names=="livestock_lag")]<-"livestock_lag1"

lag2vars<-c("livestock", "PO")
lag2names<-paste0(lag2vars, "_lag2")

cat.vars<-c("conflict_lag", "conflict", "disaster", "disaster_lag", "protected")

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

#Set exposure levels to compare
tx.vec<-c("a0", "a1")
nSamp=100

res.mean=data.frame(matrix(nrow=length(years), ncol=3))
names(res.mean)<-c("year", tx.vec)
res.mean$year<-years

res.total=data.frame(matrix(nrow=length(years), ncol=3))
names(res.total)<-c("year", tx.vec)
res.total$year<-years

res.inc=data.frame(matrix(nrow=length(years), ncol=3))
names(res.inc)<-c("year", tx.vec)
res.inc$year<-years

results.mean=list()
results.total=list()
results.inc=list()

validation=list()

for(i in 1:nSamp){
  res<-gform(tx.vec, HAT, years, type, species,
                Y.vars=Y, X.vars=X, family, test, lag1vars, 
                lag2vars, lag1names, 
                lag2names, cat.vars)
  
  if(test=="yes"){
    results.mean=c(results.mean, list(res$mean))
    results.total=c(results.total, list(res$total))
    results.inc=c(results.inc, list(res$inc))
  }else{
    validation=c(validation, list(res))
  }
  print(i)
}

if(test=="no"){
  for(i in 1:length(validation)){
    obs.outcome<-validation[[i]]$observed.outcome
    pred.outcome<-validation[[i]]$predicted.outcome
    perc.off<-(abs(obs.outcome-pred.outcome)/obs.outcome)*100
    validation[[i]]$perc.off<-perc.off
  }
  
  badly.off<-as.data.frame(matrix(NA, nrow=length(validation), ncol=nrow(validation[[1]])))
  names(badly.off)<-rownames(validation[[1]])
  
  for(i in 1:length(validation)){
    perc.off<-validation[[i]]$perc.off
    idx<-ifelse(perc.off>100,1,0)
    badly.off[i,]<-idx
  }
  
  for(i in 1:ncol(badly.off)){
    print(paste0(colnames(badly.off)[i], ": ", mean(badly.off[,i], na.rm=T)))
  }
}

if(test=="no"){
  if(species=="cattle"){
    save(file=paste0("gform_results/DRC/natural_course_cattle_ug_",type, ".Rdata"), validation)
  }else{
    save(file=paste0("gform_results/DRC/natural_course_pigs_ug_",type, ".Rdata"), validation)
  }
}

if(test=="yes"){
  if(species=="cattle"){
    save(file=paste0("gform_results/DRC/gform_cattle_mean_",type, ".Rdata"), results.mean)
    save(file=paste0("gform_results/DRC/gform_cattle_total_",type, ".Rdata"), results.total)
    save(file=paste0("gform_results/DRC/gform_cattle_inc_",type, ".Rdata"), results.inc)
  }else{
    save(file=paste0("gform_results/DRC/gform_pigs_mean_",type, ".Rdata"), results.mean)
    save(file=paste0("gform_results/DRC/gform_pigs_total_",type, ".Rdata"), results.total)
    save(file=paste0("gform_results/DRC/gform_pigs_inc_",type, ".Rdata"), results.inc)
  }
}