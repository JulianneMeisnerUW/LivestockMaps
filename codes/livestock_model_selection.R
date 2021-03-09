rm(list = ls())
library(rgdal)
library(dplyr)
library(ggplot2)
library(INLA)
library(raster)
`%ni%` = Negate(`%in%`)

country<-"DRC"

if(country=="Uganda"){
  country.short<-"ug"
  n<-10
}

if(country=="Malawi"){
  n<-12
  country.short<-"mw"
  type<-"rhod"
}

if(country=="DRC"){
  n<-2
  country.short<-"drc"
  type<-"gamb"
}

folderData<-paste0("plots/", country, "/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

########
# LOOCV
########

MSE_point<-data.frame(matrix(nrow=1, ncol=3))
colnames(MSE_point)<-c("Survey", "Model","MSE")

for(i in 1:n){
  load(paste0("results/livestock_mapping/", country, "/validation", i,"_", country, ".Rdata"))
  if(class(Validation)=="list"){
    print(paste("Survey ", i, ":", length(Validation)))
    newdat<-data.frame(matrix(nrow=length(Validation), ncol=3))
    colnames(newdat)<-c("Survey", "Model", "MSE")
    for(j in 1:length(Validation)){
      if(class(Validation[[j]])=="list"){
        newdat$MSE[j]<-Validation[[j]]$MSE_point
        newdat$Survey[j]<-as.character(Validation[[j]]$survey)
        newdat$Model[j]<-Validation[[j]]$model
      }
    }
    MSE_point<-rbind(MSE_point, newdat)
    if(i==1){
      MSE_point<-MSE_point[-1,]
    }
  }
}
MSE_point<-MSE_point[!is.na(MSE_point$Survey),]
MSE_point$Model<-substring(MSE_point$Model, 7)

MSE_point_subset<-MSE_point[which(MSE_point$MSE<1),]#Remove MSEs>1 for plotting purposes only
ndeleted<-nrow(MSE_point)-nrow(MSE_point_subset) 
deleted<-MSE_point[(rownames(MSE_point)%ni% rownames(MSE_point_subset)),] 

png(paste0("plots/", country, "/LOOCV_livestock_", country.short, ".png"), width=8, height=4, units="in", res=300)
ggplot(MSE_point_subset, aes(x=Survey, y=MSE)) + geom_point(aes(color=Model)) +geom_text(data=subset(MSE_point_subset, MSE > 0.1),aes(Survey,MSE,label=Model), nudge_y=-0.005) + theme_bw()+theme(text = element_text(size=15),
                                                                                                                                                                                                  axis.text.x = element_text(angle=45, hjust=1)) 
dev.off() #Anything>0.1 is labeled

MSE_point<- MSE_point %>% 
  mutate_all(.funs = function(x) replace(x, which(x == Inf | x == "N/A"), NA))

MSE_summary<-data.frame(matrix(nrow=6, ncol=2))
colnames(MSE_summary)<-c("Model", "MSE_point.LOOCV")
MSE_summary$Model<-c(1:6)

for(i in 1:6){
  subset<-subset(MSE_point, Model==i)
  MSE_summary$MSE_point.LOOCV[i]<-median(subset$MSE, na.rm=T)
}

#####################
# Make plots in time
#####################
if(country!="DRC"){
  years<-2000:2020
}else{
  years<-2008:2015
}

load(paste0("results/livestock_mapping/", country, "/prediction_cattle_", country.short, ".Rdata"))
  pdf(paste0("plots/", country, "/time_cattle.pdf"))
  for(m in c(1:6)){
    time_cattle<-as.data.frame(matrix(ncol=4, nrow=length(years)))
    colnames(time_cattle)<-c("Year", "Median", "Lower", "Upper")
    time_cattle$Year<-years
    for(year in years){
      idx=which(years%in%year)
      medIdx = (idx-1)*6 + 1
      upperIdx = (idx-1)*6 + 4
      lowerIdx = (idx-1)*6 + 3
      
      median<-Predict_nST.results.cattle[[m]]$prediction[[medIdx]]
      lower<-Predict_nST.results.cattle[[m]]$prediction[[lowerIdx]]
      upper<-Predict_nST.results.cattle[[m]]$prediction[[upperIdx]]
      
      time_cattle$Median[idx]<-mean(median, na.rm=T)
      time_cattle$Lower[idx]<-mean(lower, na.rm=T)
      time_cattle$Upper[idx]<-mean(upper, na.rm=T)
      print(year)
    }
    plot1<-ggplot() 
    plot2<- plot1 + geom_line(data = time_cattle, aes(x = Year, y = Median), color = "blue") +ggtitle(paste("Model", m)) 
    plot3<- plot2 + geom_line(data = time_cattle, aes(x = Year, y = Lower), color = "blue", linetype="dotted") 
    plot4<- plot3 + geom_line(data = time_cattle, aes(x = Year, y = Upper), color = "blue", linetype="dotted") + theme_bw()
    print(plot4)
  } 
dev.off()

##################################
# GLW correlation for final model
##################################
cattle<-raster("Data/External_validation/5_Ct_2010_Da.tif")#Dasymetric results, cattle
pigs<-raster("Data/External_validation/5_Pg_2010_Da.tif")#Dasymetric results, pigs

if(country=="Uganda"){
HAT<-readOGR("results/regression/Uganda/HATData", layer="ug_regress_gamb")
colnames(HAT@data)<-c("year", "focus", "cases", "density_med_c", 
                      "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                      "density_sd_p", "density_q025_p", "density_q975_p","source", "prox_1km",
                      "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                      "wealth_med", "wealth_sd", "wealth_q025", "wealth_975")
}else{
  HAT<-readOGR(paste0("results/regression/", country, "/HATData"), layer=paste0(country.short, "_regress_", type))
  colnames(HAT@data)<-c("year", "focus", "cases", "density_med_c", 
                        "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                        "density_sd_p", "density_q025_p", "density_q975_p","source", "prox_1km",
                        "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                        "wealth_med", "wealth_sd", "wealth_q025", "wealth_975")
}

HAT10<-HAT[which(HAT@data$year==2010),]
merged.trans <- spTransform(HAT10, crs(cattle))

cattle.extract<-raster::extract(cattle, merged.trans)
HAT10@data$GLW_c<-cattle.extract
merged.trans <- spTransform(HAT10, crs(pigs))
pig.extract<-raster::extract(pigs, merged.trans)
HAT10@data$GLW_p<-pig.extract

HAT10<-HAT10[which(HAT10@data$WorldPop>0),]

HAT10@data$GLW_dens_c<-HAT10@data$GLW_c/HAT10@data$WorldPop
HAT10@data$GLW_dens_p<-HAT10@data$GLW_p/HAT10@data$WorldPop

res.c<-lm(HAT10@data$GLW_dens_c~HAT10@data$density_med_c)
res.p<-lm(HAT10@data$GLW_dens_p~HAT10@data$density_med_p)
