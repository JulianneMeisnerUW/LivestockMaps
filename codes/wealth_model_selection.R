library(rgdal)
library(dplyr)
library(ggplot2)
library(raster)
library(INLA)
`%ni%` = Negate(`%in%`)
source("codes/functions_prediction.R")

country<-"Malawi"

folderData<-paste0("plots/", country)

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

if(country=="Uganda"){
  country.short<-"ug"
  n<-11
  wealth_data<-readOGR("Data/Confounder_data/Uganda/wealth", layer="wealth")
  height=4
  width=10
}

if(country=="Malawi"){
  n<-15
  country.short<-"mw"
  wealth_data<-readOGR("Data/Confounder_data/Malawi/wealth", layer="wealth")
  height=4
  width=10
}

if(country=="DRC"){
  n<-3
  country.short<-"drc"
  wealth_data<-readOGR("Data/Confounder_data/DRC/wealth", layer="wealth")
  height=4
  width=10
}

########
# LOOCV
########
MSE_point<-data.frame(matrix(nrow=1, ncol=4))
colnames(MSE_point)<-c("Survey", "Model", "MSE", "Rank")
for(i in 1:n){
  load(paste0("results/wealth/", country, "/validation_wealth", i,"_", country, ".Rdata"))
  print(paste(Validation[[1]]$survey, ":", length(Validation)))
  newdat<-data.frame(matrix(nrow=length(Validation), ncol=4))
  colnames(newdat)<-c("Survey", "Model", "MSE", "Rank")
  for(j in 1:length(Validation)){
    newdat$MSE[j]<-Validation[[j]]$MSE_point
    newdat$Survey[j]<-as.character(Validation[[j]]$survey)
    newdat$Model[j]<-substring(Validation[[j]]$model,7)
  }
  newdat$Model<-as.numeric(newdat$Model)
  ord<-order(newdat$MSE)[newdat$Model]
  newdat$Rank<-ord
  MSE_point<-rbind(MSE_point, newdat)
  if(i==1){
    MSE_point<-MSE_point[-1,]
  }
}

MSE_point$Model<-as.factor(MSE_point$Model)
MSE_point<-MSE_point[which(MSE_point$Model!=5),] 
MSE_subset<-MSE_point[which(MSE_point$MSE<1),]
ndeleted<-nrow(MSE_point)-nrow(MSE_subset) 
deleted<-MSE_point[(rownames(MSE_point)%ni% rownames(MSE_subset)),] 

MSE_subset$Survey<-as.factor(MSE_subset$Survey)
if(country=="Malawi"){
  levels(MSE_subset$Survey)<-c("DHS2000", "DHS2004", "DHS2010", "DHS2016", "IHPS2013", "IHPS2016", "IHPS2004", "IHS2010", "MICS2013", "MIS2012", "TARIS2006", "WHS2003")
}
if(country=="Uganda"){
  levels(MSE_subset$Survey)<-c("AIS2011", "DHS2001", "DHS2006", "DHS2011", "DHS2016", "MIS2009", "MIS2014", "MIS2018", "NPS2009", "NPS2010", "NPS2011")
}

png(paste0("plots/", country, "/LOOCV_wealth.png"), width=width, height=height, units="in", res=300)
ggplot(MSE_subset, aes(x=Survey, y=MSE)) + geom_point(aes(color=Model)) +geom_text(data=subset(MSE_subset, MSE > 0.8),aes(Survey,MSE,label=Model), position=position_dodge2(width=0.5))+theme_bw()
dev.off()

MSE_summary<-data.frame(matrix(nrow=4, ncol=5))
colnames(MSE_summary)<-c("Model","MSE", "Mean_rank", "Min_rank", "Max_rank")

for(i in 1:4){
  subset<-subset(MSE_point, Model==i)
  MSE_summary$MSE[i]<-median(subset$MSE, na.rm=T)
  MSE_summary$Model[i]<-i
  MSE_summary$Mean_rank[i]<-mean(subset$Rank, na.rm=T)
  MSE_summary$Min_rank[i]<-min(subset$Rank, na.rm=T)
  MSE_summary$Max_rank[i]<-max(subset$Rank, na.rm=T)
}

##################################################
#Compare to WorldPop under $2/day (for 2010-2011)
##################################################
#Uganda and Malawi only, not produced for DRC

if(country=="Malawi"){
  districts<-readOGR(paste("Data/Exposure_data/Malawi/shapefiles/mwi_admbnda_adm2_nso_20181016", sep=""),layer="mwi_admbnda_adm2_nso_20181016")
  formedMap<-districts
  gridP<-readOGR(paste("Data/Exposure_data/Malawi/Created_datasets/prediction_data", sep=""), layer = "prediction_data")
  poverty2<-raster("Data/Confounder_data/Malawi/mwi11povcons200.tif")
  HAT<-readOGR("results/regression/Malawi/HATData", layer="MW_regress", p4s=proj)
  colnames(HAT@data)<-c("District", "year", "focus", "cases", "density_med_c", 
                        "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                        "density_sd_p", "density_q025_p", "density_q975_p", "counts_med_c", "counts_sd_c", "counts_q025_c", "counts_q975_c", "counts_med_p", 
                        "counts_sd_p", "counts_q025_p", "counts_q975_p","source", 
                        "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                        "wealth_med", "wealth_sd", "wealth_q025", "wealth_975")
}

if(country=="Uganda"){
  districts<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_4")
  gridP<-readOGR(paste("Data/Exposure_data/Uganda/Created_datasets/prediction_data", sep=""), layer = "prediction_data")
  poverty2<-raster("Data/Confounder_data/Uganda/uga10povcons200.tif")
  if(country=="Uganda"){
    formedMap<-spTransform(districts, CRS("+proj=utm +zone=35N +ellps=WGS84"))
    formedMap2<-rgeos::gBuffer(formedMap, byid=TRUE, width=0)
    formedMap<-spTransform(formedMap2, CRS("+proj=longlat +datum=WGS84"))
  }
    HAT<-readOGR("results/regression/Uganda/HATData", layer="ug_regress_gamb", p4s=proj)
    colnames(HAT@data)<-c("year", "focus", "cases", "density_med_c", 
                            "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                            "density_sd_p", "density_q025_p", "density_q975_p","source", "prox_1km",
                            "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                            "wealth_med", "wealth_sd", "wealth_q025", "wealth_975")
}


spUnion = maptools::unionSpatialPolygons(formedMap, rep(1, dim(formedMap)[1]))
loc.mesh = matrix(0, nrow = 0, ncol = 2)
for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
  loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
}

HAT10<-HAT[which(HAT@data$year==2010),]
merged.trans <- spTransform(HAT10, crs(poverty2))
pov.extract<-raster::extract(poverty2, merged.trans)
HAT10@data$prop_u2<-pov.extract

mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                    max.edge = c(0.3,0.6), 
                    offset = c(1, 1))

spde_RE<-inla.spde2.pcmatern(mesh, prior.range=c(0.3,0.05), prior.sigma=c(1,0.05))

pc.u=1
pc.alpha=0.01

hyperpc1<-list(prec=list(prior="pc.prec", param=c(pc.u, pc.alpha))) 

iset <- inla.spde.make.index('i', n.spde=spde_RE$n.spde)

A<-inla.spde.make.A(mesh, loc=HAT10@coords)

stk.dat<-inla.stack(
  data=list(resp=HAT10@data$wealth_med), 
  A=list(A,1), 
  effects=list(iset, HAT10@data), 
  tag='est')

formula <- resp ~ 0 + intercept + prop_u2 + f(i, model=spde_RE) 
res<-inla (formula, data=inla.stack.data(stk.dat), family="gaussian", control.predictor=list(A=inla.stack.A(stk.dat)), 
           control.compute=list(config = TRUE))
