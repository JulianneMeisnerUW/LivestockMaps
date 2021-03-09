#------------Code adapted from Yunhan Wu----------------#

rm(list = ls())
library(rgdal)
library(spdep)
library(geosphere)
library(raster)
library(mapplots)

country<-"Uganda"

folderData<-paste0("Data/Exposure_data/", country, "/Created_datasets/Urbanicity/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

if(country=="Malawi"){
  s.frames<-c(1998, 2008)
  country.short<-"mwi"
  country.gadm<-"MWI"
  humdata.code<-"nso_20181016"
  UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
}

if(country=="Uganda"){
  s.frames<-c(1991, 2002, 2014)
  country.short<-"uga"
  country.gadm<-"UGA"
  humdata.code<-"ubos_20200824"
  UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
}

if(country=="DRC"){
  s.frames<-c(2003, 2010)
  country.short<-"cod"
  country.gadm<-"COD"
  humdata.code<-"rgc_20190911"
  UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
}

#Load survey data
data<-readOGR(paste0("Data/Confounder_data/", country, "/wealth"), layer="wealth") 

data@data$x<-data@coords[,1]
data@data$y<-data@coords[,2]

data@data$pop_den<-data@data$noaa<-data@data$prot_dist<-rep(NA, nrow(data@data))
data@data$cluster<-1:nrow(data@data)

#Extract WorldPop population density and NOAA nighttime lights
for(i in 1:length(s.frames)){
  if(i<length(s.frames)){
    subset<-data[which(data@data$year>=s.frames[i]&data@data$year<s.frames[i+1]),]
  }else{
    subset<-data[which(data@data$year>=s.frames[i]),]
  }
  if(s.frames[i]<2000){
    wp.year<-2000
  }else{
    wp.year<-s.frames[i]
  }
  
  nl.year<-s.frames[i]
  
  if(s.frames[i]<1992){
    nl.year<-1992
  }
  if(s.frames[i]>2013){
    nl.year<-2013
  }
  pop<-raster(paste0("Data/Predictor_data/Urbanicity/", country, "/", country.short, "_ppp_", wp.year, "_1km_Aggregated_UNadj.tif"))
  merged.trans <- spTransform(subset, crs(pop))  
  pop.extract<-raster::extract(pop, merged.trans)
  ids<-which(data@data$cluster%in%subset@data$cluster)
  data@data$pop_den[ids]<-pop.extract
  
  nightlights<-raster(paste('Data/Predictor_data/Nighttime_lights/nightlights_',nl.year,'.tif', sep=""))
  
  merged.trans <- spTransform(subset, crs(nightlights))  
  noaa.extract<-raster::extract(nightlights, merged.trans)
  data@data$noaa[ids]<-noaa.extract
  data@data$noaa<-ifelse(data@data$noaa==255, NA, data@data$noaa)
  
}

data@data$log_pop<-log(data@data$pop_den)
data@data$log_noaa<-log(data@data$noaa)
data@data$sqrt_pop<-sqrt(data@data$pop_den)
data@data$sqrt_noaa<-sqrt(data@data$noaa)

#Extract administrative levels
shapefile<-readOGR(paste0("Data/Predictor_data/Urbanicity/", country, "/shapeFiles_humdata"), layer=paste0(country.short, "_admbnda_adm2_", humdata.code)) 
admin<-over(data, shapefile)
if(country=="Uganda"|country=="Malawi"){
  data@data$adm0.5<-admin$ADM1_EN
  data@data$adm1<-admin$ADM2_EN
}else{
  data@data$adm0.5<-admin$ADM1_FR
  data@data$adm1<-admin$ADM2_FR
}

saveRDS(data@data, file = paste0(folderData, country.short, "_crc.rds"))

#Create 1km grid and do the same thing
for(i in 1:length(s.frames)){
  
  if(s.frames[i]<2000){
    wp.year<-2000
  }else{
    wp.year<-s.frames[i]
  }
  
  nl.year<-s.frames[i]
  
  if(s.frames[i]<1992){
    nl.year<-1992
  }
  if(s.frames[i]>2013){
    nl.year<-2013
  }
  
  pop<-raster(paste0("Data/Predictor_data/Urbanicity/", country, "/", country.short, "_ppp_", wp.year, "_1km_Aggregated_UNadj.tif"))
  nightlights<-raster(paste('Data/Predictor_data/Nighttime_lights/nightlights_',nl.year,'.tif', sep=""))
  
  urb_dat<-as.data.frame(coordinates(pop))
  colnames(urb_dat)<-c('x','y')
  urb_dat$pop_den<-extract(pop, urb_dat[c('x','y')])
  urb_dat$noaa<-extract(nightlights,urb_dat[c('x','y')])
  urb_dat$noaa<-ifelse(urb_dat$noaa==255, NA, urb_dat$noaa)
  if(country=="Uganda"){
    urb_dat$prot_dist<-extract(prot_dist,urb_dat[c('x','y')])
  }
  
  grid <- SpatialPointsDataFrame(urb_dat[c('x','y')], data=urb_dat[c('pop_den', 'noaa')], proj4string = crs(shapefile))
  
  admin<-over(grid, shapefile)
  if(country=="Uganda"|country=="Malawi"){
    urb_dat$adm0.5<-admin$ADM1_EN
    urb_dat$adm1<-admin$ADM2_EN
  }else{
    urb_dat$adm0.5<-admin$ADM1_FR
    urb_dat$adm1<-admin$ADM2_FR
  }
  
  urb_dat$complete<-complete.cases(urb_dat)

  saveRDS(urb_dat, file = paste0(folderData, "/nat_grid_", s.frames[i],".rds"))
}