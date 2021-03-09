library(INLA)
library(rgdal)
library(sp)
library(splines)
library(gdalUtils)
library(raster)
library(FedData)
library(rgeos)
library(dplyr)
set.seed(1154)
country<-"DRC"
geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
if(country=="Malawi"){
  UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
  districts<-readOGR("Data/Exposure_data/Malawi/shapefiles/mwi_admbnda_adm2_nso_20181016", layer="mwi_admbnda_adm2_nso_20181016")
  geographic<-proj4string(districts)
  proj<-proj4string(districts)
  polygon<- polygon_from_extent(eas,proj4string=geographic)
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
  UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
  geographic<-proj4string(districts)
  proj<-proj4string(districts)
  polygon<- polygon_from_extent(districts,proj4string=geographic)
  gridP<-readOGR(paste("Data/Exposure_data/Uganda/Created_datasets/prediction_data", sep=""), layer = "prediction_data")
  poverty2<-raster("Data/Confounder_data/Uganda/uga10povcons200.tif")
  HAT.g<-readOGR(paste0("results/regression/Uganda/HATData"), layer="ug_regress_gamb", p4s=proj)
  colnames(HAT.g@data)<-c("year", "focus", "cases", "density_med_c", 
                          "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                          "density_sd_p", "density_q025_p", "density_q975_p","source", "prox_1km",
                          "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                          "wealth_med", "wealth_sd", "wealth_q025", "wealth_975")
  
  HAT.r<-readOGR(paste0("results/regression/Uganda/HATData"), layer="ug_regress_rhod", p4s=proj)
  colnames(HAT.r@data)<-c("year", "focus", "cases", "density_med_c", 
                          "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                          "density_sd_p", "density_q025_p", "density_q975_p","source", "prox_1km",
                          "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                          "wealth_med", "wealth_sd", "wealth_q025", "wealth_975")
}

if(country=="DRC"){
  districts<-readOGR("Data/Exposure_data/DRC/shapefiles/gadm36_COD_shp", layer="gadm36_COD_2") 
  UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
  geographic<-proj4string(districts)
  proj<-proj4string(districts)
  polygon<- polygon_from_extent(districts,proj4string=geographic)
  gridP<-readOGR(paste("Data/Exposure_data/DRC/Created_datasets/prediction_data", sep=""), layer = "prediction_data")
  
  HAT<-readOGR(paste0("results/regression/DRC/HATData"), layer=paste0(country.short, "_regress_", type))
  colnames(HAT@data)<-c("year", "focus", "cases", "density_med_c", 
                        "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                        "density_sd_p", "density_q025_p", "density_q975_p","source", "prox_1km",
                        "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                        "wealth_med", "wealth_sd", "wealth_q025", "wealth_975")
  
}

source('codes/my_functions.R')
inla.setOption("num.threads", 4)
'%ni%'<-Negate('%in%')

if(country=="Malawi"|country=="DRC"){
  HAT@data$cluster<-1:nrow(HAT@data)

  #Now remove clusters with 0 cases that are >5 hrs travel time from FHF capable of diagnosis of rHAT
  proj<-"+proj=longlat +ellps=clrk66 +no_defs "
  if(country=="Malawi"){
    fhf <- raster("Data/Denominator_data/WHO_travel_time_FHF/rHAT/Dx.tif")
  }else{
    fhf <- raster("Data/Denominator_data/WHO_travel_time_FHF/gHAT/Dx.tif")
  }
  fhf_new<-projectRaster(fhf, crs=proj)
  fhf.crop<-crop(fhf_new, polygon)
  HAT.trans <- spTransform(HAT, crs(fhf.crop))
  HAT.extract<-raster::extract(fhf.crop, HAT.trans)
  HAT.hours<-HAT.extract/60
  HAT@data$travel_hrs<-HAT.hours
  HAT<-HAT[which(HAT@data$travel_hrs<=5),]
  
  #Add in protected area
  if(country=="Malawi"){
    proj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    protected<-readOGR("Data/Exposure_data/Malawi/Created_datasets/prediction_data_protected", layer="prediction_data_protected", p4s=proj)
    merge<-over(HAT, protected)
    merge$protected<-ifelse(is.na(merge$WDPAID)==FALSE, 1,0)
    prot<-HAT[which(merge$protected==1),]
    nprot<-HAT[which(merge$protected==0),]
    
    #Project both in UTM: EPSG:20936
    prot.test<-spTransform(prot, CRS("+init=epsg:20936"))
    nprot.test<-spTransform(nprot, CRS("+init=epsg:20936"))
    points_matrix <- gWithinDistance(prot.test, nprot.test, dist = 5000, byid = TRUE)#This is in meters
    v <- rowSums(points_matrix, na.rm=TRUE) == 0
    nprot$protected5km<-rep(0, nrow(nprot@data))
    nprot$protected5km[v]<-1
    prot$protected5km<-1
    
    HAT<-bind(prot, nprot)
  }else{
    HAT@data$protected5km<-rep(NA, nrow(HAT@data))
  }
  
  writeOGR(obj=HAT, dsn=paste0("Data/Outcome_data/", country, "/HAT_regress"), layer="HAT_regress", driver="ESRI Shapefile", overwrite_layer=TRUE)
}

if(country=="Uganda"){
  HAT.g<-HAT.g[which(HAT.g@data$prox_1km==0),]
  HAT.r<-HAT.r[which(HAT.r@data$prox_1km==0),]
  
  #Restrict to admin 2s which reported cases:
  merge.g<-over(HAT.g, districts) 
  merge.r<-over(HAT.r, districts) 
  HAT.g$district<-merge.g$NAME_2
  HAT.r$district<-merge.r$NAME_2
  
  pos.sub.r<-HAT.r[which(HAT.r@data$cases>0),]
  pos.districts.r<-unique(pos.sub.r@data$district)
  HAT.r<-HAT.r[which(HAT.r$district%in%pos.districts.r),]
  
  pos.sub.g<-HAT.g[which(HAT.g@data$cases>0),]
  pos.districts.g<-unique(pos.sub.g@data$district)
  HAT.g<-HAT.g[which(HAT.g$district%in%pos.districts.g),]
    
  HAT.g@data$cluster<-1:nrow(HAT.g@data)
  HAT.r@data$cluster<-1:nrow(HAT.r@data)
  
  #Now remove clusters with 0 cases that are >5 hrs travel time from FHF capable of diagnosis of rHAT
  proj<-"+proj=longlat +ellps=clrk66 +no_defs "
  fhf.g <- raster("Data/Denominator_data/WHO_travel_time_FHF/gHAT/Dx.tif")
  fhf_new<-projectRaster(fhf.g, crs=proj)
  fhf.crop<-crop(fhf_new, polygon)
  HAT.g.trans <- spTransform(HAT.g, crs(fhf.crop))
  HAT.g.extract<-raster::extract(fhf.crop, HAT.g.trans)
  HAT.g.hours<-HAT.g.extract/60
  HAT.g@data$travel_hrs<-HAT.g.hours
  HAT.g<-HAT.g[which(HAT.g@data$travel_hrs<=5),]
  
  proj<-"+proj=longlat +ellps=clrk66 +no_defs "
  fhf.r <- raster("Data/Denominator_data/WHO_travel_time_FHF/rHAT/Dx.tif")
  fhf_new<-projectRaster(fhf.r, crs=proj)
  fhf.crop<-crop(fhf_new, polygon)
  HAT.r.trans <- spTransform(HAT.r, crs(fhf.crop))
  HAT.r.extract<-raster::extract(fhf.crop, HAT.r.trans)
  HAT.r.hours<-HAT.r.extract/60
  HAT.r@data$travel_hrs<-HAT.r.hours 
  HAT.r<-HAT.r[which(HAT.r@data$travel_hrs<=5),]
  
  #Add in protected area
  proj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  protected<-readOGR("Data/Exposure_data/Uganda/Created_datasets/prediction_data_protected", layer="prediction_data_protected", p4s=proj)
  merge.g<-over(HAT.g, protected)
  merge.r<-over(HAT.r, protected)
  
  HAT.g@data$protected<-ifelse(is.na(merge.g$WDPAID)==FALSE, 1,0)
  HAT.r@data$protected<-ifelse(is.na(merge.r$WDPAID)==FALSE, 1,0)
  
  writeOGR(obj=HAT.g, dsn=paste0("Data/Outcome_data/", country, "/HAT_regress"), layer="HAT_regress_gamb", driver="ESRI Shapefile", overwrite_layer=TRUE)
  writeOGR(obj=HAT.r, dsn=paste0("Data/Outcome_data/", country, "/HAT_regress"), layer="HAT_regress_rhod", driver="ESRI Shapefile", overwrite_layer=TRUE)

}

