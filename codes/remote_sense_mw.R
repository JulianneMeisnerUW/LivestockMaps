library(gdalUtils)
library(rgdal)
library(raster)
library(FedData)

eas<-readOGR("Data/Exposure_data/Malawi/shapefiles/eas_bnd", layer="eas_bnd")
geographic<-proj4string(eas)
polygon<- polygon_from_extent(eas,proj4string=geographic)

HAT<-readOGR("Outcome_data/Malawi/HATData_dis", layer="HATData_dis", p4s=proj)
colnames(HAT@data)<-c("District","year","focus","cases","density_med_c", "density_sd_c",  
                      "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                      "density_q025_p", "density_q975_p", "counts_med_c", "counts_sd_c",
                      "counts_q025_c", "counts_q975_c", "counts_med_p", "counts_sd_p", 
                      "counts_q025_p", "counts_q975_p", "source", "cluster", "intercept",
                      "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                      "wealth_975","prox1km", "travel_hrs","protected5km", "disaster", "dis_num", "dis_deaths")

years<-c(2000:2018)

med.names<-c(paste("NDVI", 2000:2018, sep="_"), paste("LST", 2000:2018, sep="_"))
mediators<-data.frame(matrix(nrow=nrow(HAT@data), ncol=length(med.names)))
colnames(mediators)<-med.names                            

#--------------------------------#
# NDVI (LTDR5; an AVHRR product)
#--------------------------------#
for(year in years){
  yIdx=which(years==year)
  sds <- get_subdatasets(paste0('Data/Mediator/AVHRR/AVHRR_',year,'.hdf'))
  name <- sds[1]
  filename <- paste0('Data/Mediator/LTDR5_',year,'.tif')
  gdal_translate(sds[1], dst_dataset = filename)
  NDVI <- raster(filename)
  NDVI.crop<-crop(NDVI, polygon)
  HAT.trans <- spTransform(HAT, crs(NDVI.crop))
  HAT.extract<-raster::extract(NDVI.crop, HAT.trans)
  mediators[,yIdx]<-HAT.extract
}

#-------------------------#
# LST (MODIS MOD21)
#-------------------------#
proj<-"+proj=longlat +ellps=clrk66 +no_defs "

for(year in years){
  yIdx=which(years==year)
  sds1 <- get_subdatasets(paste0('Data/Mediator/MODIS/MODIS_',year,'_1.hdf'))
  name <- sds1[1]
  filename1 <- paste0('Data/Mediator/MOD21_',year,'_1.tif')
  gdal_translate(sds1[1], dst_dataset = filename1)
  LST1 <- raster(filename1)
  LST_new1<-projectRaster(LST1, crs=proj)
  LST.crop1<-crop(LST_new1, polygon)
  HAT.trans1 <- spTransform(HAT, crs(LST.crop1))
  HAT.extract1<-raster::extract(LST.crop1, HAT.trans1)
  
  sds2 <- get_subdatasets(paste0('Data/Mediator/MODIS/MODIS_',year,'_2.hdf'))
  name <- sds2[1]
  filename2 <- paste0('Data/Mediator/MOD21_',year,'_2.tif')
  gdal_translate(sds2[1], dst_dataset = filename2)
  LST2 <- raster(filename2)
  LST_new2<-projectRaster(LST2, crs=proj)
  LST.crop2<-crop(LST_new2, polygon)
  HAT.trans2 <- spTransform(HAT, crs(LST.crop2))
  HAT.extract2<-raster::extract(LST.crop2, HAT.trans2)
  
  temp<-ifelse(is.na(HAT.extract1)==TRUE, HAT.extract2, HAT.extract1)
  mediators[,yIdx+17]<-temp
}

NDVI_2000 <- raster('Data/Mediator/LTDR5_2000.tif')
plot(NDVI_2000)

LST_2000<-raster('Data/Mediator/MOD21_2000_1.tif')
plot(LST_2000)

HAT@data$NDVI<-rep(NA, nrow(HAT@data))
HAT@data$LST<-rep(NA, nrow(HAT@data))

#Merge with HAT data
for(i in 1:nrow(HAT@data)){
  year=HAT@data$year[i]
  yIdx=which(years==year)
  NDVI<-mediators[,yIdx]
  LST<-mediators[,yIdx+17]
  HAT@data$NDVI[i]<-NDVI[i]
  HAT@data$LST[i]<-LST[i]
}

HAText<-extent(HAT)

#------------#
#Elevation
#------------#
gmted3 <- raster("Data/Predictor_data/Elevation/30s030e_20101117_gmted_med075.tif") 
gmted3.crop<-crop(gmted3, HAText)
merged.trans <- spTransform(HAT, crs(gmted3.crop))
elevation<-raster::extract(gmted3.crop, merged.trans)
HAT@data$elevation<-elevation

writeOGR(obj=HAT, dsn="Data/Outcome_data/Malawi/HAT_final", layer="HAT_final", driver="ESRI Shapefile", overwrite_layer=TRUE)
