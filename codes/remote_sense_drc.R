rm(list = ls())
library(stringr)
library(gdalUtils)
library(rgdal)
library(raster)
library(FedData)

country<-"DRC"

districts<-readOGR("Data/Exposure_data/DRC/shapefiles/gadm36_COD_shp", layer="gadm36_COD_2") 
geographic<-proj4string(districts)
polygon<- polygon_from_extent(districts,proj4string=geographic)

HAT<-readOGR(paste0("Data/Outcome_data/", country, "/HATData_dis"), layer="HATData_dis_gamb")
colnames(HAT@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                        "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                        "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                        "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                        "wealth_975", "travel_hrs","protected5km", "adm2", "conflict", "district", "disaster", "dis_num", 
                        "dis_deaths", "interaction")

years<-c(2008:2015)

med.names<-c(paste("NDVI", years, sep="_"), paste("LST", years, sep="_"))
mediators<-data.frame(matrix(nrow=nrow(HAT@data), ncol=length(med.names)))
colnames(mediators)<-med.names                            

#--------------------------------#
# NDVI (LTDR5; an AVHRR product)
#--------------------------------#
# Get a list of sds names
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

files=c("h19v09", "h19v08", "h19v10","h20v08","h20v09", "h20v10","h21v08")
extracts<-matrix(NA, ncol=length(files), nrow=nrow(HAT@data))

for(year in years){
  yIdx=which(years==year)
  print(year)
  raster_list=list()
  for(fIdx in 1:length(files)){
    filenames<-list.files('Data/Mediator/MODIS/')
    filenames.short<-str_sub(filenames, start=10, end= 13)
    search<-which(filenames.short==year)
    filenames.year<-filenames[search]
    use<-filenames.year[which(grepl(files[fIdx], filenames.year)==TRUE)][1]
    sds <- get_subdatasets(paste0('Data/Mediator/MODIS/', use))
    name <- sds[1]
    filename <- paste0('Data/Mediator/MOD21_',year,'_', fIdx, '.tif')
    gdal_translate(sds[1], dst_dataset = filename)
    LST <- raster(filename)
    raster_list=c(raster_list, LST)
  }
  
  master_raster<-do.call(merge, raster_list)
  
  HAT.trans <- spTransform(HAT, crs(master_raster))
  LST<-raster::extract(master_raster, HAT.trans)

  mediators[,(yIdx+length(years))]<-LST
}

HAT@data$NDVI<-rep(NA, nrow(HAT@data))
HAT@data$LST<-rep(NA, nrow(HAT@data))

#Merge with HAT data
for(i in 1:length(unique(HAT@data$year))){
  year=unique(HAT@data$year)[i]
  NDVI<-mediators[,i]
  LST<-mediators[,i+length(years)]
  HAT@data$NDVI[which(HAT@data$year==year)]<-NDVI
  HAT@data$LST[which(HAT@data$year==year)]<-LST
}

HAText<-extent(HAT)

#------------#
#Elevation
#------------#
gmted2 <- raster("Data/Predictor_data/Elevation/10s000e_20101117_gmted_med075.tif") 
gmted3 <- raster("Data/Predictor_data/Elevation/30s030e_20101117_gmted_med075.tif") 
gmted4 <- raster("Data/Predictor_data/Elevation/30s000e_20101117_gmted_med075.tif") 
gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
gmted2.crop<-crop(gmted2, HAText) 
gmted3.crop<-crop(gmted3, HAText)
gmted4.crop<-crop(gmted4, HAText) 
gmted6.crop<-crop(gmted6, HAText)
merged.trans <- spTransform(HAT, crs(gmted2.crop))
elevation.1<-raster::extract(gmted2.crop, merged.trans)
elevation.2<-raster::extract(gmted3.crop, merged.trans)
elevation.4<-raster::extract(gmted4.crop, merged.trans)
elevation.3<-raster::extract(gmted6.crop, merged.trans)
elev1<-ifelse(is.na(elevation.1)==TRUE, elevation.2, elevation.1)
elev2<-ifelse(is.na(elev1)==TRUE, elevation.3, elev1)
elev3<-ifelse(is.na(elev2)==TRUE, elevation.4, elev2)
HAT@data$elevation<-elev3


writeOGR(obj=HAT, dsn="Data/Outcome_data/DRC/HAT_final", layer="HAT_final_gamb", driver="ESRI Shapefile", overwrite_layer=TRUE)
