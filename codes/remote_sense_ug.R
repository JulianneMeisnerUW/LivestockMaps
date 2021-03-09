library(stringr)
library(gdalUtils)
library(rgdal)
library(raster)
library(FedData)

districts<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_4")
geographic<-proj4string(districts)
polygon<- polygon_from_extent(districts,proj4string=geographic)

HAT.g<-readOGR("Data/Outcome_data/Uganda/HATData_dis", layer="HATData_dis_gamb", p4s=geographic)
HAT.r<-readOGR("Data/Outcome_data/Uganda/HATData_dis", layer="HATData_dis_rhod", p4s=geographic)

colnames(HAT.g@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                      "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                      "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                      "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                      "wealth_975","district", "travel_hrs","protected","parish", "conflict", 
                      "disaster", "dis_num", "dis_deaths")

colnames(HAT.r@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                        "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                        "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                        "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                        "wealth_975","district", "travel_hrs","protected","parish", "conflict",
                        "disaster", "dis_num", "dis_deaths")

years<-c(2000:2018)

med.names<-c(paste("NDVI", 2000:2018, sep="_"), paste("LST", 2000:2018, sep="_"))
mediators.g<-data.frame(matrix(nrow=nrow(HAT.g@data), ncol=length(med.names)))
colnames(mediators.g)<-med.names                            

mediators.r<-data.frame(matrix(nrow=nrow(HAT.r@data), ncol=length(med.names)))
colnames(mediators.r)<-med.names                            

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
  HAT.g.trans <- spTransform(HAT.g, crs(NDVI.crop))
  HAT.g.extract<-raster::extract(NDVI.crop, HAT.g.trans)
  mediators.g[,yIdx]<-HAT.g.extract
  
  HAT.r.trans <- spTransform(HAT.r, crs(NDVI.crop))
  HAT.r.extract<-raster::extract(NDVI.crop, HAT.r.trans)
  mediators.r[,yIdx]<-HAT.r.extract
}

#-------------------------#
# LST (MODIS MOD21)
#-------------------------#
proj<-"+proj=longlat +ellps=clrk66 +no_defs "

files=c("h20v08", "h20v09", "h21v08", "h21v09")
extracts.g<-matrix(NA, ncol=length(files), nrow=nrow(HAT.g@data))
extracts.r<-matrix(NA, ncol=length(files), nrow=nrow(HAT.r@data))

for(year in years){
  yIdx=which(years==year)
  print(year)
  for(fIdx in 1:length(files)){
    filenames<-list.files('Data/Mediator/MODIS/')
    filenames.short<-str_sub(filenames, start=10, end= 13)
    search<-which(filenames.short==year)
    filenames.year<-filenames[search]
    use<-filenames.year[which(grepl(files[fIdx], filenames.year)==TRUE)][1]
    sds <- get_subdatasets(paste0('Data/Mediator/MODIS/original/', use))
    name <- sds[1]
    filename <- paste0('Data/Mediator/MOD21_',year,'_', fIdx, '.tif')
    gdal_translate(sds[1], dst_dataset = filename)
    LST <- raster(filename)
    LST_new<-projectRaster(LST, crs=crs(HAT))
    LST.crop<-try(crop(LST_new, polygon), TRUE)
    if(class(crs(LST.crop))=="CRS"){
      HAT.g.trans <- spTransform(HAT.g, crs(LST.crop))
      HAT.g.extract<-raster::extract(LST.crop, HAT.g.trans)
      extracts.g[,fIdx]<-HAT.g.extract
      
      HAT.r.trans <- spTransform(HAT.r, crs(LST.crop))
      HAT.r.extract<-raster::extract(LST.crop, HAT.r.trans)
      extracts.r[,fIdx]<-HAT.r.extract
    }else{
      print(paste0(year, "_", fIdx))
      extracts.g[,fIdx]<-NA
      extracts.r[,fIdx]<-NA
    }
  }
  results.g<-c()
  results.r<-c()
  for(i in 1:nrow(extracts.g)){
    rows<-extracts.g[i,]
    res<-which(!is.na(rows))[1]
    lst<-extracts.g[i,res]
    results.g<-c(results.g, lst)
  }
  for(i in 1:nrow(extracts.r)){
    rows<-extracts.r[i,]
    res<-which(!is.na(rows))[1]
    lst<-extracts.r[i,res]
    results.r<-c(results.r, lst)
  }
  mediators.g[,yIdx+19]<-results.g
  mediators.r[,yIdx+19]<-results.r
}

HAT.g@data$NDVI<-rep(NA, nrow(HAT.g@data))
HAT.g@data$LST<-rep(NA, nrow(HAT.g@data))

HAT.r@data$NDVI<-rep(NA, nrow(HAT.r@data))
HAT.r@data$LST<-rep(NA, nrow(HAT.r@data))

#Merge with HAT data
for(i in 1:nrow(HAT.g@data)){
  year=HAT.g@data$year[i]
  yIdx=which(years==year)
  NDVI<-mediators.g[,yIdx]
  LST<-mediators.g[,yIdx+17]
  HAT.g@data$NDVI[i]<-NDVI[i]
  HAT.g@data$LST[i]<-LST[i]
}

#Merge with HAT data
for(i in 1:nrow(HAT.r@data)){
  year=HAT.r@data$year[i]
  yIdx=which(years==year)
  NDVI<-mediators.r[,yIdx]
  LST<-mediators.r[,yIdx+17]
  HAT.r@data$NDVI[i]<-NDVI[i]
  HAT.r@data$LST[i]<-LST[i]
}
HAText.g<-extent(HAT.g)
HAText.r<-extent(HAT.r)

#------------#
#Elevation
#------------#
gmted2 <- raster("Data/Predictor_data/Elevation/10s000e_20101117_gmted_med075.tif") 
gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
gmted6.crop<-crop(gmted6, HAText.g)
merged.trans <- spTransform(HAT.g, crs(gmted6.crop))
elevation<-raster::extract(gmted6.crop, merged.trans)
HAT.g@data$elevation<-elevation

gmted6.crop<-crop(gmted6, HAText.r)
merged.trans <- spTransform(HAT.r, crs(gmted6.crop))
elevation<-raster::extract(gmted6.crop, merged.trans)
HAT.r@data$elevation<-elevation

writeOGR(obj=HAT.g, dsn="Data/Outcome_data/Uganda/HAT_final", layer="HAT_final_gamb", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(obj=HAT.r, dsn="Data/Outcome_data/Uganda/HAT_final", layer="HAT_final_rhod", driver="ESRI Shapefile", overwrite_layer=TRUE)
