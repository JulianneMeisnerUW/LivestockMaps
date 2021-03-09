library(gdalUtils)
library(rgdal)
library(raster)
library(FedData)
outcomeFolderData<-"Data/Outcome_data/South_Sudan/"

HAT<-readOGR(paste(outcomeFolderData, "SS_HAT", sep=""), layer="SS_HAT")
HAT.active<-readOGR(paste(outcomeFolderData, "SS_HAT", sep=""), layer="SS_HAT_active")

geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")

proj<-proj4string(HAT)

e<-extent(HAT)
polygon<-polygon_from_extent(HAT)

years<-c(2000:2016)

med.names<-c(paste("NDVI", 2000:2016, sep=""), paste("LST", 2000:2016, sep="_"))
mediators<-data.frame(matrix(nrow=nrow(HAT@data), ncol=length(med.names))) 
colnames(mediators)<-med.names   
rownames(mediators)<-HAT@data$new_adm2

mediators.active<-data.frame(matrix(nrow=nrow(HAT.active@data), ncol=length(med.names))) 
colnames(mediators.active)<-med.names   
rownames(mediators.active)<-HAT.active@data$new_adm2

#-------------------------------------------------------------------------#
# Extract NDVI (LTDR5; an AVHRR product) as mean over county for each year
#-------------------------------------------------------------------------#
# Get a list of sds names
for(year in years){
  yIdx=which(years==year)
  sds <- get_subdatasets(paste0('Data/Mediator/SS/AVHRR/AVHRR_',year,'.hdf'))
  name <- sds[1]
  filename <- paste0('Data/Mediator/SS/LTDR5_',year,'.tif')
  gdal_translate(sds[1], dst_dataset = filename)
  NDVI <- raster(filename)
  NDVI.crop<-crop(NDVI, polygon)
  HAT.trans <- spTransform(HAT, crs(NDVI.crop))
  HAT.active.trans <- spTransform(HAT.active, crs(NDVI.crop))
  v1 <- as.numeric(raster::extract(NDVI.crop, HAT.trans, fun=mean, na.rm=TRUE))
  mediators[,yIdx]<-v1
  v2 <- as.numeric(raster::extract(NDVI.crop, HAT.active.trans, fun=mean, na.rm=TRUE))
  mediators.active[,yIdx]<-v2
}

#-------------------------#
# LST (MODIS MOD21)
#-------------------------#
for(year in years){
  yIdx=which(years==year)
  raslist<-list()
  for(i in 3:4){ 
    sds <- get_subdatasets(paste0('Data/Mediator/SS/MODIS21/MODIS_',year,'_0', i, '.hdf'))
    name <- sds[1]
    filename <- paste0('Data/Mediator/SS/MOD21_',year,'_', i, '.tif')
    gdal_translate(sds[1], dst_dataset = filename)
    LST<-raster(filename)
    raslist[[i-2]]<-LST
  }
  LST<-merge(raslist[[1]], raslist[[2]])
  
  LST_new<-projectRaster(LST, crs=proj)
  LST.crop<-crop(LST_new, polygon)

  HAT.trans <- spTransform(HAT, crs(LST.crop))
  v1 <- as.numeric(raster::extract(LST_new, HAT.trans, fun=mean, na.rm=TRUE))
  mediators[,yIdx+17]<-v1
  
  HAT.active.trans <- spTransform(HAT.active, crs(LST.crop))
  v2 <- as.numeric(raster::extract(LST_new, HAT.active.trans, fun=mean, na.rm=TRUE))
  mediators.active[,yIdx+17]<-v2
}

#Merge with HAT data
HAT@data<-cbind(HAT@data, mediators)

HAT.active@data<-cbind(HAT.active@data, mediators.active)

writeOGR(obj=HAT, dsn="Data/Outcome_data/South_Sudan/HAT_RS_SS", layer="HAT_RS_SS", driver="ESRI Shapefile", overwrite_layer=TRUE)

writeOGR(obj=HAT.active, dsn="Data/Outcome_data/South_Sudan/HAT_RS_SS", layer="HAT_RS_SS_active", driver="ESRI Shapefile", overwrite_layer=TRUE)
