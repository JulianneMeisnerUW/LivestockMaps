library(rgdal)
library(sp)
library(raster)
library(foreign)
library(SpatialEpi)
library(FedData)
library(survey)
library(rgeos)
library(RColorBrewer)
library(INLA)
library(uwIntroStats)
library(gpclib)
library(ggplot2)
library(splancs)
library(fields)
library(maptools)
library(parallel)
geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
inla.setOption("num.threads", 4)
inla.setOption("mkl",TRUE)
source('codes/my_functions.R')
'%ni%' <- Negate('%in%')
folderData<-"Data/Exposure_data/Malawi/"

#----------------------#
#   Read in data
#----------------------#
merged_data<-readOGR(paste(folderData, "Created_datasets/merged_data", sep=""), layer="merged_map") 
colnames(merged_data@data)<-c("members", "mem.pigs", "mem.cattle",
                              "urban", "n.cattle", "n.pigs", "survey", "year", "base_urb",
                              "district", "cluster", "water.body", "protected", "elevation")

proj.sp<-readOGR(paste(folderData, "Created_datasets/prediction_data", sep=""), layer="prediction_data")
colnames(proj.sp@data)<-c("protected", "water.body", "elevation")
proj<-proj4string(merged_data)
districts<-readOGR(paste(folderData, "/shapefiles/mwi_admbnda_adm2_nso_20181016", sep=""), layer="mwi_admbnda_adm2_nso_20181016", p4s=proj)
malawiMap<- districts

merged_data@data$urban<-as.numeric(as.character(merged_data@data$urban))

merged_data@data$mem.cattle<-as.numeric(as.character(merged_data@data$mem.cattle))
merged_data@data$mem.pigs<-as.numeric(as.character(merged_data@data$mem.pigs))
merged_data@data$members<-as.numeric(as.character(merged_data@data$members))

cattle_data<-merged_data[which(merged_data@data$mem.cattle>0),]
cattle_data@data$intercept<-rep(1, length(nrow(cattle_data@data)))
pig_data<-merged_data[which(merged_data@data$mem.pigs>0),]
pig_data@data$intercept<-rep(1, length(nrow(pig_data@data)))

all.years<-2000:2020

cattle_data@data$cluster = 1:nrow(cattle_data@data) 
rownames(cattle_data@data) = 1:nrow(cattle_data@data)
cattle_data@data<-cattle_data@data[,c("n.cattle", "mem.cattle", "intercept","urban", "protected","water.body", "elevation",
                                       "survey", "cluster", "district", "year")]
cattle_data@data$lat<-cattle_data@coords[2]
cattle_data@data$lon<-cattle_data@coords[1]

pig_data@data$cluster = 1:nrow(pig_data@data) 
rownames(pig_data@data) = 1:nrow(pig_data@data)
pig_data@data<-pig_data@data[,c("n.pigs", "mem.pigs", "intercept","urban", "protected","water.body", "elevation",
                                     "survey", "cluster", "district", "year")]
pig_data@data$lat<-pig_data@coords[2]
pig_data@data$lon<-pig_data@coords[1]

cattle_data@data$district<-as.character(cattle_data@data$district)
pig_data@data$district<-as.character(pig_data@data$district)

missing.years<-all.years[which(all.years %ni% cattle_data$year)]
for (year in missing.years){
  cattle_data = rbind(cattle_data[1,], cattle_data)
  cattle_data@data$survey[1]= 
    cattle_data@data$urban[1] = cattle_data@data$protected[1] =
    cattle_data@data$elevation[1] = cattle_data@data$water.body[1] = 
    cattle_data@data$mem.cattle[1] = cattle_data@data$district[1]=
    cattle_data@data$n.cattle[1]=NA
  cattle_data@data$year[1] = year
}

missing.years<-all.years[which(all.years %ni% pig_data$year)]
for (year in missing.years){
  pig_data = rbind(pig_data[1,], pig_data)
  pig_data@data$survey[1]= 
    pig_data@data$urban[1] = pig_data@data$protected[1] =
    pig_data@data$elevation[1] = pig_data@data$water.body[1] = 
    pig_data@data$mem.pigs[1]= pig_data@data$district[1]=
    pig_data@data$n.pigs[1]=NA
  pig_data@data$year[1] = year
}

cattle_data@data$time.unstruct = cattle_data@data$time = cattle_data@data$year - (min(cattle_data@data$year)-1)
pig_data@data$time.unstruct = pig_data@data$time = pig_data@data$year - (min(pig_data@data$year)-1)

mean.elev<-mean(cattle_data@data$elevation, na.rm=T)
sd.elev<-sd(cattle_data@data$elevation, na.rm=T)
cattle_data@data$elevation<-(cattle_data@data$elevation-mean.elev)/sd.elev

mean.elev<-mean(pig_data@data$elevation, na.rm=T)
sd.elev<-sd(pig_data@data$elevation, na.rm=T)
pig_data@data$elevation<-(pig_data@data$elevation-mean.elev)/sd.elev

polygon<- polygon_from_extent(districts)
cattle_data <- crop(cattle_data, polygon)
pig_data <- crop(pig_data, polygon)

writeOGR(obj=cattle_data, dsn="Data/Exposure_data/Malawi/Created_datasets/cattle_data", layer="cattle_data", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(obj=pig_data, dsn="Data/Exposure_data/Malawi/Created_datasets/pig_data", layer="pig_data", driver="ESRI Shapefile", overwrite_layer=TRUE)

districts_pigs<-pig_data@data$district
districts_cattle<-cattle_data@data$district
save(file="Data/Exposure_data/Malawi/Created_datasets/District_vector.Rdata", districts_cattle)
save(file="Data/Exposure_data/Malawi/Created_datasets/District_vector_pigs.Rdata", districts_pigs)

coords_cattle<-cattle_data@coords 
coords_pigs<-pig_data@coords