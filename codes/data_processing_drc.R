library(rgdal)
library(sp)
library(raster)
library(foreign)
library(readstata13)
library(rgeos)
library(SpatialEpi)
library(survey)
library(RColorBrewer)
library(FedData)
library(INLA)
library(uwIntroStats)
library(dplyr)
library(geosphere)
source('codes/my_functions.R')
geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
districts<-readOGR("Data/Exposure_data/DRC/shapefiles/gadm36_COD_shp", layer="gadm36_COD_2") 
proj<-proj4string(districts)

#----------------------#
#   Read in data
#----------------------#
#######
# DHS
#######
dhs.data <- read.csv("Data/Exposure_data/Uganda/inputs/idhs_00016.csv")
dhs.data<-subset(dhs.data, COUNTRY==180)
dhs.data<-subset(dhs.data, YEAR==2013)

dhs.data$district<-as.factor(dhs.data$GEO_CD2007_2013)

dhs.data$HHID<-paste(dhs.data$DHSID, dhs.data$HHNUMALL)
dhs.data$URBANHH<-ifelse(dhs.data$URBANHH==2,0,1)

dhs.data$CATTLENUM<-ifelse(dhs.data$COWBULLNUM>95,NA, dhs.data$COWBULLNUM)
dhs.data$PIGNUM<-ifelse(dhs.data$PIGNUM>95,NA, dhs.data$PIGNUM)

#Collapse on household
dhs.hh<-dhs.data %>% group_by(HHID, YEAR)%>%
  summarise(num.id=n(), DHSCLUST= head(DHSID,1), HHMEMBERS= head(HHMEMBERS,1), urban=head(URBANHH,1),  
            district= head (district,1), cattle.num = mean(CATTLENUM, na.rm=T), 
            pig.num = mean(PIGNUM, na.rm=T))%>% 
  arrange(HHID)

dhs.hh<-as.data.frame(dhs.hh)

dhs.hh$hh.indicator<-rep(1, nrow(dhs.hh))

#Collapse on cluster
dhs.cluster<-dhs.hh %>% group_by(DHSCLUST, YEAR)%>%
  summarise(cluster.members=sum(HHMEMBERS, na.rm=T), hh.size= mean(HHMEMBERS, na.rm=T), households=n(),
            urban=head(urban,1), district=head(district,1), 
            num.cattle = sum(cattle.num, na.rm=T), 
            num.pigs = sum(pig.num, na.rm=T))%>% 
  arrange(DHSCLUST)
dhs.cluster<-as.data.frame(dhs.cluster)

dhs.cluster$clustermembers.cattle<-dhs.cluster$clustermembers.pigs<-dhs.cluster$cluster.members

#Density = # cattle/ # in in house (not #s surveyed!), summed across households
dhs.cluster$cattle.dens<-dhs.cluster$num.cattle/dhs.cluster$clustermembers.cattle
dhs.cluster$pig.dens<-dhs.cluster$num.pigs/dhs.cluster$clustermembers.pigs

DHS13.cluster<-dhs.cluster

DHS13.cluster$survey <-rep("DHS2013", nrow(DHS13.cluster))

###########
# MICS 2010
###########
mics4<-read.spss("Data/Exposure_data/DRC/inputs/MICS_2010_geo/SPSS_dataset/hh.sav", to.data.frame=TRUE)
mics4<-dplyr::rename(mics4, cluster_num= HH1, hh_num= HH2, livestockbin=HC13, cattle=HC14A, pigs=HC14F) 

mics4$urban<-ifelse(mics4$HH6=="Rural", 0,1)

mics4$num.cattle<-ifelse(mics4$cattle=="95+", 95,
                         ifelse(mics4$cattle=="NSP"|mics4$cattle=="Manquant", NA,
                                as.numeric(as.character(mics4$cattle))))
mics4$num.pigs<-ifelse(mics4$pigs=="95+", 95,
                       ifelse(mics4$pigs=="NSP"|mics4$pigs=="Manquant", NA,
                              as.numeric(as.character(mics4$pigs))))

mics4$hh.members<-mics4$HH11

#Set hh members to NA for hh with missing cattle, missing pigs
mics4$hhmembers.cattle<-ifelse(!is.na(mics4$num.cattle), mics4$hh.members, NA)
mics4$hhmembers.pigs<-ifelse(!is.na(mics4$num.pigs), mics4$hh.members, NA)

mics4.ea<-mics4 %>% group_by(cluster_num)%>%
  summarise(ea.members=sum(hh.members, na.rm=T), eamembers.cattle=sum(hhmembers.cattle, na.rm=T),
            eamembers.pigs=sum(hhmembers.pigs, na.rm=T), 
            num.cattle=sum(num.cattle,na.rm=T), num.pigs = sum(num.pigs, na.rm=T),
            urban= head(urban,1), 
            households=n())%>% 
  arrange(cluster_num)%>%as.data.frame()

mics4.ea$cattle.dens<-mics4.ea$num.cattle/mics4.ea$eamembers.cattle
mics4.ea$pig.dens<-mics4.ea$num.pigs/mics4.ea$eamembers.pigs

mics4.ea$survey <-rep("MICS2010", nrow(mics4.ea))
mics4<-mics4.ea

mics4geo<-read.csv("Data/Exposure_data/DRC/inputs/MICS_2010_geo/MICS-RDC 2010_GPS data_1.csv")
mics4geo<-dplyr::rename(mics4geo, cluster_num= Numero.de.la.grappe, lat=GP8, lon=GP9, village_name=Nom.quartier..village, urban= Milieu) 
latpaste<-as.numeric(paste0("-", mics4geo$lat))
mics4geo$lat<-ifelse(mics4geo$Nom.province=="Equateur"|mics4geo$Nom.province=="Nord-Kivu"|mics4geo$Nom.province=="Orientale", mics4geo$lat, latpaste)

mics4geo$issue<-ifelse(mics4geo$Nom.province=="Equateur"|mics4geo$Nom.province=="Nord-Kivu"|mics4geo$Nom.province=="Orientale",1,0)
to_solve<-subset(mics4geo, issue==1)

settlements<-readOGR("Data/Exposure_data/DRC/shapefiles/Settlements-shp", layer="Settlements")
settlements$nom_loc1<-toupper(settlements$nom_loc1)

solve.merge<-merge(to_solve[,c("village_name", "cluster_num")], settlements[,c("nom_loc1", "latitude", "longitude")], by.x="village_name", by.y="nom_loc1")
solve.merge$lat<-solve.merge$coords.x2
solve.merge$lon<-solve.merge$coords.x1
solve.merge<-solve.merge %>% group_by(village_name, cluster_num)%>%
  summarise(lat=mean(lat, na.rm=T),
            lon=mean(lon, na.rm=T))%>% 
  arrange(village_name, cluster_num)%>%as.data.frame()
solved.clusters<-solve.merge$cluster_num

for(i in 1:nrow(mics4geo)){
  num<-mics4geo$cluster_num[i]
  if(num%in%solved.clusters){
    idx<-which(solve.merge$cluster_num==num)
    mics4geo$lat[i]<-solve.merge$lat[idx]
    mics4geo$lon[i]<-solve.merge$lon[idx]
    mics4geo$issue[i]<-0
  }
}

mics4geo$district<-mics4geo$Nom.province

mics4<-merge(mics4.ea, mics4geo[,c("cluster_num", "lat", "lon", "district","issue")], by="cluster_num")

mics4<-mics4[!is.na(mics4$lon),]

############################
# Make column names match up
############################
DHS13.cluster$geotype<-rep("cluster", nrow(DHS13.cluster))
DHS13.cluster$year = 2013
DHS13.cluster$ea.members<-DHS13.cluster$cluster.members
DHS13.cluster$eamembers.cattle<-DHS13.cluster$clustermembers.cattle
DHS13.cluster$eamembers.pigs<-DHS13.cluster$clustermembers.pigs
DHS13<-DHS13.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
DHS13$baseline_urban<-rep(NA, length(nrow(DHS13)))

mics4.coords<-cbind(mics4$lon, mics4$lat)
mics4$geotype<-rep("cluster", nrow(mics4))
mics4$year = 2010
mics4<-mics4[,c("cluster_num", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]

##########################
# Shapefiles
##########################
shapedhs<-readOGR("Data/Exposure_data/DRC/inputs/DHS_2013/CDGE61FL", layer="CDGE61FL") 

geographic<-proj4string(districts)
polygon<- polygon_from_extent(districts,proj4string=geographic)
dcrop <- crop(shapedhs, polygon)

dhs13<-merge(dcrop, DHS13, by.x="DHSID", by.y="DHSCLUST")

proj<-proj4string(dhs13)
mics4.sp<-SpatialPointsDataFrame(coords = mics4.coords, data = mics4,
                                 proj4string = CRS(proj))

library(raster)

#Combine surveys
merged<-bind(dhs13, mics4.sp)
merged@data$cluster = rownames(merged@data)

##########################
# Add predictors
##########################
WCA1<-readOGR("Data/Exposure_data/DRC/inputs/WCA/DRC_SHAPEFILE", layer="DRCL1_01_CATTLE") 
WCA3<-readOGR("Data/Exposure_data/DRC/inputs/WCA/DRC_SHAPEFILE", layer="DRCL3_01_CATTLE_PIGS") 

proj<-proj4string(WCA3)

#Districts
merge1<-over(merged, WCA1)
merge1$province<-merge1$NAME
merge3<-over(merged, WCA3)
merge3$district<-merge3$NAME

merged@data$province<-merge1$province
merged@data$district<-merge3$district

idxNA = which (is.na(merged@data$district))
missing_district = merged[idxNA,]
plot(districts)
plot(missing_district, add=TRUE, col="red")

dists1 <- dist2Line(missing_district,WCA1)
dists3 <- dist2Line(missing_district,WCA3)
missing_district@data$district<-as.factor(WCA3@data$NAME[dists3[,4]])
missing_district@data$province<-as.factor(WCA1@data$NAME[dists1[,4]])

merged@data$district<-ifelse(is.na(merged@data$district), as.character(missing_district@data$district), as.character(merged@data$district))
merged@data$province<-ifelse(is.na(merged@data$province), as.character(missing_district@data$province), as.character(merged@data$province))

#Protected areas
protected<-readOGR("Data/Predictor_data/Protected_areas/WDPA_Jan2019-shapefile", layer="WDPA_Jan2019-shapefile-polygons", p4s=proj)
protected.crop<-crop(protected, polygon)
merge<-over(merged, protected.crop)
merge$protected<-ifelse(is.na(merge$WDPAID)==FALSE, 1,0)
merged@data$protected<-merge$protected

#Elevation
gmted2 <- raster("Data/Predictor_data/Elevation/10s000e_20101117_gmted_med075.tif") 
gmted3 <- raster("Data/Predictor_data/Elevation/30s030e_20101117_gmted_med075.tif") 
gmted4 <- raster("Data/Predictor_data/Elevation/30s000e_20101117_gmted_med075.tif") 
gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
gmted2.crop<-crop(gmted2, polygon) 
gmted3.crop<-crop(gmted3, polygon)
gmted4.crop<-crop(gmted4, polygon) 
gmted6.crop<-crop(gmted6, polygon)
merged.trans <- spTransform(merged, crs(gmted2.crop))
elevation.1<-raster::extract(gmted2.crop, merged.trans)
elevation.2<-raster::extract(gmted3.crop, merged.trans)
elevation.4<-raster::extract(gmted4.crop, merged.trans)
elevation.3<-raster::extract(gmted6.crop, merged.trans)
elev1<-ifelse(is.na(elevation.1)==TRUE, elevation.2, elevation.1)
elev2<-ifelse(is.na(elev1)==TRUE, elevation.3, elev1)
elev3<-ifelse(is.na(elev2)==TRUE, elevation.4, elev2)
merged@data$elevation<-elev3

#Lakes and wetlands
lakes1<-readOGR("Data/Predictor_data/Lakes_wetlands/GLWD-level1", layer="glwd_1", p4s=proj)
lakes2<-readOGR("Data/Predictor_data/Lakes_wetlands/GLWD-level2", layer="glwd_2", p4s=proj)
lakes1.crop<-crop(lakes1, polygon)
lakes2.crop<-crop(lakes2, polygon)
bind.lakes<-bind(lakes1.crop, lakes2.crop)
merge<-over(merged, bind.lakes)
merge$water.body<-ifelse(is.na(merge$TYPE)==FALSE, 1,0)
merged@data$water.body<-merge$water.body

merged@data$baseline_urban<-as.factor(merged@data$baseline_urban)
merged@data<-merged@data[,c("ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", 
                            "num.pigs", "survey", "year", "baseline_urban", "district", "cluster", "water.body",
                            "protected", "elevation")]
writeOGR(obj=merged, dsn="Data/Exposure_data/DRC/Created_datasets/merged_data", layer="merged_map", driver="ESRI Shapefile", overwrite_layer = TRUE)

#----------------------#
# Outcome data
#----------------------#
HAT<-read.csv("Data/Outcome_data/DRC/HAT_drc.csv")

HAT<-HAT[which(!is.na(HAT$Latitude)),]
HAT<-HAT[which(!is.na(HAT$Longitude)),]
HAT_coords<-cbind(HAT$Longitude, HAT$Latitude)

HAT$Province<-as.factor(HAT$Province)
HAT$Health_Zone<-as.factor(HAT$ZONE_DE_SANTE)
HAT$District<-as.factor(HAT$DISTRICT)
HAT$Territory<-as.factor(HAT$TERRITOIRE)

HAT.sp <- SpatialPointsDataFrame(coords = HAT_coords, data = HAT,
                                 proj4string = crs(geographic))

#(1) Protected areas
merge.protected<-over(HAT.sp, protected)
protectedv2<-ifelse(is.na(merge.protected$WDPAID)==FALSE, 1,0)

#(2) Elevation
merge.trans <- spTransform(HAT.sp, crs(gmted3))
gmted2.crop<-crop(gmted2, HAT.sp) 
gmted6.crop<-crop(gmted6, HAT.sp)
elevation.1<-raster::extract(gmted2.crop, merge.trans)
elevation.2<-raster::extract(gmted6.crop, merge.trans)
elev<-ifelse(is.na(elevation.1)==TRUE, elevation.2, elevation.1)

#(3) Lakes and wetlands
bind.lakes<-bind(lakes1, lakes2)
lakes.crop<-crop(bind.lakes, districts)
merge.proj<-over(HAT.sp, bind.lakes)
water.body<-ifelse(is.na(merge.proj$TYPE)==FALSE, 1,0)

HATsp.data<-HAT.sp@data
HATsp.data$protected<-protectedv2
HATsp.data$elevation<-elev
HATsp.data$water.body<-water.body
HATsp.data<-HATsp.data[,c("Longitude", "Latitude", "Location_name", "Year", "Date", "Census", "People_screened", "Surveillance_type", "New_HAT_cases", "P1", 
                  "P2", "P_na", "Parasite", "Health_Zone", "protected", "elevation", "water.body", "District", "Territory")]

HAT.sp <- SpatialPointsDataFrame(coords = coordinates(HAT.sp), data = HATsp.data,
                                   proj4string = CRS(proj))


writeOGR(obj=HAT.sp, dsn="Data/Outcome_data/DRC/HAT_Data", layer="HAT_data", driver="ESRI Shapefile", overwrite_layer=TRUE)

#----------------------#
#Random field prediction
#----------------------#
mapExt=extent(districts)
nx=1170
ny=1090
nPred = nx*ny
x = seq(mapExt@xmin-0.1, mapExt@xmax+0.1, length.out = nx)
xs = rep(x, each = ny)
y = seq(mapExt@ymin-0.1, mapExt@ymax+0.1, length.out = ny)
ys = rep(y, nx)
loc.pred = cbind(xs, ys)
xMat = matrix(xs, ncol = nx)
yMat = matrix(ys, ncol = nx)
gridP = data.frame(Longitude = as.vector(loc.pred[,1]),
                   Latitude = as.vector(loc.pred[,2]))
coordinates(gridP) = ~ Longitude + Latitude
proj4string(gridP) = proj4string(districts)

#(1) Protected areas
merge.protected<-over(gridP, protected)
protectedv2<-ifelse(is.na(merge.protected$WDPAID)==FALSE, 1,0)

#(2) Elevation
merge.trans <- spTransform(gridP, crs(gmted2))
elevation.1<-raster::extract(gmted2, merge.trans)
elevation.2<-raster::extract(gmted3, merge.trans)
elevation.3<-raster::extract(gmted4, merge.trans)
elevation.4<-raster::extract(gmted6, merge.trans)
elev1<-ifelse(is.na(elevation.1)==TRUE, elevation.2, elevation.1)
elev2<-ifelse(is.na(elev1)==TRUE, elevation.3, elev1)
elev3<-ifelse(is.na(elev2)==TRUE, elevation.4, elev2)
elevation<-elevation

#(3) Lakes and wetlands
bind.lakes<-bind(lakes1, lakes2)
lakes.crop<-crop(bind.lakes, districts)
merge.proj<-over(gridP, bind.lakes)
water.body<-ifelse(is.na(merge.proj$TYPE)==FALSE, 1,0)

gridP.data=data.frame(protected=protectedv2, 
                      water.body= water.body,
                      elevation =elevation)

gridP.sp <- SpatialPointsDataFrame(coords = coordinates(gridP), data = gridP.data,
                                   proj4string = CRS(proj))

writeOGR(obj=gridP.sp, dsn="Data/Exposure_data/DRC/Created_datasets/prediction_data", layer="prediction_data", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=protected.crop, dsn="Data/Exposure_data/Created_datasets/prediction_data_protected", layer="prediction_data_protected", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeRaster(gmted2.crop, "Data/Exposure_data/DRC/Created_datasets/prediction_data_elevation1", format="GTiff", overwrite=TRUE)
writeRaster(gmted6.crop, "Data/Exposure_data/DRC/Created_datasets/prediction_data_elevation2", format="GTiff", overwrite=TRUE)
writeOGR(obj=lakes.crop, dsn="Data/Exposure_data/DRC/Created_datasets/prediction_data_lakes", layer="prediction_data_lakes", driver="ESRI Shapefile", overwrite_layer=TRUE)

#---------------#
#  LandScan
#---------------#
library(raster)
for (year in 2000:2018){
  dpath<-paste0("Data/denominator data/LandScan Global ",year, "/", "lspop", year)
  x <- new("GDALReadOnlyDataset", dpath)
  getDriver(x)
  getDriverLongName(getDriver(x))
  xx<-asSGDF_GROD(x)
  r <- raster(xx)
  landScancrop2<-crop(r, districts)
  filename<-paste0("Data/denominator data/LandScan/DRC/LandScan_", year)
  writeRaster(landScancrop2, filename, format="GTiff", overwrite=TRUE)
}