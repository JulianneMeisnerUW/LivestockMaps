library(dplyr)
library(gdalUtils)
library(rgdal)
library(raster)
library(FedData)
outcomeFolderData<-"Data/Outcome_data/South_Sudan/"
confounderFolderData<-"Data/Confounder_data/South_Sudan/"

wealth<-readOGR(paste(confounderFolderData, "census_wealth", sep=""), layer="census_wealth")
wealth.active<-wealth
geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")

proj<-proj4string(wealth)

years=c(2000:2014)

#---------------------------------------------------------#
# Load outcome data and convert to spatial points dataframe
#---------------------------------------------------------#
HAT<-read.csv(paste0(outcomeFolderData,"HAT_South_Sudan_2000_2014.csv"))

HAT<-HAT[which(!is.na(HAT$Latitude)),]
HAT_coords<-cbind(HAT$Longitude, HAT$Latitude)

HAT.sp <- SpatialPointsDataFrame(coords = HAT_coords, data = HAT,
                                 proj4string = geographic)

#Assign to counties
merge.HAT<-over(HAT.sp, wealth)
HAT.sp@data$county<-merge.HAT$new_adm2

HAT.active<-HAT.sp[which(HAT.sp@data$Surveillance_type=="Active"),]

HAT.sp@data<-HAT.sp@data %>% dplyr::group_by(county, Year)%>%
  dplyr::summarise(cases=sum(New_HAT_cases,na.rm=T))%>% 
  dplyr::arrange(county, Year) %>%
  as.data.frame()
HAT.sp@data$county<-droplevels(HAT.sp@data$county)
study.counties<-unique(HAT.sp@data$county)

HAT.active@data<-HAT.active@data %>% group_by(county, Year)%>%
  summarise(cases=sum(New_HAT_cases,na.rm=T), 
            n=sum(People_screened, na.rm=T))%>% 
  arrange(county, Year) %>%
  as.data.frame()
HAT.active@data$county<-droplevels(HAT.active@data$county)
active.counties<-unique(HAT.active@data$county)

outcome.names<-paste("cs", 2000:2014, sep="")
outcome<-data.frame(matrix(nrow=nrow(wealth@data), ncol=length(outcome.names))) 
colnames(outcome)<-outcome.names 
rownames(outcome)<-wealth@data$new_adm2

#Make separate column for each year
for(year in years){
  yIdx<-which(years==year)
  subset<-HAT.sp@data[which(HAT.sp@data$Year==year),]
  for(i in 1:nrow(subset)){
    county<-subset$county[i]
    cIdx<-which(rownames(outcome)==county)
    outcome[cIdx,yIdx]<-subset$cases[i]
  }
}

outcome[is.na(outcome)] <- 0

#Merge with data
wealth@data<-cbind(wealth@data, outcome)

#Restrict to counties with positive cases over the study period
wealth<-wealth[which(wealth@data$new_adm2%in%study.counties),]

#Repeat for active surveillance only
outcome.names.active<-c(paste("n",2000:2014, sep=""),paste("cs",2000:2014, sep=""))
outcome.active<-data.frame(matrix(nrow=nrow(wealth@data), ncol=length(outcome.names.active))) 
colnames(outcome.active)<-outcome.names.active
rownames(outcome.active)<-wealth@data$new_adm2

#Make separate column for each year
for(year in years){
  yIdx<-which(years==year)
  subset<-HAT.active@data[which(HAT.active@data$Year==year),]
  for(i in 1:nrow(subset)){
    county<-subset$county[i]
    cIdx<-which(rownames(outcome.active)==county)
    outcome.active[cIdx,yIdx+15]<-subset$cases[i]
    outcome.active[cIdx,yIdx]<-subset$n[i]
  }
}

outcome.active[is.na(outcome.active)] <- 0

#Merge with data
wealth.active<-wealth.active[which(wealth.active@data$new_adm2%in%study.counties),]
wealth.active@data<-cbind(wealth.active@data, outcome.active)

#Restrict to counties with active surveillance cases over the study period
wealth.active<-wealth.active[which(wealth.active@data$new_adm2%in%active.counties),]
#-----------#
# Landscan
#-----------#
#for (year in 2000:2018){
#  dpath<-paste0("Data/denominator data/LandScan Global ",year, "/", "lspop", year)
#  x <- new("GDALReadOnlyDataset", dpath)
#  getDriver(x)
#  getDriverLongName(getDriver(x))
#  xx<-asSGDF_GROD(x)
#  r <- raster(xx)
#  landScancrop2<-crop(r, wealth)
#  filename<-paste0("Data/denominator data/LandScan/SS/LandScan_", year)
#  writeRaster(landScancrop2, filename, format="GTiff", overwrite=TRUE)
#}

LSpop.names<-paste("LS", 2000:2014, sep="")
LSpop<-data.frame(matrix(nrow=nrow(wealth@data), ncol=length(LSpop.names))) 
colnames(LSpop)<-LSpop.names 
rownames(LSpop)<-wealth@data$new_adm2

for(year in years){
  yIdx<-which(years==year)
  popRasterL = raster(paste('Data/Denominator_data/LandScan/SS/LandScan_', year, '.tif', sep = "")) 
  PL.crop<-crop(popRasterL, wealth)
  crs(PL.crop)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  trans.LP <- spTransform(wealth, crs(PL.crop))
  merge.LP<-extract(PL.crop, trans.LP, fun=sum, na.rm=TRUE)
  LS<-as.numeric(round(merge.LP,0))#Rounded to the nearest person
  LSpop[,yIdx]<-LS
}

#Merge with data
wealth@data<-cbind(wealth@data, LSpop)

#-----------#
# WorldPop
#-----------#
WPpop.names<-paste("WP", 2007:2014, sep="")
WPpop<-data.frame(matrix(nrow=nrow(wealth@data), ncol=length(WPpop.names))) 
colnames(WPpop)<-WPpop.names 
rownames(WPpop)<-wealth@data$new_adm2

WPyears<-c(2007:2014) #WorldPop rasters were corrupted for 2000-2006
for(year in WPyears){
  yIdx<-which(WPyears==year)
  popRasterW = raster(paste('Data/Denominator_data/WorldPop/SS/ssd_ppp_', year, '.tif', sep = "")) 
  PW.crop<-crop(popRasterW, wealth)
  trans.WP <- spTransform(wealth, crs(popRasterW))
  merge.WP<-extract(PW.crop, trans.WP, fun=sum, na.rm=TRUE)
  WP<-as.numeric(round(merge.WP,0)) 
  WPpop[,yIdx]<-WP
}

#Merge with data
wealth@data<-cbind(wealth@data, WPpop)

#-----------#
# Elevation
#-----------#
gmted2 <- raster("Data/Predictor_data/Elevation/10s000e_20101117_gmted_med075.tif") 
gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
gmted2.crop<-crop(gmted2, wealth) 
gmted6.crop<-crop(gmted6, wealth)

merged.trans <- spTransform(wealth, crs(gmted2.crop))
elevation<-merge(gmted2.crop, gmted6.crop)
elevation<-as.numeric(raster::extract(elevation, merged.trans, fun=mean, na.rm=T))
wealth@data$elevation<-elevation

merged.trans <- spTransform(wealth.active, crs(gmted2.crop))
elevation<-merge(gmted2.crop, gmted6.crop)
elevation<-as.numeric(raster::extract(elevation, merged.trans, fun=mean, na.rm=T))
wealth.active@data$elevation<-elevation

writeOGR(obj=wealth, dsn=paste0(outcomeFolderData, "SS_HAT"), layer="SS_HAT", driver="ESRI Shapefile", overwrite_layer = TRUE)

writeOGR(obj=wealth.active, dsn=paste0(outcomeFolderData, "SS_HAT"), layer="SS_HAT_active", driver="ESRI Shapefile", overwrite_layer = TRUE)
