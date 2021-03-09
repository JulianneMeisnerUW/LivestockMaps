library(stringr)
library(stringi)
library(rgdal)
library(cobalt)
library(rgeos)
library(sp)
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
outcomeFolderData<-"Data/Outcome_data/"
districts<-readOGR("Data/Exposure_data/DRC/shapefiles/gadm36_COD_shp", layer="gadm36_COD_2") 
geographic<-proj4string(districts)

HAT<-readOGR(paste0("Data/Outcome_data/", country, "/HAT_regress"), layer="HAT_regress")
colnames(HAT@data)<-c("year", "focus", "cases", "density_med_c", 
                      "density_sd_c", "density_q025_c", "density_q975_c", "density_med_p", 
                      "density_sd_p", "density_q025_p", "density_q975_p","source", "prox_1km",
                      "cluster","intercept",  "WorldPop", "LandScan", "urban", "nightlights", 
                      "wealth_med", "wealth_sd", "wealth_q025", "wealth_975", "travel_hrs", "protected_5km")

'%ni%' <- Negate('%in%')
#-----------#
# Conflict
#-----------#
conflict<-read.csv("Data/Confounder_data/ged191.csv", header=T)
conflict_drc<-subset(conflict, country=="DR Congo (Zaire)")
conflict_drc$start<-as.numeric(str_sub(conflict_drc$date_start, end=4))
conflict_drc$end<-as.numeric(str_sub(conflict_drc$date_end, end=4))

#Convert to shapefile
conflict_coords<-cbind(conflict_drc$longitude, conflict_drc$latitude)
conflict.sp <- SpatialPointsDataFrame(coords = conflict_coords, data = conflict_drc,
                                      proj4string = CRS(geographic))

#Assign to Parish
merge.conflict<-over(conflict.sp, districts)
merge.HAT<-over(HAT, districts)
conflict.sp@data$adm2<-merge.conflict$NAME_2
HAT@data$adm2<-merge.HAT$NAME_2

study.adm2<-unique(HAT@data$adm2)

conflict.sp<-conflict.sp[which(conflict.sp@data$adm2%in%study.adm2),]

years<-2008:2015

HAT@data$conflict<-rep(0, nrow(HAT@data))

for(year in years){
  yIdx<-which(years==year)+1
  subset.conflict<-conflict.sp[which(conflict.sp@data$year==year),]
  adm2<-droplevels(subset.conflict$adm2)
  
  idx<-which(HAT@data$adm2%in%adm2 & HAT@data$year==year)
  subset.HAT<-HAT[idx,]
  
  HAT@data$conflict[idx]<-1
  print(year)
}

#-----------#
# Natural
#-----------#
natural<-read.csv("Data/Confounder_data/Julianne_Meisner_2019-11-15.csv", header=F)
names(natural) <- lapply(natural[2, ], as.character)
natural <- natural[-1,] 
natural <- natural[-1,] 
natural_drc<-subset(natural, Country=="Congo (the Democratic Republic of the)")

natural_drc$start<-str_sub(natural_drc$`Start date`, start= -4)
natural_drc$end<-str_sub(natural_drc$`End date`, start= -4)
natural_drc<-subset(natural_drc, end>=2008)

natural_drc$affected<-natural_drc$`Total affected`
natural_drc$affected<-ifelse(natural_drc$affected=="Total affected", NA, natural_drc$affected)
natural_drc$affected<-as.numeric(as.character(natural_drc$affected))
natural_drc<-natural_drc[which(natural_drc$affected>10),]
natural_drc$deaths<-natural_drc$`Total deaths`
natural_drc$type<-natural_drc$`Disaster type`
natural_drc$location<-droplevels(natural_drc$Location)

natural_drc$location<-as.character(natural_drc$location)

HAT@data$district<-merge.HAT$NAME_1
district<-unique(HAT@data$district)
district<-district[!is.na(district)]
district<-stri_trans_general(district, "latin-ascii")
HAT@data$district<-stri_trans_general(HAT@data$district, "latin-ascii")

natural_drc$district<-rep(NA, nrow(natural_drc))
Encoding(natural_drc$location) <- 'latin1'
natural_drc$location<-stri_trans_general(natural_drc$location, "latin-ascii")

keep=c()
for(i in 1:nrow(natural_drc)){
  loc=natural_drc$location[i]
  keep.district=c()
  for(j in 1:length(district)){
      r<-grepl(district[j], loc, fixed=TRUE)
      if(r=="TRUE"){
        keep.district=c(keep.district, i)
        natural_drc$district[i]<-as.character(district[j])
      }
  }
  if(length(keep.district)>0){
    keep=c(keep, i)
  }
}

natural_drc<-natural_drc[keep,]

HAT@data$disaster<-rep(NA, nrow(HAT@data))
HAT@data$dis_num<-rep(NA, nrow(HAT@data))
HAT@data$dis_deaths<-rep(NA, nrow(HAT@data))

HAT@data$interaction<-interaction(HAT@data$year, HAT@data$district)

for(i in 1:length(unique(HAT@data$interaction))){
  district_year<-unique(HAT@data$interaction)[i]
  year<-str_sub(district_year, end= 4)
  district<-str_sub(district_year, start=6)
  
  sub.HAT<-HAT[which(HAT@data$interaction==district_year),]
  
  sub<-natural_drc[which(natural_drc$start<=year&natural_drc$end>=year),]
  rows.district<-sub[which(grepl(district,sub$district)==TRUE),]
  
  if(nrow(rows.district)==0){
    HAT@data$disaster[which(HAT@data$interaction==district_year)]<-"None"
    HAT@data$dis_num[which(HAT@data$interaction==district_year)]<-"0"
    HAT@data$dis_deaths[which(HAT@data$interaction==district_year)]<-"0"
  }else{
    HAT@data$disaster[which(HAT@data$interaction==district_year)]<-as.character(rows.district$type)
    HAT@data$dis_num[which(HAT@data$interaction==district_year)]<-rows.district$affected
    HAT@data$dis_deaths[which(HAT@data$interaction==district_year)]<-rows.district$deaths
  }
}

writeOGR(obj=HAT, dsn="Data/Outcome_data/DRC/HATData_dis", layer="HATData_dis_gamb", driver="ESRI Shapefile", overwrite_layer=TRUE)
