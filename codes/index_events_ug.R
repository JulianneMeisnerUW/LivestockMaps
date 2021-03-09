library(stringr)
library(rgdal)
library(cobalt)
library(rgeos)
library(sp)
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
districts<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_4")
geographic<-proj4string(districts)

#Read in outcome data to figure out which locations I need
HAT.r<-readOGR("Data/Outcome_data/Uganda/HAT_regress", layer="HAT_regress_rhod", p4s=geographic)
colnames(HAT.r@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                      "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                      "density_q025_p", "density_q975_p", "source","prox1km", "cluster", "intercept",
                      "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                      "wealth_975", "district", "travel_hrs","protected")

HAT.g<-readOGR("Data/Outcome_data/Uganda/HAT_regress", layer="HAT_regress_gamb", p4s=geographic)
colnames(HAT.g@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                        "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                        "density_q025_p", "density_q975_p", "source","prox1km", "cluster", "intercept",
                        "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                        "wealth_975", "district", "travel_hrs","protected")

#-----------#
# Conflict
#-----------#
conflict<-read.csv("Data/Confounder_data/ged191.csv", header=T)
conflict_uganda<-subset(conflict, country=="Uganda")
conflict_uganda$start<-as.numeric(str_sub(conflict_uganda$date_start, end=4))
conflict_uganda$end<-as.numeric(str_sub(conflict_uganda$date_end, end=4))

#Convert to shapefile
conflict_coords<-cbind(conflict_uganda$longitude, conflict_uganda$latitude)
conflict.sp <- SpatialPointsDataFrame(coords = conflict_coords, data = conflict_uganda,
                                      proj4string = CRS(geographic))

#Assign to Parish
merge.conflict<-over(conflict.sp, districts)
merge.HAT.g<-over(HAT.g, districts)
merge.HAT.r<-over(HAT.r, districts)
conflict.sp@data$parish<-merge.conflict$NAME_4
HAT.g@data$parish<-merge.HAT.g$NAME_4
HAT.r@data$parish<-merge.HAT.r$NAME_4
  
study.parish.g<-unique(HAT.g@data$parish)
study.parish.r<-unique(HAT.r@data$parish)

conflict.sp.g<-conflict.sp[which(conflict.sp@data$parish%in%study.parish.g),]
conflict.sp.r<-conflict.sp[which(conflict.sp@data$parish%in%study.parish.r),]

years<-2000:2018

conflict.names<-paste("war", 2000:2018, sep="")
conflict.g<-data.frame(matrix(nrow=nrow(HAT.g@data), ncol=length(conflict.names)+1)) 
colnames(conflict.g)<-c("parish", conflict.names) 
conflict.g$parish<-HAT.g@data$parish

conflict.r<-data.frame(matrix(nrow=nrow(HAT.r@data), ncol=length(conflict.names)+1)) 
colnames(conflict.r)<-c("parish", conflict.names) 
conflict.r$parish<-HAT.r@data$parish

for(year in years){
  yIdx<-which(years==year)
  subset.g<-conflict.sp.g[which(conflict.sp.g@data$year==year),]
  subset.g$parish<-droplevels(subset.g$parish)
  tab.g<-table(subset.g@data$parish)
  
  subset.r<-conflict.sp.r[which(conflict.sp.r@data$year==year),]
  subset.r$parish<-droplevels(subset.r$parish)
  tab.r<-table(subset.r@data$parish)
  
  if(dim(tab.g)==0){
    conflict.g[,yIdx]<-0
  }else{
    for(i in 1:dim(tab.g)){
      parish<-names(tab.g)[i]
      cIdx.g<-which(conflict.g$parish==parish)
      if(tab.g[i]>0){
        conflict.g[cIdx.g,yIdx]<-1
      }else{
        conflict.g[cIdx.g,yIdx]<-0
      }
    }
  }
  
  if(dim(tab.r)==0){
    conflict.r[,yIdx]<-0
  }else{
    for(i in 1:dim(tab.r)){
      parish<-names(tab.r)[i]
      cIdx.r<-which(conflict.r$parish==parish)
      if(tab.r[i]>0){
        conflict.r[cIdx.r,yIdx]<-1
      }else{
        conflict.r[cIdx.r,yIdx]<-0
      }
    }
  }
  print(year)
}
conflict.g[is.na(conflict.g)] = 0
conflict.r[is.na(conflict.r)] = 0
conflict.g<-conflict.g[,-1]
conflict.r<-conflict.r[,-1]

#Merge with HAT data
for(i in 1:nrow(HAT.g@data)){
  year=HAT.g@data$year[i]
  yIdx=which(years==year)
  conflict<-conflict.g[,yIdx]
  HAT.g@data$conflict[i]<-conflict[i]
}
#Merge with HAT data
for(i in 1:nrow(HAT.r@data)){
  year=HAT.r@data$year[i]
  yIdx=which(years==year)
  conflict<-conflict.r[,yIdx]
  HAT.r@data$conflict[i]<-conflict[i]
}

#-----------#
# Natural
#-----------#
natural<-read.csv("Data/Confounder_data/Julianne_Meisner_2019-11-15.csv", header=F)
names(natural) <- lapply(natural[2, ], as.character)
natural <- natural[-1,] 
natural <- natural[-1,] 
natural_uganda<-subset(natural, Country=="Uganda")

natural_uganda$start<-str_sub(natural_uganda$`Start date`, start= -4)#Just do time in years
natural_uganda$end<-str_sub(natural_uganda$`End date`, start= -4)
natural_uganda<-subset(natural_uganda, end>=2000)#Subset to study years

natural_uganda$affected<-natural_uganda$`Total affected`
natural_uganda$affected<-ifelse(natural_uganda$affected=="Total affected", NA, natural_uganda$affected)
natural_uganda$affected<-as.numeric(as.character(natural_uganda$affected))
natural_uganda<-natural_uganda[which(natural_uganda$affected>10),]#subset to natural disasters affecting 10 or more peoples
natural_uganda$deaths<-natural_uganda$`Total deaths`
natural_uganda$type<-natural_uganda$`Disaster type`
natural_uganda$location<-droplevels(natural_uganda$Location)

natural_uganda$location<-as.character(natural_uganda$location)

#--------------------------------------------#
#Subset to districts with cases
#--------------------------------------------#
HAT.g@data$district<-merge.HAT.g$NAME_1
district.g<-unique(HAT.g@data$district)
HAT.r@data$district<-merge.HAT.r$NAME_1
district.r<-unique(HAT.r@data$district)

natural_uganda$district<-rep(NA, nrow(natural_uganda))

keep.g=c()
keep.r=c()
for(i in 1:nrow(natural_uganda)){
  loc=natural_uganda$location[i]
  keep.district.g=keep.district.r=c()
  for(j in 1:length(district.g)){
      r<-grepl(district.g[j], loc, fixed=TRUE)
      if(r=="TRUE"){
        keep.district.g=c(keep.district.g, i)
        natural_uganda$district[i]<-as.character(district.g[j])
    }
  }
  if(length(keep.district.g)>0){
    keep.g=c(keep.g, i)
  }
  
  for(j in 1:length(district.r)){
    r<-grepl(district.r[j], loc, fixed=TRUE)
    if(r=="TRUE"){
      keep.district.r=c(keep.district.r, i)
      natural_uganda$district[i]<-as.character(district.r[j])
    }
  }
  if(length(keep.district.r)>0){
    keep.r=c(keep.r, i)
  }
}

natural_uganda.g<-natural_uganda[keep.g,]
natural_uganda.r<-natural_uganda[keep.r,]

#Check things
HAT.g@data$disaster<-rep(NA, nrow(HAT.g@data))
HAT.g@data$dis_num<-rep(NA, nrow(HAT.g@data))
HAT.g@data$dis_deaths<-rep(NA, nrow(HAT.g@data))

HAT.r@data$disaster<-rep(NA, nrow(HAT.r@data))
HAT.r@data$dis_num<-rep(NA, nrow(HAT.r@data))
HAT.r@data$dis_deaths<-rep(NA, nrow(HAT.r@data))

for(i in 1:nrow(HAT.g@data)){
  district<-HAT.g@data$district[i]
  year<-HAT.g@data$year[i]
  sub<-natural_uganda.g[which(natural_uganda.g$start<=year&natural_uganda.g$end>=year),]
  rows.district<-sub[which(grepl(district,sub$district)==TRUE),]
  if(nrow(rows.district)==0){
    HAT.g@data$disaster[i]<-"None"
    HAT.g@data$dis_num[i]<-0
    HAT.g@data$dis_deaths[i]<-0
  }else{
    HAT.g@data$disaster[i]<-as.character(rows.district$type)
    HAT.g@data$dis_num[i]<-rows.district$affected
    HAT.g@data$dis_deaths[i]<-rows.district$deaths
  }
}

for(i in 1:nrow(HAT.r@data)){
  district<-HAT.r@data$district[i]
  year<-HAT.r@data$year[i]
  sub<-natural_uganda.r[which(natural_uganda.r$start<=year&natural_uganda.r$end>=year),]
  rows.district<-sub[which(grepl(district,sub$district)==TRUE),]
  if(nrow(rows.district)==0){
    HAT.r@data$disaster[i]<-"None"
    HAT.r@data$dis_num[i]<-0
    HAT.r@data$dis_deaths[i]<-0
  }else{
    HAT.r@data$disaster[i]<-as.character(rows.district$type)
    HAT.r@data$dis_num[i]<-rows.district$affected
    HAT.r@data$dis_deaths[i]<-rows.district$deaths
  }
}

writeOGR(obj=HAT.g, dsn="Data/Outcome_data/Uganda/HATData_dis", layer="HATData_dis_gamb", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(obj=HAT.r, dsn="Data/Outcome_data/Uganda/HATData_dis", layer="HATData_dis_rhod", driver="ESRI Shapefile", overwrite_layer=TRUE)
