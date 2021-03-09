library(stringr)
library(rgdal)
library(cobalt)
library(sp)
'%ni%' <- Negate('%in%')
outcomeFolderData<-"Data/Outcome Data/"
folderData<-"Data/Exposure Data/South Sudan/output/"

HAT<-readOGR("HAT_RS_SS", layer="HAT_RS_SS")
HAT.active<-readOGR("HAT_RS_SS", layer="HAT_RS_SS_active")

geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")

proj<-proj4string(HAT)

census<-readOGR(paste(folderData, "census_data", sep=""), layer="census_data")
#-----------#
# Conflict
#-----------#
conflict<-read.csv("Data/Confounder data/ged191.csv", header=T)#Global
conflict_SS<-conflict[which(conflict$country=="South Sudan"|conflict$country=="Sudan"),]
conflict_SS$start<-as.numeric(str_sub(conflict_SS$date_start, end=4))#Just do time in years
conflict_SS$end<-as.numeric(str_sub(conflict_SS$date_end, end=4))#Just do time in years

#Convert to shapefile
conflict_coords<-cbind(conflict_SS$longitude, conflict_SS$latitude)
conflict.sp <- SpatialPointsDataFrame(coords = conflict_coords, data = conflict_SS,
                                 proj4string = geographic)

#Assign to counties
merge.conflict<-over(conflict.sp, census)
conflict.sp@data$county<-merge.conflict$new_adm2

study.counties<-HAT@data$new_adm2
study.states<-unique(HAT@data$ADM1_EN)

conflict.sp<-conflict.sp[which(conflict.sp@data$county%in%study.counties),]

years<-2000:2014

conflict.names<-paste("war", 2000:2014, sep="")
conflict<-data.frame(matrix(nrow=nrow(HAT@data), ncol=length(conflict.names))) 
colnames(conflict)<-conflict.names 
rownames(conflict)<-HAT@data$new_adm2

conflict.active<-data.frame(matrix(nrow=nrow(HAT.active@data), ncol=length(conflict.names))) 
colnames(conflict.active)<-conflict.names 
rownames(conflict.active)<-HAT.active@data$new_adm2

for(year in years){
  yIdx<-which(years==year)
  subset<-conflict.sp[which(conflict.sp@data$year==year),]
  tab<-table(subset@data$county)
  for(i in 1:dim(tab)){
    county<-names(tab)[i]
    cIdx<-which(rownames(conflict)==county)
    if(tab[i]>0){
      conflict[cIdx,yIdx]<-1
    }else{
      conflict[cIdx,yIdx]<-0
    }
  }
}

HAT@data<-cbind(HAT@data, conflict)

for(year in years){
  yIdx<-which(years==year)
  subset<-conflict.sp[which(conflict.sp@data$year==year),]
  tab<-table(subset@data$county)
  for(i in 1:dim(tab)){
    county<-names(tab)[i]
    cIdx<-which(rownames(conflict.active)==county)
    if(tab[i]>0){
      conflict.active[cIdx,yIdx]<-1
    }else{
      conflict.active[cIdx,yIdx]<-0
    }
  }
}

HAT.active@data<-cbind(HAT.active@data, conflict.active)

#-----------#
# Natural
#-----------#
natural<-read.csv("Data/Confounder data/emdat_public_2020_03_03.csv", header=T)
natural_SS<-natural[which(natural$Country=="Sudan (the)"|natural$Country=="South Sudan"),]

natural_SS$start<-natural_SS$Start.year#Just do time in years
natural_SS$end<-natural_SS$End.year#Just do time in years
natural_SS<-subset(natural_SS, end>=2000)#Subset to study years

natural_SS$affected<-natural_SS$Total.affected
natural_SS$deaths<-natural_SS$Total.deaths
natural_SS$type<-natural_SS$Type

remove<-which(natural_SS$Start.year<2008&natural_SS$Country.name=="Sudan (the)") #Started recording South Sudan in 2008
natural_SS<-natural_SS[-remove,]

#Subset data that have coordinate values
natural_SS$Latitude<-as.numeric(as.character(natural_SS$Latitude))
have_coord<-natural_SS[which(!is.na(natural_SS$Latitude)),]
Lat<-as.numeric(as.character(have_coord$Latitude))
Lon<-as.numeric(as.character(have_coord$Longitude))
natural_coords<-cbind(Lon, Lat)
natural.sp <- SpatialPointsDataFrame(coords = natural_coords, data = have_coord,
                                     proj4string = geographic)

merge.natural<-over(natural.sp, census)
natural.sp@data$Location<-merge.natural$new_adm2

#And data without coordinate values
match<-natural_SS[which(is.na(natural_SS$Latitude)),]
location<-match$Location
new_location<-rep(NA, 73)
new_location[1]<-"Bor, Nahr Atiem districts (Jonglei province), Aweil district (Northern Bahr el Ghazal province), Sobat district (Upper Nile province), Warrap province"                                                                       
new_location[2]<-"Northern Bahr el Ghazal, Lakes, Jonglei, Upper Nile, Warrap provincess"
new_location[3]<-"Sudan"
new_location[4]<-"Aweil Centre"
new_location[5]<-"Gogrial West, Gogrial East"
new_location[6]<-"Sudan"
new_location[7]<-NA
new_location[8]<-NA
new_location[9]<-"Sudan"
new_location[10]<-"Sudan"
new_location[11]<-"Sudan"
new_location[12]<-"Unity, Northern Bahr el Ghazal, Jonglei, Upper Nile, Eastern Equatoria, Warrap provinces"
new_location[13]<-"Sudan"
new_location[14]<-"Jonglei, Lakes, Warrap, Northern Bahr el Ghazal, Unity, Central Equatoria, Upper Nile provinces"
new_location[15]<-"Jonglei, Upper Nile, Unity, Northern Bahr el Ghazal provinces"
new_location[16]<-"Sudan"
new_location[17]<-"Sudan"
new_location[18]<-"Sudan"
new_location[19]<-"Sudan"
new_location[20]<-"Sudan"
new_location[21]<-"Sudan"
new_location[22]<-"Sudan"
new_location[23]<-"Sudan"
new_location[24]<-"Abyei"
new_location[25]<-"Sudan"
new_location[26]<-"Sudan"
new_location[27]<-"Sudan"
new_location[28]<-"Sudan"
new_location[29]<-"Sudan"
new_location[30]<-"Sudan"
new_location[31]<-"Sudan"
new_location[32]<-"Sudan"
new_location[33]<-"Northern Bahr el Ghazal, Western Bahr el Ghazal, Warrap"
new_location[34]<-"Northern Bahr el Ghazal, Warrap, Upper Nile, Western Bahr el Ghazal"
new_location[35]<-NA
new_location[36]<-"Sudan"
new_location[37]<-"Sudan"
new_location[38]<-"Sudan"
new_location[39]<-"Sudan"
new_location[40]<-"Sudan"
new_location[41]<-"Sudan"
new_location[42]<-"Sudan"
new_location[43]<-"Sudan"
new_location[44]<-"Juba, Jonglei, Upper Nile"
new_location[45]<-"Juba, Torit"
new_location[46]<-"Sudan"
new_location[47]<-"Juba"
new_location[48]<-"Sudan"
new_location[49]<-"Sudan"
new_location[50]<-"Sudan"
new_location[51]<-"Sudan"
new_location[52]<-"Sudan"
new_location[53]<-"Juba, Bor"
new_location[54]<-"Juba"
new_location[55]<-"Maridi"
new_location[56]<-"Western Bahr el Ghazal, Northern Bahr el Ghazal, Warrap, Unity, Upper Nile, Jonglei, Lakes, Western Equatoria, Central Equatoria, Eastern Equatoria"
new_location[57]<-"Sudan"
new_location[58]<-"Sudan"
new_location[59]<-"Sudan"
new_location[60]<-NA
new_location[61]<-"Torit, Magwi,  Awerial, Juba, Terekeka, Duk, Fangak, Rubkona, Leer, Pigi (Eastern Nile)"
new_location[62]<-"Unity"
new_location[63]<-"Mayom (Unity)"
new_location[64]<-"Sudan"
new_location[65]<-"Yirol"
new_location[66]<-"Sudan"
new_location[67]<-"Sudan"
new_location[68]<-"Lafon, Torit and Kapoeta South counties (Eastern Equatoria); Ayod, Akobo, Bor South, Duk, Twic East, Pibor, Pochalla and Uror counties (Jonglei); Aweil Center, Aweil North (Northern Bahr el Ghazal); Abiemnhom, Mayom, Mayendit, Panyijiar (Unity); Maban (Upper Nile);  Gogrial East, Gogrial West, Tonj North (Warrap); Juba, Terekeka (Central Equatoria)"
new_location[69]<-"Sudan"
new_location[70]<-"Sudan"
new_location[71]<-"Sudan"
new_location[72]<-"Sudan"
new_location[73]<-"Sudan"

match$Location<-new_location

#Bind together
natural_SS<-rbind(match, natural.sp@data)

#Remove Sudan and missing location
natural_SS<-natural_SS[which(natural_SS$Location!="Sudan"),]
natural_SS<-natural_SS[which(!is.na(natural_SS$Location)),]

county=c()
state=c()
keep.county=c()
keep.state=c()
for(i in 1:nrow(natural_SS)){
  loc=natural_SS$Location[i]
  for(j in 1:length(study.counties)){
    if(grepl(study.counties[j], loc, fixed=TRUE)==TRUE){
      print(paste("Match for county",study.counties[j], ", data row", i, sep=" "))
      keep.county=c(keep.county,i)
      county=c(county,as.character(study.counties[j]))
    }
  }
  for(k in 1:length(study.states)){
    if(grepl(study.states[k], loc, fixed=TRUE)==TRUE){
      print(paste("Match for state", study.states[k], ", data row", i, sep=" "))
      keep.state=c(keep.state, i)
      state=c(state,as.character(study.states[k]))
    }
  }
}

keep<-which(keep.state%ni%keep.county)
keep.state[keep]
state[keep]

natural_SS$og_row<-1:nrow(natural_SS)

natural.names<-paste("nat", 2000:2014, sep="")
natural<-data.frame(matrix(nrow=nrow(HAT@data), ncol=length(natural.names))) 
colnames(natural)<-natural.names 
rownames(natural)<-HAT@data$new_adm2
natural$county<-HAT@data$new_adm2
natural$state<-HAT@data$ADM1_EN

for(year in years){
  yIdx<-which(years==year)
  subset<-natural_SS[which(natural_SS$Start.year<=year&natural_SS$End.year>=year),]
  if(nrow(subset)>0){
    for(i in 1:nrow(subset)){
      disaster<-subset$type[i]
      index<-subset$og_row[i]
      if(index%in%keep.county|index%in%keep.state){
        if(index%in%keep.county){
          idx<-which(keep.county==index)
          location<-county[idx]
          cIdx<-which(natural$county%in%location)
          natural[cIdx,yIdx]<-as.character(disaster)
        }else{
          idx<-which(keep.state==index)
          location<-state[idx]
          sIdx<-which(natural$state%in%location)
          natural[sIdx,yIdx]<-as.character(disaster)
        }
      }
    }
  }
}

natural[is.na(natural)] <- "None"
remove<-which(names(natural)%in%c("county", "state"))
natural<-natural[,-remove]

HAT@data<-cbind(HAT@data, natural)

natural.active<-data.frame(matrix(nrow=nrow(HAT.active@data), ncol=length(natural.names))) 
colnames(natural.active)<-natural.names 
rownames(natural.active)<-HAT.active@data$new_adm2
natural.active$county<-HAT.active@data$new_adm2
natural.active$state<-HAT.active@data$ADM1_EN

for(year in years){
  yIdx<-which(years==year)
  subset<-natural_SS[which(natural_SS$Start.year<=year&natural_SS$End.year>=year),]
  if(nrow(subset)>0){
    for(i in 1:nrow(subset)){
      disaster<-subset$type[i]
      index<-subset$og_row[i]
      if(index%in%keep.county|index%in%keep.state){
        if(index%in%keep.county){
          idx<-which(keep.county==index)
          location<-county[idx]
          cIdx<-which(natural.active$county%in%location)
          natural.active[cIdx,yIdx]<-as.character(disaster)
        }else{
          idx<-which(keep.state==index)
          location<-state[idx]
          sIdx<-which(natural.active$state%in%location)
          natural.active[sIdx,yIdx]<-as.character(disaster)
        }
      }
    }
  }
}

natural.active[is.na(natural.active)] <- "None"
remove<-which(names(natural.active)%in%c("county", "state"))
natural.active<-natural.active[,-remove]

HAT.active@data<-cbind(HAT.active@data, natural)

writeOGR(obj=HAT, dsn="HAT_final_SS", layer="HAT_final_SS", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(obj=HAT.active, dsn="HAT_final_SS", layer="HAT_final_SS_active", driver="ESRI Shapefile", overwrite_layer=TRUE)


