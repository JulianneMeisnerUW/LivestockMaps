library(readstata13)
library(dplyr)
library(mlr)
library(uwIntroStats)
library(psych)
library(Hmisc)
library(gtools)
library(questionr)
library(foreign)
library(geosphere)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(SpatialEpi)
library(survey)
library(RColorBrewer)
library(FedData)
library(INLA)
library(dplyr)
library(ipumsr)
library(sp)
library(rgdal)

source('codes/my_functions.R')

binary<-function(x){
  if(nlevels(as.factor(x))>2){
    ifelse(x==0,0,
           ifelse(x>=1,1,NA))
  }else{
    ifelse(x==1,1,0)
  }
}

ridmiss<-function(x,mIdx,yIdx){
  if(nlevels(as.factor(x))>2){
    ifelse(x==mIdx,NA, 
           ifelse(x==yIdx,1,0))
  }else{
    ifelse(x==1,1,0)
  }
}

geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
districts<-readOGR("Data/Exposure_data/DRC/inputs/shapefiles/gadm36_COD_shp", layer="gadm36_COD_2") 
proj<-proj4string(districts)
WCA1<-readOGR("Data/Exposure_data/DRC/shapefiles/WCA/DRC_SHAPEFILE", layer="DRCL1_01_CATTLE") 
WCA3<-readOGR("Data/Exposure_data/DRC/shapefiles/WCA/DRC_SHAPEFILE", layer="DRCL3_01_CATTLE_PIGS") 

nonzeroSD <- function(vector) {
  sd<-sd(vector, na.rm=T)
  is.nonzero<-ifelse(sd>0.001,TRUE,FALSE)
  return(is.nonzero)
}

propmiss <- function(dataframe) {
  m <- sapply(dataframe, function(x) {
    data.frame(
      nmiss=sum(is.na(x)), 
      n=length(x), 
      propmiss=sum(is.na(x))/length(x)
    )
  })
  d <- data.frame(t(m))
  d <- sapply(d, unlist)
  d <- as.data.frame(d)
  d$variable <- row.names(d)
  row.names(d) <- NULL
  d <- cbind(d[ncol(d)],d[-ncol(d)])
  return(d[order(d$propmiss), ])
}

create_var<-function(x, dat, y, id){
  idx<-which(names(dat)==id)
  nRow<-length(unique(dat[,idx]))
  finaldat<-data.frame(matrix(nrow=nRow, ncol=2))
  for(i in 1:length(x)){
    if(i==1){
      row<-dat[which(dat$goods==x[i]),]
      cIdx<-which(names(row)==y)
      newvar<-ifelse(row[,cIdx]==1,1,0)
      idx<-which(names(row)==id)
      newdat<-as.data.frame(cbind(row[,idx], newvar))
      names(newdat)<-c(id, x[i])
      newdat[,2]<-as.integer(newdat[,2])
      newdat.hh<-aggregate(newdat[,2], by=list(newdat[,1]), FUN=max, na.rm=TRUE)
      finaldat[,i]<-newdat.hh[,1]
      finaldat[,i+1]<-newdat.hh[,2]
      colnames(finaldat)<-colnames(newdat)
    }else{
      row<-dat[which(dat$goods==x[i]),]
      cIdx<-which(names(row)==y)
      newvar<-ifelse(row[,cIdx]==1,1,0)
      idx<-which(names(row)==id)
      newdat<-as.data.frame(cbind(row[,idx], newvar))
      names(newdat)<-c(id, x[i])
      newdat[,2]<-as.integer(newdat[,2])
      newdat.hh<-aggregate(newdat[,2], by=list(newdat[,1]), FUN=max, na.rm=TRUE)
      colnames(newdat.hh)<-colnames(newdat)
      finaldat<-merge(finaldat, newdat.hh, by="HHID")
    }
  }
  return(finaldat)
}

'%ni%' <- Negate('%in%')

needed<-c("cluster", "survey", "year", "domestic", "land", "house", "urban", "fuel", "water", "toilet", "toilet_share", "floor", "walls", 
          "roof", "electricity", "bank", "landarea", "chairs", "vehicle", "clock",  
           "fridge", "mobile", "landline", "tv", "computer", "radio", "sewingmachine", "canoe",
          "bed",  "lamp", "boat", "stove_eg", "generator", "axehoe", "possessions", "district")

################################################
# DHS
################################################
dhs.data <- read.csv("Data/Exposure_data/Uganda/inputs/DHS/idhs_00016.csv")
dhs.data<-subset(dhs.data, COUNTRY==180)

dhs.data$district<-as.factor(dhs.data$GEO_CD2007_2013)
dhs.data$HHID<-paste(dhs.data$DHSID, dhs.data$HHNUMALL)
dhs.data$urban<-ifelse(dhs.data$URBANHH==2,0,1)

dhs.data$members<-dhs.data$DEJURENO
dhs.data$rooms<-ifelse(dhs.data$SLEEPROOMS==0,1, ifelse(dhs.data$SLEEPROOMS==97|dhs.data$SLEEPROOMS==98, NA, dhs.data$SLEEPROOMS))
dhs.data$bank<-ifelse(dhs.data$BANKACC==1,1,
                 ifelse(dhs.data$BANKACC==8,NA, 0))
dhs.data$survey<-paste("DHS", dhs.data$YEAR, sep="")
dhs.data$year<-dhs.data$YEAR
dhs.data$domestic<-ifelse(dhs.data$HHRELATE==41,1,0)
dhs.data$land<-ifelse(dhs.data$AGLANDYN==1,1, 
                 ifelse(dhs.data$AGLANDYN==8,NA, 0))
dhs.data$toilet_share<-ifelse(dhs.data$TOILETSHAREYN==1,1,
                 ifelse(dhs.data$TOILETSHAREYN==8|dhs.data$TOILETSHAREYN==9,NA, 0))
dhs.data$electricity<-ifelse(dhs.data$ELECTRCHH==1,1,
                         ifelse(dhs.data$ELECTRCHH==8|dhs.data$ELECTRCHH==9,NA, 0))
dhs.data$cluster<-dhs.data$DHSID

dhs.data$fuel<-as.factor(ifelse(dhs.data$COOKFUEL==800|dhs.data$COOKFUEL==995|dhs.data$COOKFUEL==998,NA, dhs.data$COOKFUEL))
levels(dhs.data$fuel)<-c("100: Electricity","240: Kerosene", "300: Biogas",
                     "410: Coal, lignite", "510: Wood", "520: Charcoal", "540: Straw_shrub_grass",
                    "600: Dung","700: Agricultural_crop_based")
levels(dhs.data$fuel)<-c(4,3,3,2,1,2,0,0,0)

dhs.data$toilet<-as.factor(ifelse(dhs.data$TOILETTYPE==5000|dhs.data$TOILETTYPE==9998,NA, dhs.data$TOILETTYPE))
levels(dhs.data$toilet)<-c("0: None", "1100: Flush, unspecified", "1210: Flush to piped sewer system",
                           "1250: Flush, septic tank", "1410: Flush to pit latrine", "1420: Flush to somewhere else",
                           "1430: Flush, don't know where", "2100: Composting toilet", "3210: Pit latrine, no slab",
                           "3300: Pit latrine, with slab", "3400: VIP", "4100: Bucket toilet", "4300: Hanging")
levels(dhs.data$toilet)<-c(0,3,3,3,3,3,3,2,1,2,2,1,1)

dhs.data$water<-as.factor(ifelse(dhs.data$DRINKWTR==6000|dhs.data$DRINKWTR==2310|dhs.data$DRINKWTR==2320|dhs.data$DRINKWTR==9998,NA, dhs.data$DRINKWTR))
levels(dhs.data$water)<-c("1110: Piped into own dwelling", "1120: Piped into own yard/plot", 
                          "1210: Public tab/standpipe", "1220: Piped into neighbor's dwelling/yard", 
                          "2100: Unprotected/open well", "2111: Open well in own dwelling", 
                          "2112: Open well in own yard/plot", "2120: Open public well", 
                          "2200: Protected well", "2211: Protected well in own dwelling",
                          "2212: Protected well in own yard/plot", "2220: Protected public well",
                          "2230: Tube well or borehole", "3110: Spring", "3120: Unprotected spring/surface water",
                          "3200: River/dam/lake/ponds/streams/canal/irrigation channel", "3210: River, stream",
                          "3220: Pond, lake", "3230: Dam", "4000: Rainwater", "5100: Tanker truck",
                          "5200: Cart with small tank", "5400: Bottled water")
levels(dhs.data$water)<-c(4,4,3,3,0,0,0,0,2,2,2,2,2,2,0,0,0,0,0,2,1,1,1)

dhs.data$roof<-as.factor(ifelse(dhs.data$ROOF==400|dhs.data$ROOF==998,NA, dhs.data$ROOF))
levels(dhs.data$roof)<-c("0: No roof","111: Thatch/palm leaf", "120: Earth",
                         "210: Rustic mat", "230: Palm/bamboo", "240: Wood planks",
                         "310: Metal", "320: Wood", "331: Cement", "332: Concrete",
                         "343: Tiles/slate", "353: Zinc/cement fiber", "360: roofing singles")
levels(dhs.data$roof)<-c(0, 0, 0, 1, 1, 1, 2, 3, 3, 3, 3, 3, 3)

dhs.data$walls<-as.factor(ifelse(dhs.data$WALL==400|dhs.data$WALL==998,NA, dhs.data$WALL))
levels(dhs.data$walls)<-c("0: No walls", "111: Cane/palm/trunks", "125: Dirt",
                          "210: Bamboo with mud", "220: Stone with mud",
                          "230: Rough wood", "232: Reused wood", 
                          "260: Uncovered adobe", "310: Cement/concrete",
                          "320: Bricks", "330: Cement blocks", "340: Wood planks/shingles",
                          "351: Stone with lime/cement", "360: covered adobe")
levels(dhs.data$walls)<-c(0,0,0,0,3,2,2,1,3,3,3,2,3,3)

dhs.data$floor<-as.factor(ifelse(dhs.data$FLOOR==400|dhs.data$FLOOR==998,NA, dhs.data$FLOOR))
levels(dhs.data$floor)<-c("111: Earth, sand", "120: Dung-based", "211: Wood planks",
                          "213: Wood/palm/bamboo", "220: Palm/bamboo", "320: Vinyl/asphalt strips/linoleum",
                          "331: Ceramic tiles", "340: Cement/concrete", "350: carpet", "390: other finished")
levels(dhs.data$floor)<-c(0,0,1,1,1,2,2,2,2,2)

dhs.data2<- dhs.data %>% rename(landline= HHPHONEHH,
                     mobile = MOBPHONE,
                     computer=PC, 
                     bicycle=BIKEHH,
                     car=CARHH,
                     motorcycle_scooter=MOTORCYCLHH,
                     boat=BOATWMOTOR,
                     canoe=CANOE,
                     bed=BED,
                     radio=RADIOHH,
                     animal_cart=DRAWNCART,
                     fridge = FRIDGEHH,
                     stove_eg=STOVE,
                     tv=TVHH,
                     clock=WATCHCLOCK,
                     chairs=CHAIR,
                     axehoe=AXEHOEYN,
                     lamp=LAMP, 
                     sewingmachine=SEWMACHINE,
                     generator=GENERATYN)

ridmiss.dat<-as.data.frame(apply(dhs.data2[,c("landline", "mobile", "computer", 
                                    "bicycle", "car", "motorcycle_scooter", "bed",
                                    "boat", "radio", "tv","canoe", "fridge", "stove_eg", "axehoe",
                                    "clock",  "chairs", "lamp", "sewingmachine", "generator")], MARGIN=2, ridmiss, yIdx="1", mIdx="8"))

dhs.data<-cbind(dhs.data, ridmiss.dat)

dhs.data$vehicle<-ifelse(dhs.data$car==1,2,
                          ifelse(dhs.data$motorcycle_scooter==1,1,
                                 ifelse(dhs.data$bicycle==1,0,NA)))

possess<-c("landline", "mobile", "computer", 
           "bicycle", "car", "motorcycle_scooter", "bed",
           "boat", "fridge", "radio", "tv", "canoe", "stove_eg", "axehoe",
           "clock", "chairs", "lamp", "sewingmachine",  "generator")

dhs.data$possessions<-(rowSums(dhs.data[,c(possess)], na.rm=TRUE))/length(possess)

dhs.data$memroom<-dhs.data$rooms/dhs.data$members

missing.vars<-which(needed%ni%names(dhs.data))
for(i in 1:length(missing.vars)){
  current_names<-names(dhs.data)
  newvar<-rep(NA, nrow(dhs.data))
  dhs.data<-cbind(dhs.data, newvar)
  colnames(dhs.data)<-c(current_names, needed[missing.vars[i]])
}

DHS07<-subset(dhs.data, year==2007)
DHS13<-subset(dhs.data, year==2013)

shapedhs07<-readOGR("Data/Exposure_data/DRC/inputs/DHS_2007/CDGE52FL", layer="CDGE52FL") 
shapedhs13<-readOGR("Data/Exposure_data/DRC/inputs/DHS_2013/CDGE61FL", layer="CDGE61FL") 

dhs07<-merge(DHS07, shapedhs07,  by="DHSID")
dhs13<-merge(DHS13, shapedhs13, by="DHSID")


DHS07.coords<-data.frame(coords.x1=dhs07$coords.x1, coords.x2=dhs07$coords.x2)
DHS13.coords<-data.frame(coords.x1=dhs13$coords.x1, coords.x2=dhs13$coords.x2)

###########
# MICS 2010
###########
mics4<-read.spss("Data/Exposure_data/DRC/inputs/MICS_2010_geo/SPSS_dataset/hh.sav", to.data.frame=TRUE)
mics4<-dplyr::rename(mics4, cluster_num= HH1, hh_num= HH2, members=HH11, electricity=HC8A, 
                     radio=HC8B, tv=HC8C, landline=HC8D, fridge=HC8E, lamp1=HC8F, bed=HC8G,
                     lamp2=HC8H, watch=HC9A, mobile=HC9B, bicycle=HC9C, motorcycle_scooter=HC9D,
                     animal_cart=HC9E, car_truck=HC9F, boat=HC9G, canoe=HC9H,computer=HC9I,
                     house=HC10,land=HC11, landarea=HC12, fuel=HC6, 
                     water=WS1, toilet=WS8, toilet_share=WS9, rooms=HC2, roof=HC4, 
                     floor=HC3, walls=HC5, generator=HC8F, axehoe=HC8I)

mics4$urban<-ifelse(mics4$HH6=="Rural", 0,1)
mics4$lamp<-max(mics4$lamp1, mics4$lamp2, na.rm=T)

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

mics4<-merge(mics4, mics4geo[,c("cluster_num", "lat", "lon", "district","issue")], by="cluster_num")

mics4$landarea<-ifelse(mics4$landarea=="Plus de 95", 95,
                       ifelse(mics4$landarea=="NSP"|mics4$landarea=="Manquant",NA, as.numeric(mics4$landarea)))
mics4$landarea<-as.numeric(as.character(mics4$landarea))*2.47105

mics4$land<-ifelse(mics4$land=="Oui",1, 
                   ifelse(mics4$land=="Non",0,NA))
mics4$landarea<-ifelse(mics4$land==1,mics4$landarea, 
                       ifelse(mics4$land==0, 0,NA))

mics4$survey=rep("MICS2010", nrow(mics4))
mics4$year<-rep(2010, nrow(mics4))

levels(mics4$toilet)<-c(3,3,3,3,3,2,2,1,2,2,1,1,0,NA,NA)
mics4$toilet_share<-ifelse(mics4$toilet_share=="Oui",1,
                           ifelse(mics4$toilet_share=="Non", 0, NA))

levels(mics4$water)<-c(4,4,3,3,2,2,0,2,0,2,1,1,0,1,1,NA,NA)

levels(mics4$fuel)<-c(4,3,3,3,3,2,2,1,0,0,0,NA,NA,NA)

levels(mics4$roof)<-c(0,0,1,1,1,2,1,2,2,3,3,3,NA,NA)

levels(mics4$walls)<-c(0,0,0,0,3,1,2,0,2,3,3,3,3,3,2,NA,NA)

levels(mics4$floor)<-c(0,0,1,1,2,2,2,2,2,NA,NA)

mics4$electricity<-ifelse(mics4$electricity=="Oui",1, 
                          ifelse(mics4$electricity=="Non",0,NA))

mics4$house<-ifelse(mics4$house=="Propriétaire",1,
                    ifelse(mics4$house=="Manquant", NA,0))

goods<-mics4[,c("hh_num", "radio", "tv", "landline", "fridge", 
                "lamp", "bed", "axehoe",
                "computer", "watch", "mobile", "bicycle", "motorcycle_scooter",
                "animal_cart", "car_truck", "boat", "canoe")]

mics4<-mics4[,c("hh_num", "cluster_num", "members", "electricity", "floor", "fuel", "walls", "roof", "water", "toilet", "toilet_share", "rooms", "house", "land", "landarea", "district", "urban", "lat", "lon",  "survey","year")]

ridmiss.dat<-as.data.frame(apply(goods[,c("radio", "tv", "landline", "fridge", 
                                          "lamp", "bed",
                                          "computer", "watch", "mobile", "bicycle", "motorcycle_scooter",
                                          "animal_cart", "car_truck", "boat", "canoe", "axehoe")], MARGIN=2, ridmiss, yIdx="Oui", mIdx="Manquant"))

mics4<-cbind(mics4, ridmiss.dat)

mics4$vehicle<-ifelse(mics4$car_truck==1,2,
                      ifelse(mics4$motorcycle_scooter==1,1,
                             ifelse(mics4$bicycle==1|mics4$animal_cart==1,0,NA)))

possess<-c("radio", "tv", "landline", "fridge", 
           "lamp", "bed",
           "computer", "watch", "mobile", "bicycle", "motorcycle_scooter",
           "animal_cart", "car_truck", "boat", "canoe", "axehoe")

mics4$possessions<-(rowSums(mics4[,c(possess)], na.rm=TRUE))/length(possess)

mics4$rooms<-ifelse(mics4$rooms=="Manquant", NA, as.numeric(as.character(mics4$rooms)))

mics4$memroom<-as.numeric(as.character(mics4$rooms))/as.numeric(as.character(mics4$members))

keepvars<-names(mics4)[names(mics4)%ni%c("lamp", "members", "rooms")]

codeline=paste0(keepvars[1]," =head(",keepvars[1],", 1)")
for(i in 2:length(keepvars)){
  codeline_new=paste0(keepvars[i]," =head(",keepvars[i],", 1)")
  codeline=paste(codeline,codeline_new, sep=", ")
}

mics4$HHID<-paste(mics4$cluster_num, mics4$hh_num, sep="_")

mics4$clock<-mics4$watch

missing.vars<-which(needed%ni%names(mics4))
for(i in 1:length(missing.vars)){
  current_names<-names(mics4)
  newvar<-rep(NA, nrow(mics4))
  mics4<-cbind(mics4, newvar)
  colnames(mics4)<-c(current_names, needed[missing.vars[i]])
}

MICS10.coords<-data.frame(coords.x1=mics4$lon, coords.x2=mics4$lat)

MICS10<-mics4

########################
# Bind into one dataset
########################
merged_wealth<-rbind(DHS07[,c(needed)], 
                     DHS13[,c(needed)], MICS10[,c(needed)])

merged_wealth$toilet_share<-ifelse(merged_wealth$toilet_share==1,0,1)

merged_coords<-rbind(DHS07.coords, DHS13.coords, MICS10.coords)
  
merged_wealth$urban<-as.numeric(as.character(merged_wealth$urban))
merged_wealth$toilet<-as.numeric(as.character(merged_wealth$toilet))
merged_wealth$water<-as.numeric(as.character(merged_wealth$water))
merged_wealth$fuel<-as.numeric(as.character(merged_wealth$fuel))
merged_wealth$floor<-as.numeric(as.character(merged_wealth$floor))
merged_wealth$walls<-as.numeric(as.character(merged_wealth$walls))
merged_wealth$roof<-as.numeric(as.character(merged_wealth$roof))

#Make sure all variables for factor analysis are numeric or integer
is.fact<-sapply(merged_wealth, is.factor)
summary(is.fact)

is.char<-sapply(merged_wealth, is.character)
summary(is.char)

is.log<-sapply(merged_wealth, is.logical) 
summary(is.log)

is.num<-sapply(merged_wealth, is.numeric)
summary(is.num)

merged_wealth$hhid<-1:nrow(merged_wealth)
merged_coords$hhid<-1:nrow(merged_coords)

save(merged_wealth, file="Data/Confounder_data/DRC/merged_wealth.csv")

#################
#Factor analysis 
#################
cluster<-merged_wealth$cluster
district<-merged_wealth$district

remove<-which(names(merged_wealth)%in%c("cluster", "district"))
merged_wealth<-merged_wealth[,-c(remove)]

#(1) Prepare data
#(1a) Subset by period
wealth1<-merged_wealth[which(merged_wealth$year<=2007),]
wealth2<-merged_wealth[which(merged_wealth$year>2007&merged_wealth$year<=2012),]
wealth3<-merged_wealth[which(merged_wealth$year>2012),]

#(1b) Remove variables with high misssingness
datasets<-list(wealth1=wealth1, wealth2=wealth2, wealth3=wealth3)
newdat<-list()
for(i in 1:length(datasets)){
  data<-datasets[[i]]
  miss<-propmiss(data)
  lowmiss<-miss[which(miss$propmiss<=0.05),] 
  names.old<-names(data)
  names.new<-lowmiss$variable
  index<-which(names.old %in% names.new)
  dat<-data[,index]
  newdat<-c(newdat, list(dat))
}

#(1c) Remove variables with 0 variance
newdat_sd<-list()
for(i in 1:length(newdat)){
  data<-newdat[[i]]
  is.nonzeroSD<-sapply(data, nonzeroSD)
  scratch<-as.data.frame(is.nonzeroSD)
  scratch$nonzeroSD<-ifelse(scratch$is.nonzeroSD==TRUE,1,0)
  scratch2<-scratch[which(scratch$nonzeroSD==1),]
  names.old<-names(data)
  names.new<-rownames(scratch2) 
  index<-which(names.old %in% names.new)
  new_dat<-data[,index]
  new_dat$survey<-data$survey
  newdat_sd<-c(newdat_sd, list(new_dat))
}

#(1e) Discard variables not associated with wealth
discard<-c("hhid", "survey", "year")
newdat_keep<-list()
keepvars<-list()
for(i in 1:length(newdat_sd)){
  data<-newdat_sd[[i]]
  names<-names(data)
  keep<-which(names %ni% discard)
  keepvars<-c(keepvars, keep)
  data_fa<-data[,keep]
  newdat_keep<-c(newdat_keep, list(data_fa))
}

#(2) Create correlation matrix
matrices<-list()
for(i in 1:length(newdat_keep)){
  data=newdat_keep[[i]]
  matrix<-mixedCor(data= data, ncat=2)#Variables with >2 categories = continuous
  matrices<-c(matrices, list(matrix))
}

corrmat_warnings<-warnings()

#(3) Do factor analysis (from documentation: principal components extraction using correlation method with one factor
#extracted, substitution of mean for missing values, estimation of the factor scores using the regression
#method)
fits<-list()
scores<-list()
for(i in 1:length(matrices)){
  matrix<-matrices[[i]]
  data<-newdat_keep[[i]]
  fit <- fa(matrix$rho, 1, n.obs=nrow(data), fm="pa", scores="regression", missing="TRUE", weight="NULL") 
  fits<-c(fits, list(fit))
  score<-as.data.frame(factor.scores(data,fit, impute="mean")$scores)$PA1
  scores<-c(scores, list(score))
}
fa_warnings<-warnings()

wealth1$comb.score<-scores[[1]]
wealth2$comb.score<-scores[[2]]
wealth3$comb.score<-scores[[3]]

wealthscoredat<-rbind(wealth1[,c("hhid", "comb.score")], wealth2[,c("hhid", "comb.score")],wealth3[,c("hhid", "comb.score")])
merged_wealth<-merge(merged_wealth, wealthscoredat, by="hhid")

merged_wealth$district<-district
merged_wealth$cluster<-cluster

#####################
#Check wealth scores
#####################
discard<-c("hhid", "survey", "year", "cluster")
names<-names(merged_wealth)
keep<-which(names %ni% discard)

wealth_check<-merged_wealth[,keep]

for(n in 1:ncol(wealth_check)){
  variable<-names(wealth_check)[n]
  idx<-which(names(wealth_check)==variable)
  miss<-summary(is.na(wealth_check[,idx]))
  if(miss[2]<nrow(wealth_check)){
    model<-lm(wealth_check$comb.score~wealth_check[,idx])
    res<-model$coefficients[2]
    if(is.na(res)==FALSE & res<0){
      miss<-(as.numeric((summary(is.na(wealth_check[,idx]))[3]))/nrow(wealth_check))*100
      message(paste0(variable, ":(miss ",round(miss),"%)"))
      print(res)
      
    }
  }else{
    res="100% missing"
    print(paste(variable, res))
  }
}

#####################
#Collapse on cluster
#####################
wealth<-merge(merged_wealth, merged_coords, by="hhid")
wealth$location<-paste(wealth$coords.x1, wealth$coords.x2)
wealth$cluster<-ifelse(is.na(wealth$cluster), wealth$location, wealth$cluster)
wealth.cluster<-wealth%>%group_by(cluster, survey)%>%
  summarise(coords.x1= head (coords.x1, 1), coords.x2 = head (coords.x2,1), 
            comb.score = mean(comb.score, na.rm=T), 
            year=head(year,1), 
            urban=head(urban, 1), district=head(district,1))%>% 
  arrange(cluster, survey)
wealth.cluster<-as.data.frame(wealth.cluster)
wealth.cluster$cluster<-1:nrow(wealth.cluster)

wealth.cluster<-wealth.cluster[!is.na(wealth.cluster$coords.x1), ]

wealth.coords<-data.frame(coords.x1=wealth.cluster$coords.x1, coords.x2=wealth.cluster$coords.x2)

wealth.sp <- SpatialPointsDataFrame(coords = wealth.coords, data = wealth.cluster,
                                   proj4string = CRS(proj))

ext<-polygon_from_extent(districts)
  
wealth.sp <- crop(wealth.sp, ext)

wealth.sp@data$nightlights<-rep(NA, nrow(wealth.sp@data))

#Add nightlights
years<-unique(wealth.sp@data$year)
for(year in years){
  if(year<=2013){
    nightYear=year
  }else{
    nightYear=2013
  }
  nightRaster = raster(paste('Data/Predictor_data/Nighttime_lights/nightlights_',nightYear,'.tif', sep=""))
  subset<-wealth.sp[which(wealth.sp@data$year==year),]
  nightlights.crop<-crop(nightRaster, ext)
  merged.trans <- spTransform(subset, crs(nightlights.crop))
  nightlights.extract<-raster::extract(nightlights.crop, merged.trans)
  idx<-which(wealth.sp@data$cluster%in%subset@data$cluster)
  wealth.sp@data$nightlights[idx]<-nightlights.extract
  wealth.sp@data$nightlights<-ifelse(wealth.sp@data$nightlights==255, NA, wealth.sp@data$nightlights)
}

#Take care of areas with missing district data
merge1<-over(wealth.sp, WCA1)
merge1$province<-merge1$NAME
merge3<-over(wealth.sp, WCA3)
merge3$district<-merge3$NAME
wealth.sp@data$province<-merge1$province
wealth.sp@data$district<-merge3$district

idxNA = which (is.na(wealth.sp@data$district))

missing_district = wealth.sp[idxNA,]
plot(districts)
plot(missing_district, add=TRUE, col="red")

dists <- dist2Line(missing_district,WCA3)
missing_district@data$district<-as.factor(WCA3@data$NAME[dists[,4]])

wealth.sp@data$district<-ifelse(is.na(wealth.sp@data$district), as.character(missing_district@data$district), as.character(wealth.sp@data$district))

wealth.sp@data$district<-as.factor(wealth.sp@data$district)
wealth.sp@data$province<-as.factor(wealth.sp@data$province)
wealth.sp@data$prop_u2<-rep(NA, nrow(wealth.sp@data))#WorldPop doesn't do this for DRC

#---------------------------------------------#
# Fill in NAs for missing years for prediction
#---------------------------------------------#
all.years<-2000:2020
missing.years<-all.years[which(all.years %ni% wealth.sp@data$year)]
for (year in missing.years){
  wealth.sp@data = rbind(wealth.sp@data[1,], wealth.sp@data)
  wealth.sp@data$comb.score[1]= wealth.sp@data$survey[1]= 
    wealth.sp@data$urban[1] = wealth.sp@data$district[1]=
    wealth.sp@data$cluster[1] = wealth.sp@data$nightlights[1] =
    wealth.sp@data$prop_u2[1] = NA
  wealth.sp@data$year[1] = year
}

min.year<-min(wealth.sp@data$year)-1
wealth.sp@data$time<-wealth.sp@data$time.unstruct<-wealth.sp@data$year-min.year

rownames(wealth.sp@data)<-1:nrow(wealth.sp@data)

wealthdata<-wealth.sp@data
wealthcoords<-data.frame(coords.x1=wealthdata$coords.x1, coords.x2=wealthdata$coords.x2)

wealth <- SpatialPointsDataFrame(coords = wealthcoords, data = wealthdata,
                                    proj4string = CRS(proj))

writeOGR(obj=wealth, dsn="Data/Confounder_data/DRC/wealth", layer="wealth", driver="ESRI Shapefile", overwrite_layer = TRUE)