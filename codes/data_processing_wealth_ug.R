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
UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
districts<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_4")
proj<-proj4string(districts)

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
          "roof", "electricity", "bank", "landarea", "chairs", "table", "vehicle", "clock", 
          "fridge", "mobile", "landline", "tv", "computer", "radio", 
          "bed", "sofa", "boat",  "solar_panel", "generator", "possessions", "cassette", "cupboard", "district")

################################################
# DHS
################################################
dhs.data <- read.csv("Data/Exposure_data/Uganda/inputs/DHS/idhs_00016.csv")
dhs.data<-subset(dhs.data, COUNTRY==800)
dhs.data$district_06<-as.factor(dhs.data$GEO_UG2006_2016)
dhs.data$district_01<-as.factor(dhs.data$GEO_UG2001)

dhs.data$hhid<-paste(dhs.data$DHSID, dhs.data$HHNUMALL)
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
levels(dhs.data$fuel)<-c("Electricity","LPG_natural gas","LPG", "Kerosene", "Biogas",
                     "Wood", "Charcoal", "Firewood_straw", "Straw_shrub_grass",
                    "Dung","Agricultural_crop_based")
levels(dhs.data$fuel)<-c(4,3,3,3,3,1,2,1,0,0,0)

dhs.data$toilet<-as.factor(ifelse(dhs.data$TOILETTYPE==5000|dhs.data$TOILETTYPE==9998,NA, dhs.data$TOILETTYPE))
levels(dhs.data$toilet)<-c("NO FACILITY","Unspecified type of flush toilet",
                      "Own flush toilet (unspecified type)","Shared flush toilet (unspecified type)",
                      "Flush to piped sewer system", "Flush to septic tank","Flush to pit latrine","Flush to somewhere else",
                      "Flush, dont know where","Composting toilet", "Ecosan toilet", 
                      "Traditional pit toilet or latrine",
                      "Pit latrine without slab or open pit","Pit latrine with slab","Ventilated improved pit latrine", "Covered pit latrine, no slab", "Covered pit latrine, with slab",
                      "Bucket toilet","Hanging latrine over water source")
levels(dhs.data$toilet)<-c(0,3,3,3,3,3,3,3,3,2,2,1,1,2,2,2,2,1,1)

dhs.data$water<-as.factor(ifelse(dhs.data$DRINKWTR==6000|dhs.data$DRINKWTR==2310|dhs.data$DRINKWTR==2320|dhs.data$DRINKWTR==9998,NA, dhs.data$DRINKWTR))
levels(dhs.data$water)<-c("Piped into own dwelling/yard/plot","Piped into own dwelling","Piped into own yard/plot", 
                          "Public tap/standpipe","Piped into neighbors dwelling/yard", "Unprotected/open well",
                     "Open well in own dwelling/yard/plot","Open public well",
                     "Protected well","Protected well in own dwelling/yard/plot","Protected well", "Tube well or borehole",
                      "Borehole in yard/plot",
                     "Public borehole", "Spring", "Protected spring", 
                     "Unprotected spring", "River/dam/lake/ponds/streams/canal/irrigation channel", "River/dam/lake/ponds/streams/canal/irrigation channel",
                     "Pond, lake", "Dam", "Channeled by gravity flow scheme", "Rainwater", "Tanker truck", "Cart with small tank", "Water vendor from unknown source", "Bottled water", "Sachet water")
levels(dhs.data$water)<-c(3,3,3,2,2,0,0,0,2,2,2,2,2,0,2,0,0,0,0,2,2,1,1,1,1,1)

dhs.data$roof<-as.factor(ifelse(dhs.data$ROOF==400|dhs.data$ROOF==998,NA, dhs.data$ROOF))
levels(dhs.data$roof)<-c("No roof","Thatch/palm leaf",
                    "Grass, thatch","Mud", "Rustic mat", "Plastic/olythene sheet", 
                    "Wood planks", "Cardboard", "Tin cans", "Iron sheets", "Corrugated iron", 
                    "Tin", "Absetos", "Wood", "Cement", "Concrete", "Tiles", "Roofing shingles")
levels(dhs.data$roof)<-c(0,0,0,0,1,1,1,1,1,2,2,2,3,3,3,3,3,3)

dhs.data$walls<-as.factor(ifelse(dhs.data$WALL==400|dhs.data$WALL==998,NA, dhs.data$WALL))
levels(dhs.data$walls)<-c("No walls","Dirt", "Cane/palm/trunks",
                          "thatch or straw", "bamboo with mud", "stone with mud", 
                          "plywood", "reused wood", "timber", "poles and mud", 
                          "cardboard", "unburnt bricks", "unburnt brick and plaster", "unburnt brick with cement",
                          "cement/concrete", "burnt bricks with mud", "burnt bricks with cement", 
                          "cement blocks", "wood planks/shingles", "stone", "stone with lime/cement")
levels(dhs.data$walls)<-c(0,0,0,0,0,3,2,2,0,0,0,1,1,1,3,3,3,3,2,3,3)

dhs.data$floor<-as.factor(ifelse(dhs.data$FLOOR==400|dhs.data$FLOOR==998,NA, dhs.data$FLOOR))
levels(dhs.data$floor)<-c("Earth, sand","Dung","Earth and dung", "Wood planks","Palm/bamboo",
                     "Parquet/polished wood","Vinyl/asphalt strips/linoleum",
                     "Tiles/mosaic","Ceramic tiles","Cement/concrete","Carpet","Stone", "Bricks")
levels(dhs.data$floor)<-c(0,0,0,1,1,2,2,2,2,2,2,2,2)

dhs.data$TABLECHARYN<-ifelse(dhs.data$TABLE==1|dhs.data$CHAIR==1,1,0)

dhs.data2<- dhs.data %>% rename(landline= HHPHONEHH,
                     mobile = MOBPHONE,
                     computer=PC, 
                     bicycle=BIKEHH,
                     car=CARHH,
                     motorcycle_scooter=MOTORCYCLHH,
                     boat=BOATWMOTOR,
                     boatnomotor=BOATNOMOTOR,
                     boat01=BOATANY, 
                     fridge=FRIDGEHH,
                     bed=BED,
                     radio=RADIOHH,
                     animal_cart=DRAWNCART,
                     tv=TVHH,
                     watch=WATCHCLOCK,
                     clock=CLOCKONLY,
                     sofa=SOFA,
                     chairs=CHAIR,
                     table=TABLE, 
                     lamp=LAMP, 
                     sewingmachine=SEWMACHINE,
                     tape_cd_dvd_hifi=DVDVCD, 
                     cassette=CASSETTEYN,
                     generator=GENERATYN, 
                     cupboard=CUPBOARDYN, 
                     axehoe=AXEHOEYN)

dhs.data2$boat<-ifelse(dhs.data2$boat==1|dhs.data2$boatnomotor==1|dhs.data2$boat01==1,1,0)
dhs.data2$clock<-ifelse(dhs.data2$clock==1|dhs.data2$watch==1,1,0)

ridmiss.dat<-as.data.frame(apply(dhs.data2[,c("landline", "mobile", "computer", 
                                    "bicycle", "car", "motorcycle_scooter", "bed",
                                    "boat", "fridge", "radio", "tv", "cupboard", "axehoe", "cassette", 
                                    "clock", "sofa", "table", "chairs", "lamp", "sewingmachine", "tape_cd_dvd_hifi", "generator")], MARGIN=2, ridmiss, yIdx="1", mIdx="8"))

dhs.data<-cbind(dhs.data, ridmiss.dat)

dhs.data$vehicle<-ifelse(dhs.data$car==1,2,
                          ifelse(dhs.data$motorcycle_scooter==1,1,
                                 ifelse(dhs.data$bicycle==1,0,NA)))

possess<-c("landline", "mobile", "computer", 
           "bicycle", "car", "motorcycle_scooter", "bed",
           "boat", "fridge", "radio", "tv","cupboard", "axehoe", "cassette", 
           "clock", "sofa", "table", "chairs", "lamp", "sewingmachine", "tape_cd_dvd_hifi", "generator")

dhs.data$possessions<-(rowSums(dhs.data[,c(possess)], na.rm=TRUE))/length(possess)

dhs.data$memroom<-dhs.data$rooms/dhs.data$members

missing.vars<-which(needed%ni%names(dhs.data))
for(i in 1:length(missing.vars)){
  current_names<-names(dhs.data)
  newvar<-rep(NA, nrow(dhs.data))
  dhs.data<-cbind(dhs.data, newvar)
  colnames(dhs.data)<-c(current_names, needed[missing.vars[i]])
}

DHS01<-subset(dhs.data, year==2001)
DHS06<-subset(dhs.data, year==2006)
DHS11<-subset(dhs.data, year==2011)
DHS16<-subset(dhs.data, year==2016)

mshapedhs01<-readOGR("Data/Exposure_data/Uganda/inputs/DHS/shapefiles/UGGE43FL_00", layer="UGGE43FL") 
mshapedhs06<-readOGR("Data/Exposure_data/Uganda/inputs/DHS/shapefiles/UGGE53FL_06", layer="UGGE53FL") 
mshapedhs11<-readOGR("Data/Exposure_data/Uganda/inputs/DHS/shapefiles/UGGE61FL_11", layer="UGGE61FL") 
mshapedhs16<-readOGR("Data/Exposure_data/Uganda/inputs/DHS/shapefiles/UGGE7AFL_16", layer="UGGE7AFL") 

dhs01<-merge(DHS01, mshapedhs01,  by="DHSID")
dhs06<-merge(DHS06, mshapedhs06, by="DHSID")
dhs11<-merge(DHS11, mshapedhs11,  by="DHSID")
dhs16<-merge(DHS16, mshapedhs16,  by="DHSID")

DHS01.coords<-data.frame(coords.x1=dhs01$coords.x1, coords.x2=dhs01$coords.x2)
DHS06.coords<-data.frame(coords.x1=dhs06$coords.x1, coords.x2=dhs06$coords.x2)
DHS11.coords<-data.frame(coords.x1=dhs11$coords.x1, coords.x2=dhs11$coords.x2)
DHS16.coords<-data.frame(coords.x1=dhs16$coords.x1, coords.x2=dhs16$coords.x2)

###########
# MIS 2009
###########
household09<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2009/UGHR5HDT/UGHR5HFL.DTA", convert.dates = T,convert.factors = F)
individual09<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2009/UGIR5HDT/UGIR5HFL.DTA", convert.dates = T,convert.factors = F)
geographic09<-readOGR("Data/Exposure_data/Uganda/inputs/MIS_2009/UGGE5AFL", layer="UGGE5AFL")

household09<-dplyr::rename(household09, DHSCLUST= hv001, hh_id=hv002, 
                           members=hv012, rooms=hv216, land=hv244, landarea=hv245, 
                           bank=hv247,electricity=hv206, radio=hv207, tv=hv208, watch=hv243b,
                           fridge=hv209, bicycle=hv210, motorcycle_scooter=hv211,
                           car_truck=hv212, floor=hv213, walls=hv214, roof=hv215,
                           landline=hv221, fuel=hv226, toilet_share=hv225,
                           mobile=hv243a, watch=hv243b, animal_cart=hv243c, 
                           boat=hv243d, toilet=hv205, water=hv201, 
                           cassette=sh104c, table=sh104h, chairs=sh104i, 
                           sofa=sh104j, bed=sh104k, cupboard=sh104l, 
                           clock=sh104m, boatnmotor=sh112g)

household09$HHID<-paste(household09$DHSCLUST, household09$hh_id, sep="_")

household09<-merge(household09, geographic09@data[,c(4,13,15)], by = "DHSCLUST")
household09$urban<-ifelse(household09$URBAN_RURA=="R",0,1)
household09$district<-household09$DHSREGNA

levels(household09$district)<-c("Central 1", "Central 2", "East Central", "Kampala", "Mid Eastern", "Mid Northern", "Mid Western", "Northeast", "South", "West Nile")

household09$boat<-ifelse(household09$boat==1|household09$boatnmotor==1,1,0)

household09<-household09[,c("DHSCLUST", "hh_id", 
         "members", "rooms", "land", "landarea", 
         "bank","electricity", "radio", "tv",
         "fridge", "bicycle", "motorcycle_scooter",
         "car_truck", "floor", "walls", "roof",
         "landline", "fuel", "toilet_share", 
         "mobile", "clock", "animal_cart", "cassette", "table", 
         "chairs", "sofa", "bed", "cupboard", "watch",
         "boat", "toilet", "water", "district", "urban")]

conversion_hectares<-as.numeric(as.character(household09$landarea))*2.47105
household09$landarea<-ifelse(household09$land==1, conversion_hectares, NA)

household09$survey=rep("MIS2009", nrow(household09))
household09$year<-rep(2009, nrow(household09))

household09$toilet<-as.factor(ifelse(household09$toilet==96|household09$toilet==99, NA, household09$toilet))
levels(household09$toilet)<-c("Flush toilet", "VIP", "Covered pit latrine no slab", "Covered pit latrine with slab", "Uncovered pit latrine no slab", 
                              "Uncovered pit latrine with slab", "None", "Composting toilet")
levels(household09$toilet)<-c(3,2,2,2,1,2,0,2)

household09$water<-as.factor(ifelse(household09$water==96|household09$water==99, NA, household09$water))
levels(household09$water)<-c("11: Piped into dwelling", "12: Piped to yard/plot", "13: Public tap/standpipe", 
                             "22: Open well in yard/compound", "23: Open public well", "33: Protected well in yard/compound", 
                             "34: Protected public well", "35: Borehole", "41: Protected spring", "42: Unprotected spring", 
                             "44: River/stream", "45: Pond/lake", "46: Dam", "51: Rainwater", "61: Tanker truck", "71: Bottled water")
levels(household09$water)<-c(4,4,3,0,0,2,2,2,2,0,0,0,0,2,1,1)

household09$fuel<-as.factor(ifelse(household09$fuel==96|household09$fuel==95|household09$fuel==99, NA, household09$fuel))
levels(household09$fuel)<-c("Electricity","LPG/ natural gas","Biogas", "Kerosene", "Charcoal", "Wood", "Straw/shrubs/grass")
levels(household09$fuel)<-c(4,3,3,3,2,1,0)

household09$roof<-as.factor(ifelse(household09$roof==96|household09$roof==99, NA, household09$roof))
levels(household09$roof)<-c("11: Thatched", "12: Mud", "21: Wood/planks", "22: Iron sheets", "23: Asbestos", "24: Tiles", "25: Tin", "26: Cement")
levels(household09$roof)<-c(0,0,1,2,3,3,2,3)

household09$walls<-as.factor(ifelse(household09$walls==96|household09$walls==99, NA, household09$walls))
levels(household09$walls)<-c("11: Thatched/straw", "21: Mud and poles", "22: Un-burnt bricks", "23: Unburnt bricks with plaster", "24: Burnt bricks with mud", "31: Cement bricks", "32: Stone", "33: Timber", "34: Burnt bricks with cement")
levels(household09$walls)<-c(0,0,1,1,3,3,3,2,3)

household09$floor<-as.factor(ifelse(household09$floor==96|household09$floor==99, NA, household09$floor))
levels(household09$floor)<-c("11: Earth/sand", "12: Earth and dung", "33: Mosaic or tiles", "34: Bricks", "35: Cement", "36: Stones")
levels(household09$floor)<-c(0,0,2,2,2,2)

household09$rooms<-ifelse(household09$rooms==0,1,household09$rooms)

household09$vehicle<-ifelse(household09$car_truck==1,2,
                      ifelse(household09$motorcycle_scooter==1,1,
                             ifelse(household09$bicycle==1|household09$animal_cart==1,0,NA)))

possess<-c("radio", "tv",
           "fridge", "bicycle", "motorcycle_scooter",
           "car_truck", "landline", "mobile", "clock", "animal_cart",  
           "cassette", "table", 
           "chairs", "sofa", "bed", "cupboard", "watch",
           "boat")

household09$possessions<-(rowSums(household09[,c(possess)], na.rm=TRUE))/length(possess)

household09$memroom<-as.numeric(as.character(household09$rooms))/as.numeric(as.character(household09$members))

missing.vars<-which(needed%ni%names(household09))
for(i in 1:length(missing.vars)){
  current_names<-names(household09)
  newvar<-rep(NA, nrow(household09))
  household09<-cbind(household09, newvar)
  colnames(household09)<-c(current_names, needed[missing.vars[i]])
}

household<-merge(household09, geographic09[,c("DHSCLUST", "LATNUM", "LONGNUM")], by="DHSCLUST")

MIS09.coords<-data.frame(coords.x1=household$LONGNUM, coords.x2=household$LATNUM)

MIS09<-household

###########
# MIS 2014
###########
household14<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2014/UGHR72DT/UGHR72FL.DTA", convert.dates = T,convert.factors = F)
individual14<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2014/UGIR72DT/UGIR72FL.DTA", convert.dates = T,convert.factors = F)
geographic14<-readOGR("Data/Exposure_data/Uganda/inputs/MIS_2014/UGGE71FL", layer="UGGE71FL")

household14<-dplyr::rename(household14, DHSCLUST= hv001, hh_id=hv002, 
                           members=hv012, rooms=hv216, land=hv244, landarea=hv245, 
                           bank=hv247,electricity=hv206, radio=hv207, tv=hv208,
                           fridge=hv209, bicycle=hv210, motorcycle_scooter=hv211,
                           car_truck=hv212, floor=hv213, walls=hv214, roof=hv215,
                           landline=hv221, fuel=hv226, toilet_share=hv225,
                           mobile=hv243a, watch=hv243b, animal_cart=hv243c, 
                           boat=hv243d, toilet=hv205, water=hv201, 
                           cassette=sh107c, table=sh107h, chairs=sh107i, 
                           sofa=sh107j, bed=sh107k, cupboard=sh107l, 
                           clock=sh107m, boatnmotor=sh113g)

household14$HHID<-paste(household14$DHSCLUST, household14$hh_id, sep="_")

household14<-merge(household14, geographic14@data[,c(4,13,15)], by = "DHSCLUST")
household14$urban<-ifelse(household14$URBAN_RURA=="R",0,1)
household14$district<-household14$DHSREGNA

levels(household14$district)<-c("Central 1", "Central 2", "East Central", "Kampala", "Mid Eastern", "Mid Northern", "Mid Western", "Northeast", "South Western", "West Nile")
household14$boat<-ifelse(household14$boat==1|household14$boatnmotor==1,1,0)

household14<-household14[,c("DHSCLUST", "hh_id", 
                            "members", "rooms", "land", "landarea", 
                            "bank","electricity", "radio", "tv",
                            "fridge", "bicycle", "motorcycle_scooter",
                            "car_truck", "floor", "walls", "roof",
                            "landline", "fuel", "toilet_share",
                            "mobile", "clock", "animal_cart", "watch", "cassette", "table", 
                            "chairs", "sofa", "bed", "cupboard",
                            "boat", "toilet", "water", "district", "urban")]

conversion_hectares<-as.numeric(as.character(household14$landarea))*2.47105
household14$landarea<-ifelse(household14$land==1, conversion_hectares, NA)

household14$survey=rep("MIS2014", nrow(household14))
household14$year<-rep(2014, nrow(household14))

household14$toilet<-as.factor(ifelse(household14$toilet==96, NA, household14$toilet))
levels(household14$toilet)<-c("11: Flush to piped sewer system", "12: Flush to septic tank", "13: Flush to pit latrine",
                              "14: Flush to somewhere else", "21: VIP", "22: Covered pit latrine with slab", "23: Covered pit latrine without slab", "24: Uncovered pit latrine with slab", 
                              "25: Uncovered pit latrine without slab", "31: No facility/bush/field", "41: Composting toilet", "42: Bucket toilet",
                              "43: Hanging latrine/toilet")
levels(household14$toilet)<-c(3,3,3,3,2,2,2,2,1,0,2,1,1)

household14$water<-as.factor(ifelse(household14$water==96, NA, household14$water))
levels(household14$water)<-c("11: Piped into dwelling", "12: Piped to yard/plot", "13: Public tap/standpipe", 
                             "21: Borehole in yard/plot", "22: Public borehole", "31: Protected well", "32: Unprotected well",
                             "41: Protected spring", "42: Unprotectd spring", "43: River/dam/lake/ponds/stream/canal/irrigation channel", 
                             "44: Gravity flow scheme", "51: Rainwater", "61: Tanker truck", "62: Vendor", "63: Cart with small tank", "71: Bottled water")
levels(household14$water)<-c(4,4, 3, 2, 2, 2, 0, 2, 0, 0, 2, 2, 1,1,1,1)

household14$fuel<-as.factor(ifelse(household14$fuel==95|household14$fuel==96, NA, household14$fuel))
levels(household14$fuel)<-c("1: Electricity","2: LPG","3: Natural gas", "4: Biogas", "5: Kerosene", "6: Coal/lignite",
                            "7: Charcoal", "8: Wood", "9: Straw/shrubs/grass", "10: Agricultural crop")
levels(household14$fuel)<-c(4,3,3,3,3,2,2,1,0,0)

household14$roof<-as.factor(ifelse(household14$roof==96, NA, household14$roof))
levels(household14$roof)<-c("12: Thatched", "13: Mud", "21: Tin", "22: Palm", "23: Wood planks", 
                            "31: Iron sheets", "32: Wood", "33: Cement fiber", "34: Tiles", 
                            "35: Cement", "36: Roofing shingles", "37: Asbestos")
levels(household14$roof)<-c(0,0,2,0,1,2,3,3,3,3,3,3)

household14$walls<-as.factor(ifelse(household14$walls==96, NA, household14$walls))
levels(household14$walls)<-c("11: No walls","12: Thatched/straw", "13: Dirt", "21: Mud and poles", 
                             "22: Stone with mud", "23: Reused wood", "24: Unburnt bricks", 
                             "25: Unburnt bricks with plaster", "26: Unburnt bricks with mud", "31: Cement", 
                             "32: Stone with lime/cement", "33: Burnt bricks with cement", "34: Cement blocks",
                             "35: Wood planks/shingles")
levels(household14$walls)<-c(0,0,0,1,3,2,1,1,1,3,3,3,3,2)

household14$floor<-as.factor(ifelse(household14$floor==96, NA, household14$floor))
levels(household14$floor)<-c("11: Earth, sand","12: sand and dung","21: Wood planks", "31: Parquet/polished wood", 
                             "32: Mosaic or tile", "33: Cement", "34: Stones", "35: Bricks")
levels(household14$floor)<-c(0,0,1,2,2,2,2,2)

household14$rooms<-ifelse(household14$rooms==0,1,household14$rooms)

household14$vehicle<-ifelse(household14$car_truck==1,2,
                            ifelse(household14$motorcycle_scooter==1,1,
                                   ifelse(household14$bicycle==1|household14$animal_cart==1,0,NA)))

possess<-c("radio", "tv",
           "fridge", "bicycle", "motorcycle_scooter",
           "car_truck", 
           "landline", 
           "mobile", "clock", "animal_cart", 
           "boat","cassette", "table", 
           "chairs", "sofa", "bed", "cupboard", "watch",
           "boat")

household14$possessions<-(rowSums(household14[,c(possess)], na.rm=TRUE))/length(possess)

household14$memroom<-as.numeric(as.character(household14$rooms))/as.numeric(as.character(household14$members))

missing.vars<-which(needed%ni%names(household14))
for(i in 1:length(missing.vars)){
  current_names<-names(household14)
  newvar<-rep(NA, nrow(household14))
  household14<-cbind(household14, newvar)
  colnames(household14)<-c(current_names, needed[missing.vars[i]])
}

household<-merge(household14, geographic14[,c("DHSCLUST", "LATNUM", "LONGNUM")], by="DHSCLUST")

MIS14.coords<-data.frame(coords.x1=household$LONGNUM, coords.x2=household$LATNUM)

MIS14<-household

###########
# MIS 2018
###########
household18<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2018/UGHR7IDT/UGHR7IFL.DTA", convert.dates = T,convert.factors = F)
individual18<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2018/UGIR7IDT/UGIR7IFL.DTA", convert.dates = T,convert.factors = F)
geographic18<-readOGR("Data/Exposure_data/Uganda/inputs/MIS_2018/shapefiles/UGGE7IFL", layer="UGGE7IFL")

household18<-dplyr::rename(household18, DHSCLUST= hv001, hh_id=hv002, 
                           members=hv012, rooms=hv216, land=hv244, landarea=hv245, 
                           bank=hv247,electricity=hv206, radio=hv207, tv=hv208,
                           fridge=hv209, bicycle=hv210, motorcycle_scooter=hv211,
                           car_truck=hv212, floor=hv213, walls=hv214, roof=hv215,
                           landline=hv221, fuel=hv226, toilet_share=hv225,
                           mobile=hv243a, watch=hv243b, animal_cart=hv243c, 
                           boat=hv243d, toilet=hv205, water=hv201,
                           table=sh114h, chairs=sh114i, 
                           sofa=sh114j, bed=sh114k, cupboard=sh114l, 
                           clock=sh114m)

household18$HHID<-paste(household18$DHSCLUST, household18$hh_id, sep="_")

household18<-merge(household18, geographic18@data[,c(4,13,15)], by = "DHSCLUST")
household18$urban<-ifelse(household18$URBAN_RURA=="R",0,1)
household18$district<-household18$DHSREGNA

household18<-household18[,c("DHSCLUST", "hh_id", 
                            "members", "rooms", "land", "landarea", 
                            "bank","electricity", "radio", "tv",
                            "fridge", "bicycle", "motorcycle_scooter",
                            "car_truck", "floor", "walls", "roof",
                            "landline", "fuel", "toilet_share",
                            "mobile", "clock", "animal_cart",  "table", "chairs", 
                            "sofa", "bed", "cupboard", 
                            "boat", "toilet", "water", "district", "urban")]

conversion_hectares<-as.numeric(as.character(household18$landarea))*2.47105
household18$landarea<-ifelse(household18$land==1, conversion_hectares, NA)

household18$survey=rep("MIS2018", nrow(household18))
household18$year<-rep(2018, nrow(household18))

household18$toilet<-household18$water<-household18$fuel<-household18$roof<-household18$walls<-household18$floor<-rep(NA, nrow(household18))

household18$rooms<-ifelse(household18$rooms==0,1,household18$rooms)

household18$vehicle<-ifelse(household18$car_truck==1,2,
                            ifelse(household18$motorcycle_scooter==1,1,
                                   ifelse(household18$bicycle==1|household18$animal_cart==1,0,NA)))

possess<-c("radio", "tv",
           "fridge", "bicycle", "motorcycle_scooter",
           "car_truck", 
           "landline", "table", "chairs", 
           "sofa", "bed", "cupboard", 
           "boat", 
           "mobile", "clock", "animal_cart")

household18$possessions<-(rowSums(household18[,c(possess)], na.rm=TRUE))/length(possess)

household18$memroom<-as.numeric(as.character(household18$rooms))/as.numeric(as.character(household18$members))

missing.vars<-which(needed%ni%names(household18))
for(i in 1:length(missing.vars)){
  current_names<-names(household18)
  newvar<-rep(NA, nrow(household18))
  household18<-cbind(household18, newvar)
  colnames(household18)<-c(current_names, needed[missing.vars[i]])
}

household<-merge(household18, geographic14[,c("DHSCLUST", "LATNUM", "LONGNUM")], by="DHSCLUST")

MIS18.coords<-data.frame(coords.x1=household$LONGNUM, coords.x2=household$LATNUM)

MIS18<-household

###########
# AIS 2011
###########
household11<-read.dta("Data/Exposure_data/Uganda/inputs/AIS_2011/UGHR6ADT/UGHR6AFL.DTA", convert.dates = T,convert.factors = F)
individual11<-read.dta("Data/Exposure_data/Uganda/inputs/AIS_2011/UGIR6ADT/UGIR6AFL.DTA", convert.dates = T,convert.factors = F)
geographic11<-readOGR("Data/Exposure_data/Uganda/inputs/AIS_2011/shapefiles/UGGE6AFL", layer="UGGE6AFL")

household11<-dplyr::rename(household11, DHSCLUST= hv001, hh_id=hv002, 
                           members=hv012, rooms=hv216, land=hv244, landarea=hv245, 
                           bank=hv247,electricity=hv206, radio=hv207, tv=hv208,
                           fridge=hv209, bicycle=hv210, motorcycle_scooter=hv211,
                           car_truck=hv212, floor=hv213, walls=hv214, roof=hv215,
                           landline=hv221, fuel=hv226, toilet_share=hv225,
                           mobile=hv243a, watch=hv243b, animal_cart=hv243c, 
                           boat=hv243d, toilet=hv205, water=hv201,
                           cassette=sh104c, table=sh104h, chairs=sh104i, 
                           sofa=sh104j, bed=sh104k, cupboard=sh104l, 
                           clock=sh104m, boatnmotor=sh111g)

household11$HHID<-paste(household11$DHSCLUST, household11$hh_id, sep="_")

household11<-merge(household11, geographic11@data[,c(4,13,15)], by = "DHSCLUST")
household11$urban<-ifelse(household11$URBAN_RURA=="R",0,1)
household11$district<-household11$DHSREGNA

levels(household11$district)<-c("Central 1", "Central 2", "East Central", "Kampala", "Mid Eastern", "Mid Northern", "Mid Western", "Northeast", "South Western", "West Nile")
household11$boat<-ifelse(household11$boat==1|household11$boatnmotor==1,1,0)

household11<-household11[,c("DHSCLUST", "hh_id", 
                            "members", "rooms", "land", "landarea", 
                            "bank","electricity", "radio", "tv",
                            "fridge", "bicycle", "motorcycle_scooter",
                            "car_truck", "floor", "walls", "roof",
                            "landline", "fuel", "toilet_share",
                            "mobile", "watch", "animal_cart", 
                            "cassette", "table", "chairs", 
                            "sofa", "bed", "cupboard", 
                            "clock",
                            "boat", "toilet", "water", "district", "urban")]

conversion_hectares<-as.numeric(as.character(household11$landarea))*2.47105
household11$landarea<-ifelse(household11$land==1, conversion_hectares, NA)

household11$survey=rep("AIS2011", nrow(household11))
household11$year<-rep(2011, nrow(household11))

household11$toilet<-as.factor(ifelse(household11$toilet==96, NA, household11$toilet))
levels(household11$toilet)<-c("11: Flush to piped sewer system",  "21: VIP", "22: Pit latrine with slab", "23: Pit latrine without slab", "31: No facility/bush/field", "41: Composting toilet")
levels(household11$toilet)<-c(3,2,2,1,0,2)

household11$water<-as.factor(ifelse(household11$water==96, NA, household11$water))
levels(household11$water)<-c("11: Piped into dwelling", "12: Piped to yard/plot", "13: Public tap/standpipe", 
                             "31: Protected well", "32: Unprotected well",
                             "41: Protected spring", "42: Unprotectd spring", "43: River/dam/lake/ponds/stream/canal/irrigation channel", 
                              "51: Rainwater", "61: Tanker truck",  "71: Bottled water", "81: gravity flow scheme")
levels(household11$water)<-c(4,4, 3, 2, 0,2,0,0,2,1,1,2)

household11$fuel<-as.factor(ifelse(household11$fuel==95|household11$fuel==96, NA, household11$fuel))
levels(household11$fuel)<-c("1: Electricity","3: Natural gas", "5: Kerosene", 
                            "7: Charcoal", "8: Wood", "9: Straw/shrubs/grass", "10: Agricultural crop")
levels(household11$fuel)<-c(4,3,3,2,1,0,0)

household11$roof<-as.factor(ifelse(household11$roof==96, NA, household11$roof))
levels(household11$roof)<-c("11: thatched", "12: Mud", 
                            "31: Wood/planks", "32: Iron sheets", "33: Asbestos", "34: Tiles", 
                            "35: Tin", "36: Cement")
levels(household11$roof)<-c(0,0,2,2,3,3,3,3)

household11$walls<-as.factor(ifelse(household11$walls==96, NA, household11$walls))
levels(household11$walls)<-c("11: Thatched/straw", "21: Mud and poles", "22: Unburnt bricks", "23: Unburnt bricks with plaster", 
                             "24: Burnt bricks iwth mud", "31: Cement blocks", "32: Stone", "33: Timber", "34: Burnt bricks with cement")
levels(household11$walls)<-c(0,1,1,1,3,3,3,2,3)

household11$floor<-as.factor(ifelse(household11$floor==96, NA, household11$floor))
levels(household11$floor)<-c("11: Earth, sand","12: earth and dung","31: Parquet/polished wood", 
                              "33: Mosaic or tiles", "34: Bricks", "35: Cement", "36: Stones")
levels(household11$floor)<-c(0,0,2,2,2,2,2)

household11$rooms<-ifelse(household11$rooms==0,1,household11$rooms)

household11$vehicle<-ifelse(household11$car_truck==1,2,
                            ifelse(household11$motorcycle_scooter==1,1,
                                   ifelse(household11$bicycle==1|household11$animal_cart==1,0,NA)))

possess<-c("radio", "tv",
           "fridge", "bicycle", "motorcycle_scooter",
           "car_truck", 
           "landline", "watch", "cassette", "table", "chairs", 
           "sofa", "bed", "cupboard", 
           "boat", 
           "mobile", "clock", "animal_cart", 
           "boat")

household11$possessions<-(rowSums(household11[,c(possess)], na.rm=TRUE))/length(possess)

household11$memroom<-as.numeric(as.character(household11$rooms))/as.numeric(as.character(household11$members))

missing.vars<-which(needed%ni%names(household11))
for(i in 1:length(missing.vars)){
  current_names<-names(household11)
  newvar<-rep(NA, nrow(household11))
  household11<-cbind(household11, newvar)
  colnames(household11)<-c(current_names, needed[missing.vars[i]])
}

household<-merge(household11, geographic11[,c("DHSCLUST", "LATNUM", "LONGNUM")], by="DHSCLUST")

AIS11.coords<-data.frame(coords.x1=household$LONGNUM, coords.x2=household$LATNUM)

AIS11<-household

#############################
# National Panel Survey 2009
#############################
household09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC2.dta", convert.dates = T,convert.factors = F)
household09b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC1.dta", convert.dates = T,convert.factors = F)
household09c<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC9.dta", convert.dates = T,convert.factors = F)
land09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/AGSEC2A.dta", convert.dates = T,convert.factors = F)
fuel09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC10.dta", convert.dates = T,convert.factors = F)
bank09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC13.dta", convert.dates = T,convert.factors = F)
electricity09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC10A.dta", convert.dates = T,convert.factors = F)
goods09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC14.dta", convert.dates = T,convert.factors = F)
geographic09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/UNPS_Geovars_0910.dta", convert.dates = T,convert.factors = F)

land09$size<-land09$a2aq4
land.hh<-land09 %>% group_by(HHID)%>%
  summarise(landarea=sum(size, na.rm=T))%>% 
  arrange(HHID)%>%as.data.frame()

bank09$bank<-ifelse(bank09$h13q20==1,1,0)

electricity09$electricity<-ifelse(electricity09$h10q1==1,1,0)

household09c<-dplyr::rename(household09c, rooms=h9q03, roof=h9q04,
                            walls=h9q05, floor=h9q06, water=h9q07,
                            toilet=h9q22)

household09c$toilet_share<-ifelse(household09c$toilet==2|household09c$toilet==4|household09c$toilet==7,1,0)
household09c$toilet<-as.factor(ifelse(household09c$toilet==9, NA, household09c$toilet))
levels(household09c$toilet)<-c(2,2,2,2,1,3,3,0)

household09c$water<-as.factor(ifelse(household09c$water==96, NA, household09c$water))
levels(household09c$water)<-c(4,3,2,2,0,0,1,2,2)

household09c$roof<-as.factor(ifelse(household09c$roof==96, NA, household09c$roof))
levels(household09c$roof)<-c(0,0,1,2,3,3,2,3)

household09c$walls<-as.factor(ifelse(household09c$walls==96, NA, household09c$walls))
levels(household09c$walls)<-c(0,0,2,1,3,3,3,3)

household09c$floor<-as.factor(ifelse(household09c$floor==96, NA, household09c$floor))
levels(household09c$floor)<-c(0,0,2,2,2,2,1)

household09b$district<-as.factor(household09b$region)
levels(household09b$district)<-c("Kampala", "Central", "Eastern", "Northern", "Western")
household09b$DHSID<-household09b$comm

household09<-dplyr::rename(household09, relationship=h2q4)
household09$domestic<-ifelse(household09$relationship==9,1,0)
household09<-household09 %>% group_by(HHID)%>%
  summarise(members=n(),
            domestic=max(domestic, na.rm=T))%>% 
  arrange(HHID)%>% as.data.frame()

fuel09$fuel.lev<-as.factor(fuel09$h10q13)
levels(fuel09$fuel.lev)<-c("Firewood", "Dung", "Crop residue",
                           "Kerosene", "LPG", "Charcoal", "Solar", "Electricity")
fuel09$paste<-paste(fuel09$fuel.lev, fuel09$h10q14, sep="_")
fuel09$fuel<-ifelse(fuel09$paste=="Firewood_1", 1,
                    ifelse(fuel09$paste=="Dung_1", 0,
                           ifelse(fuel09$paste=="Crop residue_1", 0,
                                  ifelse(fuel09$paste=="Kerosene_1", 3,
                                         ifelse(fuel09$paste=="LPG_1", 3,
                                                ifelse(fuel09$paste=="Charcoal_1", 2,
                                                       ifelse(fuel09$paste=="Solar_1", 3, 
                                                              ifelse(fuel09$paste=="Electricity_1",4, NA))))))))

goods09$goods<-as.factor(goods09$h14q2)
levels(goods09$goods)<-c("house", "Other buildings", "land", 
                           "table", "fridge", "tv",
                           "radio", "generator", "solar_panel",
                           "bicycle", "motorcycle", "car", "boat", "other_transport",
                           "watch", "mobile", "computer", "internet", "other_electronics",
                           "other_assets", "other_assets", "other_assets")
newvars.g<-levels(goods09$goods)
goods09$own<-ifelse(goods09$h14q3==1,1,0)
goods2<-create_var(newvars.g, goods09, "own", "HHID")
goods2<- goods2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))
hhid.keep<-goods2$HHID
goods2<-as.data.frame(apply(goods2[,newvars.g], MARGIN=2, binary))
goods2$HHID<-hhid.keep

nps09<-merge(land.hh, bank09, by="HHID", all=TRUE)
nps09<-merge(nps09, electricity09, by="HHID", all=TRUE)
nps09<-merge(nps09, household09c, by="HHID", all=TRUE)
nps09<-merge(nps09, household09, by="HHID", all=TRUE)
nps09<-merge(nps09, fuel2, by="HHID", all=TRUE)
nps09<-merge(nps09, goods2, by="HHID", all=TRUE)
nps09<-merge(nps09, household09b, by="HHID", all=TRUE)

keep<-c("HHID", "electricity", "landarea", "bank", "water", "rooms", "roof", "fuel", "toilet","walls", "floor", "water", "toilet_share", "members",
"domestic", "fuel", newvars.g, "district", "DHSID") 
nps09<-nps09[,which(colnames(nps09)%in%keep)]

nps09$vehicle<-ifelse(nps09$car==1,2,
                          ifelse(nps09$motorcycle==1,1,
                                 ifelse(nps09$bicycle==1,0,NA)))

possess<-levels(goods09$goods)

nps09$possessions<-(rowSums(nps09[,c(possess)], na.rm=TRUE))/length(possess)

nps09$survey<-rep("NPS09", nrow(nps09))
nps09$year<-rep(2009, nrow(nps09))

nps09$memroom<-nps09$rooms/nps09$members

NPS09<-merge(nps09, geographic09[c(1,2,3,7)], by = "HHID")

NPS09.coords<-data.frame(coords.x1=NPS09$lon_mod, coords.x2=NPS09$lat_mod)

missing.vars<-which(needed%ni%names(NPS09))
for(i in 1:length(missing.vars)){
  current_names<-names(NPS09)
  newvar<-rep(NA, nrow(NPS09))
  NPS09<-cbind(NPS09, newvar)
  colnames(NPS09)<-c(current_names, needed[missing.vars[i]])
}

#############################
# National Panel Survey 2010
#############################
household10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC2.dta", convert.dates = T,convert.factors = F)
household10b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC1.dta", convert.dates = T,convert.factors = F)
household10c<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC9A.dta", convert.dates = T,convert.factors = F)

land10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/AGSEC2A.dta", convert.dates = T,convert.factors = F)
fuel10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC10B.dta", convert.dates = T,convert.factors = F)
bank10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC13A.dta", convert.dates = T,convert.factors = F)
electricity10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC10A.dta", convert.dates = T,convert.factors = F)
goods10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC14.dta", convert.dates = T,convert.factors = F)
geographic10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/UNPS_Geovars_1011.dta", convert.dates = T,convert.factors = F)

land10$size<-land10$a2aq4
land.hh<-land10 %>% group_by(HHID)%>%
  summarise(landarea=sum(size, na.rm=T))%>% 
  arrange(HHID)%>%as.data.frame()

bank10$bank<-ifelse(bank10$h13q20==1,1,0)

electricity10$electricity<-ifelse(electricity10$h10q1==1,1,0)

household10c<-dplyr::rename(household10c, rooms=h9q3, roof=h9q4,
                            walls=h9q5, floor=h9q6, water=h9q7,
                            toilet=h9q22)

household10c$toilet_share<-ifelse(household10c$toilet==2|household10c$toilet==4|household10c$toilet==7,1,0)
household10c$toilet<-as.factor(ifelse(household10c$toilet==9, NA, household10c$toilet))
levels(household10c$toilet)<-c(2,2,2,2,1,3,3,0)

household10c$water<-as.factor(ifelse(household10c$water==96, NA, household10c$water))
levels(household10c$water)<-c(4,3,2,2,0,0,1,2,2)

household10c$roof<-as.factor(ifelse(household10c$roof==96, NA, household10c$roof))
levels(household10c$roof)<-c(0,0,1,2,3,3,2,3)

household10c$walls<-as.factor(ifelse(household10c$walls==96, NA, household10c$walls))
levels(household10c$walls)<-c(0,0,2,1,3,3,3,3)

household10c$floor<-as.factor(ifelse(household10c$floor==96, NA, household10c$floor))
levels(household10c$floor)<-c(0,0,2,2,2,2,1)

household10b$district<-as.factor(household10b$region)
levels(household10b$district)<-c("Kampala", "Central", "Eastern", "Northern", "Western")
household10b$DHSID<-household10b$comm

household10<-dplyr::rename(household10, relationship=h2q4)
household10$domestic<-ifelse(household10$relationship==9,1,0)
household10<-household10 %>% group_by(HHID)%>%
  summarise(members=n(),
            domestic=max(domestic, na.rm=T))%>% 
  arrange(HHID)%>% as.data.frame()

fuel10$fuel.lev<-as.factor(fuel10$h10q13_1)
levels(fuel10$fuel.lev)<-c("Firewood", "Dung", "Crop residue",
                           "Kerosene", "LPG", "Charcoal", "Solar", "Electricity")
fuel10$paste<-paste(fuel10$fuel.lev, fuel10$h10q13_2, sep="_")
fuel10$fuel<-ifelse(fuel10$paste=="Firewood_1", 1,
                    ifelse(fuel10$paste=="Dung_1", 0,
                           ifelse(fuel10$paste=="Crop residue_1", 0,
                                  ifelse(fuel10$paste=="Kerosene_1", 3,
                                         ifelse(fuel10$paste=="LPG_1", 3,
                                                ifelse(fuel10$paste=="Charcoal_1", 2,
                                                       ifelse(fuel10$paste=="Solar_1", 3, 
                                                              ifelse(fuel10$paste=="Electricity_1",4, NA))))))))

goods10$goods<-as.factor(goods10$h14q2)
levels(goods10$goods)<-c("house", "Other buildings", "land", 
                         "table", "fridge", "tv",
                         "radio", "generator", "solar_panel",
                         "bicycle", "motorcycle", "car", "boat", "other_transport",
                         "watch", "mobile", "computer", "internet", "other_electronics",
                         "other_assets", "other_assets", "other_assets")
newvars.g<-levels(goods10$goods)
goods10$own<-ifelse(goods10$h14q3==1,1,0)
goods2<-create_var(newvars.g, goods09, "own", "HHID")
goods2<- goods2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))
hhid.keep<-goods2$HHID
goods2<-as.data.frame(apply(goods2[,newvars.g], MARGIN=2, binary))
goods2$HHID<-hhid.keep

nps10<-merge(land.hh, bank10, by="HHID", all=TRUE)
nps10<-merge(nps09, electricity10, by="HHID", all=TRUE)
nps10<-merge(nps09, household10c, by="HHID", all=TRUE)
nps10<-merge(nps09, household10, by="HHID", all=TRUE)
nps10<-merge(nps09, fuel2, by="HHID", all=TRUE)
nps10<-merge(nps09, goods2, by="HHID", all=TRUE)
nps10<-merge(nps09, household10b, by="HHID", all=TRUE)

keep<-c("HHID", "electricity", "landarea", "bank", "water", "rooms", "roof", "fuel", "toilet", "walls", "floor", "water", "toilet_share", "members",
        "domestic", "fuel", newvars.g, "district", "DHSID") 
nps10<-nps10[,which(colnames(nps10)%in%keep)]

nps10$vehicle<-ifelse(nps10$car==1,2,
                      ifelse(nps10$motorcycle==1,1,
                             ifelse(nps10$bicycle==1,0,NA)))

possess<-levels(goods10$goods)

nps10$possessions<-(rowSums(nps10[,c(possess)], na.rm=TRUE))/length(possess)

nps10$survey<-rep("NPS10", nrow(nps10))
nps10$year<-rep(2010, nrow(nps10))

nps10$memroom<-nps10$rooms/nps10$members

NPS10<-merge(nps10, geographic10[c(1,2,3,8)], by = "HHID")

NPS10.coords<-data.frame(coords.x1=NPS10$lon_mod, coords.x2=NPS10$lat_mod)

missing.vars<-which(needed%ni%names(NPS10))
for(i in 1:length(missing.vars)){
  current_names<-names(NPS10)
  newvar<-rep(NA, nrow(NPS10))
  NPS10<-cbind(NPS10, newvar)
  colnames(NPS10)<-c(current_names, needed[missing.vars[i]])
}

#############################
# National Panel Survey 2011
#############################
household11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC2.dta", convert.dates = T,convert.factors = F)
household11b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC1.dta", convert.dates = T,convert.factors = F)
household11c<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC9A.dta", convert.dates = T,convert.factors = F)

geographic11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/UNPS_Geovars_1112.dta", convert.dates = T,convert.factors = F)

land11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/AGSEC2A.dta", convert.dates = T,convert.factors = F)
fuel11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC10C.dta", convert.dates = T,convert.factors = F)
electricity11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC10A.dta", convert.dates = T,convert.factors = F)
goods11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC14.dta", convert.dates = T,convert.factors = F)

land11$size<-land11$a2aq4
land.hh<-land11 %>% group_by(HHID)%>%
  summarise(landarea=sum(size, na.rm=T))%>% 
  arrange(HHID)%>%as.data.frame()

electricity11$electricity<-ifelse(electricity11$h10q1==1,1,0)

household11c<-dplyr::rename(household11c, rooms=h9q3, roof=h9q4,
                            walls=h9q5, floor=h9q6, water=h9q7,
                            toilet=h9q22)

household11c$toilet_share<-ifelse(household11c$toilet==2|household11c$toilet==4|household11c$toilet==7,1,0)
household11c$toilet<-as.factor(ifelse(household11c$toilet==9, NA, household11c$toilet))
levels(household11c$toilet)<-c(2,2,2,2,1,3,3,0)

household11c$water<-as.factor(ifelse(household11c$water==96, NA, household11c$water))
levels(household11c$water)<-c(4,3,2,2,0,0,1,2,2)

household11c$roof<-as.factor(ifelse(household11c$roof==96, NA, household11c$roof))
levels(household11c$roof)<-c(0,0,1,2,3,3,2,3)

household11c$walls<-as.factor(ifelse(household11c$walls==96, NA, household11c$walls))
levels(household11c$walls)<-c(0,0,2,1,3,3,3,3)

household11c$floor<-as.factor(ifelse(household11c$floor==96, NA, household11c$floor))
levels(household11c$floor)<-c(0,0,2,2,2,2,1)

household11b$district<-as.factor(household11b$region)
levels(household11b$district)<-c("Central", "Eastern", "Northern", "Western")
household11b$DHSID<-household11b$comm

household11<-dplyr::rename(household11, relationship=h2q4)
household11$domestic<-ifelse(household11$relationship==9,1,0)
household11<-household11 %>% group_by(HHID)%>%
  summarise(members=n(),
            domestic=max(domestic, na.rm=T))%>% 
  arrange(HHID)%>% as.data.frame()

fuel11$fuel.lev<-as.factor(fuel11$h10q13_1)
levels(fuel11$fuel.lev)<-c("Firewood", "Dung", "Crop residue",
                       "Kerosene", "LPG", "Charcoal", "Solar", "Electricity")
fuel11$paste<-paste(fuel11$fuel.lev, fuel11$h10q13_2, sep="_")
fuel11$fuel<-ifelse(fuel11$paste=="Firewood_1", 1,
                    ifelse(fuel11$paste=="Dung_1", 0,
                           ifelse(fuel11$paste=="Crop residue_1", 0,
                                  ifelse(fuel11$paste=="Kerosene_1", 3,
                                         ifelse(fuel11$paste=="LPG_1", 3,
                                                ifelse(fuel11$paste=="Charcoal_1", 2,
                                                       ifelse(fuel11$paste=="Solar_1", 3, 
                                                              ifelse(fuel11$paste=="Electricity_1",4, NA))))))))


goods11$goods<-as.factor(goods11$h14q2)
levels(goods11$goods)<-c("house", "Other buildings", "land", 
                         "table", "fridge", "tv",
                         "radio", "generator", "solar_panel",
                         "bicycle", "motorcycle", "car", "boat", "other_transport",
                         "watch", "mobile", "computer", "internet", "other_electronics",
                         "other_assets", "other_assets", "other_assets")
newvars.g<-levels(goods11$goods)
goods11$own<-ifelse(goods11$h14q3==1,1,0)
goods2<-create_var(newvars.g, goods09, "own", "HHID")
goods2<- goods2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))
hhid.keep<-goods2$HHID
goods2<-as.data.frame(apply(goods2[,newvars.g], MARGIN=2, binary))
goods2$HHID<-hhid.keep

nps11<-merge(land.hh, electricity11, by="HHID", all=TRUE)
nps11<-merge(nps11, household11c, by="HHID", all=TRUE)
nps11<-merge(nps11, household11, by="HHID", all=TRUE)
nps11<-merge(nps11, goods2, by="HHID", all=TRUE)
nps11<-merge(nps11, household11b, by="HHID", all=TRUE)

keep<-c("HHID", "electricity", "landarea", "water", "rooms", "roof", "fuel", "toilet", "walls", "floor", "water", "toilet_share", "members",
        "domestic", "fuel", newvars.g, "district", "DHSID") 
nps11<-nps11[,which(colnames(nps11)%in%keep)]

nps11$vehicle<-ifelse(nps11$car==1,2,
                      ifelse(nps11$motorcycle==1,1,
                             ifelse(nps11$bicycle==1,0,NA)))

possess<-levels(goods11$goods)

nps11$possessions<-(rowSums(nps11[,c(possess)], na.rm=TRUE))/length(possess)

nps11$survey<-rep("NPS11", nrow(nps11))
nps11$year<-rep(2011, nrow(nps11))

nps11$memroom<-nps11$rooms/nps11$members

NPS11<-merge(nps11, geographic11[c(1,2,3,7)], by = "HHID")

NPS11.coords<-data.frame(coords.x1=NPS11$lon_mod, coords.x2=NPS11$lat_mod)

missing.vars<-which(needed%ni%names(NPS11))
for(i in 1:length(missing.vars)){
  current_names<-names(NPS11)
  newvar<-rep(NA, nrow(NPS11))
  NPS11<-cbind(NPS11, newvar)
  colnames(NPS11)<-c(current_names, needed[missing.vars[i]])
}

########################
# Bind into one dataset
########################
merged_wealth<-rbind(dhs01[,c(needed)], dhs06[,c(needed)], 
                     dhs11[,c(needed)], dhs16[,c(needed)], MIS09[,c(needed)], 
                     AIS11[,c(needed)], 
                     MIS14[,c(needed)], MIS18[,c(needed)], NPS09[,c(needed)],
                     NPS10[,c(needed)], NPS11[,c(needed)])

merged_wealth$toilet_share<-ifelse(merged_wealth$toilet_share==1,0,1)

merged_coords<-rbind(DHS01.coords, DHS06.coords, 
                    DHS11.coords, DHS16.coords, MIS09.coords, MIS14.coords, MIS18.coords,
                    AIS11.coords, NPS09.coords, NPS10.coords, NPS11.coords)
  
merged_wealth$urban<-as.numeric(as.character(merged_wealth$urban))
merged_wealth$toilet<-as.numeric(as.character(merged_wealth$toilet))
merged_wealth$water<-as.numeric(as.character(merged_wealth$water))
merged_wealth$fuel<-as.numeric(as.character(merged_wealth$fuel))
merged_wealth$floor<-as.numeric(as.character(merged_wealth$floor))
merged_wealth$walls<-as.numeric(as.character(merged_wealth$walls))
merged_wealth$roof<-as.numeric(as.character(merged_wealth$roof))
merged_wealth$cluster<-as.numeric(as.character(merged_wealth$cluster))

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

save(merged_wealth, file="Data/Confounder_data/Uganda/merged_wealth.csv")

#################
#Factor analysis 
#################
clust_list<-merged_wealth[,1]
merged_wealth<-merged_wealth[,-1] #remove cluster id
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
discard<-c("hhid", "survey", "year", "district", "cluster")
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

#(1f) Handle factors with missing levels
newdat_keep[[1]]$floor<-ifelse(newdat_keep[[1]]$floor==2,1,0)
newdat_keep[[2]]$floor<-ifelse(newdat_keep[[2]]$floor==2,1,0)
newdat_keep[[1]]$walls<-ifelse(newdat_keep[[1]]$walls==3,2,newdat_keep[[1]]$walls)

#(2) Create correlation matrix
matrices<-list()
for(i in 1:length(newdat_keep)){
  data=newdat_keep[[i]]
  matrix<-mixedCor(data= data, ncat=2)
  matrices<-c(matrices, list(matrix))
}

matrices2<-list()
for(i in 1:length(newdat_keep)){
  data=newdat_keep[[i]]
  matrix<-mixedCor(data= data, ncat=2)
  matrices2<-c(matrices2, list(matrix))
}
corrmat_warnings<-warnings()
#(3) Do factor analysis (from documentation: principal components extraction using correlation method with one factor
#extracted, substitution of mean for missing values, estimation of the factor scores using the regression
#method)
fits<-list()
scores<-list()
for(i in 1:length(matrices2)){
  matrix<-matrices2[[i]]
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
    if(res<0){
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
wealth<-cbind(wealth, clust_list)
wealth$cluster<-wealth$clust_list
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
  nightRaster = raster(paste('Data/Predictor_data/Nighttime lights/nightlights_',nightYear,'.tif', sep=""))
  subset<-wealth.sp[which(wealth.sp@data$year==year),]
  nightlights.crop<-crop(nightRaster, ext)
  merged.trans <- spTransform(subset, crs(nightlights.crop))
  nightlights.extract<-raster::extract(nightlights.crop, merged.trans)
  idx<-which(wealth.sp@data$cluster%in%subset@data$cluster)
  wealth.sp@data$nightlights[idx]<-nightlights.extract
  wealth.sp@data$nightlights<-ifelse(wealth.sp@data$nightlights==255, NA, wealth.sp@data$nightlights)
}

#Take care of areas with missing district data
merge<-over(wealth.sp, districts)
merge$districts<-merge$NAME_4
wealth.sp@data$district<-ifelse(is.na(wealth.sp@data$district), as.character(merge$districts), as.character(wealth.sp@data$district))
idxNA = which (is.na(wealth.sp@data$district))
missing_district = wealth.sp[idxNA,]
plot(districts)
plot(missing_district, add=TRUE, col="red")

dists <- dist2Line(missing_district,districts)
missing_district@data$district<-as.factor(districts@data$NAME_4[dists[,4]])
wealth.sp@data$district<-ifelse(is.na(wealth.sp@data$district), as.character(missing_district@data$district), as.character(wealth.sp@data$district))

wealth.sp@data$district<-as.factor(wealth.sp@data$district)

poverty2<-raster("Data/Confounder_data/Uganda/uga10povcons200.tif")
merged.trans <- spTransform(wealth.sp, crs(poverty2))
pov.extract<-raster::extract(poverty2, merged.trans)
wealth.sp@data$prop_u2<-pov.extract

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

writeOGR(obj=wealth, dsn="Data/Confounder_data/Uganda/wealth", layer="wealth", driver="ESRI Shapefile", overwrite_layer = TRUE)
