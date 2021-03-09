library(readstata13)
library(dplyr)
library(mlr)
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
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
districts<-readOGR("Data/Exposure_data/Malawi/shapefiles/mwi_admbnda_adm2_nso_20181016", layer="mwi_admbnda_adm2_nso_20181016")
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
  finaldat<-data.frame(matrix(nrow=nRow, ncol=length(x)+1))
  colnames(finaldat)<-c(id, x)
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
    }else{
      row<-dat[which(dat$goods==x[i]),]
      cIdx<-which(names(row)==y)
      newvar<-ifelse(row[,cIdx]==1,1,0)
      idx<-which(names(row)==id)
      newdat<-as.data.frame(cbind(row[,idx], newvar))
      names(newdat)<-c(id, x[i])
      newdat[,2]<-as.integer(newdat[,2])
      newdat.hh<-aggregate(newdat[,2], by=list(newdat[,1]), FUN=max, na.rm=TRUE)
      finaldat[,i+1]<-newdat.hh[,2]
    }
  }
  return(finaldat)
}

'%ni%' <- Negate('%in%')

needed<-c("cluster", "survey", "year", "domestic", "land", "house", "urban", "fuel", "water", "toilet", "toilet_share", "floor", "walls", 
          "roof", "electricity", "bank", "landarea", "chairs", "table", "vehicle", "clock", "bucket", "clotheswash", 
          "dishwash", "fridge", "mobile", "landline", "tv", "computer", "radio", "sewingmachine", 
          "bed", "sofa", "lamp", "boat", "fishing_net", "mortar_pestle", "fan", "aircon", "tape_cd_dvd_hifi", "vcr", "stove_kp", "stove_eg", 
          "beer_drum", "coffee_table", "bureau", "desk", "iron", "satellite", "solar_panel", "generator", "possessions", "plough", "district")
################################################
# WHS
################################################
WHS<-read.dta13("Data/Exposure_data/Malawi/inputs/WHS_2003/MWI_2003_WHS_v01_M_Stata8/WHS-Malawi_F2.dta") 

all.neg <- function(x) -1*abs(x)
WHS<-WHS[which(!is.na(WHS$q0200_2)),]
WHS$q0200_2 <- all.neg(WHS$q0200_2)
lat<-paste(WHS$q0200_2, WHS$q0200_3, sep=".")
lon<-paste(WHS$q0201_2, WHS$q0201_3, sep=".")
WHS$lat<-as.numeric(as.character(lat))
WHS$lon<-as.numeric(as.character(lon))

WHS<- WHS %>% rename(rooms= q0700,
                     chairs = q0701,
                     table = q0702,
                     electricity = q0704,
                     clock = q0706, 
                     bucket = q0707,
                     clotheswash = q0708,
                     dishwash = q0709,
                     fridge = q0710,
                     landline = q0711,
                     mobile= q0712,
                     tv = q0713,
                     computer = q0714, 
                     radio = q0716, 
                     sewingmachine=q0717,
                     plough= q0719, 
                     bank=q0818, 
                     car=q0703,
                     bicycle=q0705,
                     animal_cart=q0715)

binvars<-as.data.frame(apply(WHS[,c("chairs", "table", "electricity", "clock", "bucket", "clotheswash", "dishwash", 
                                    "fridge", "landline", "mobile", "computer", "radio", "tv",
                                    "sewingmachine", "plough", "bank", "car", "bicycle", "animal_cart")], MARGIN=2, binary))

binvars$id<-WHS$id
binvars$lat<-WHS$lat
binvars$lon<-WHS$lon
binvars$rooms<-WHS$rooms
WHS<-binvars

WHS3<-read.dta13("Data/Exposure_data/Malawi/inputs/WHS_2003/MWI_2003_WHS_v01_M_Stata8/WHS-Malawi_F3.dta") 
WHS3$i<-rep(1, nrow(WHS3))
WHS3.hh<-WHS3 %>% group_by(id)%>%
  summarise(household.members=sum(i, na.rm=T))%>% 
  arrange(id)
WHS3.hh<-as.data.frame(WHS3.hh)

WHS5<-read.dta13("Data/Exposure_data/Malawi/inputs/WHS_2003/MWI_2003_WHS_v01_M_Stata8/WHS-Malawi_F5.dta") 
WHS5<- WHS5 %>% rename(toilet=q4045, water = q4042, floor = q4040, walls = q4041, fuel = q4047)

WHS5$toilet<-ifelse(WHS5$toilet==8, NA, WHS5$toilet)
WHS5$toilet<-as.factor(WHS5$toilet)
levels(WHS5$toilet)<-c("piped", "piped", "flush_latrine", "dry_covered_latrine", "dry_uncovered_latrine", "bucket", "open")
levels(WHS5$toilet)<-c(3,3,3,2,1,1,0)

WHS5$water<-as.factor(WHS5$water)
levels(WHS5$water)<-c("piped", "public_standpipe", "protected_well", "protected_well", "unprotected_well_spring", "raintank", "surface")
levels(WHS5$water)<-c(4,3,2,1,2,0)
  
WHS5$walls<-ifelse(WHS5$walls==6, NA, WHS5$walls)
WHS5$walls<-as.factor(WHS5$walls)
levels(WHS5$walls)<-c("cement_stone_brick_wood", "mudbrick", "thatch", "plastic_sheet", "metal_sheet")
levels(WHS5$walls)<-c(3,1,0,2,2)

WHS5$fuel<-ifelse(WHS5$fuel==10, NA, WHS5$fuel)
WHS5$fuel<-as.factor(WHS5$fuel)
levels(WHS5$fuel)<-c("gas", "electricity", "kerosene", "coal", "charcoal", "wood", "crop")
levels(WHS5$fuel)<-c(3,4,3,2,2,1,0)

WHS5used<-WHS5[,c("id", "toilet", "water", "walls", "fuel")]

WHS<-merge(WHS, WHS3.hh, by="id")

WHS<-merge(WHS, WHS5used, by="id")

WHS$rooms<-ifelse(WHS$rooms==0,1,WHS$rooms)
WHS$memroom<-WHS$rooms/WHS$household.members 

WHS.coords<-data.frame(coords.x1=WHS$lon, coords.x2=WHS$lat)
urbanNum<-getUrbanicity(year=2003, loc=WHS.coords, regMap=districts, country="Malawi")
urbanFac = factor(x = urbanNum,
                  levels = c(0, 1))
WHS$urban<-urbanFac

#Possessions
possess<-c("chairs", "table", "clock", "bucket", 
           "clotheswash", "dishwash", "fridge", "landline", "mobile", "tv", "computer", 
           "radio", "sewingmachine", "plough", "car", "bicycle", "animal_cart")

WHS$possessions<-(rowSums(WHS[,c(possess)], na.rm=TRUE))/length(possess)

WHS$survey<-rep("WHS", nrow(WHS))
WHS$year<-rep(2003, nrow(WHS))

missing.vars<-which(needed%ni%names(WHS))
for(i in 1:length(missing.vars)){
  current_names<-names(WHS)
  newvar<-rep(NA, nrow(WHS))
  WHS<-cbind(WHS, newvar)
  colnames(WHS)<-c(current_names, needed[missing.vars[i]])
}
  
################################################
# Technology Adoption and Risk Initiation Survey
################################################
matching<-read.dta13("Data/Exposure_data/Malawi/inputs/TARIHS_2006/MTARI_matching.dta", convert.dates = T,convert.factors = F)
qno<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/_qno.csv")
main<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/bmain_puf.csv")
land<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/blanduc_puf.csv")
bank<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/bsavins_puf.csv")
household<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/bhroster_puf.csv")
household$i<-rep(1, length(nrow(household)))
household$domestic<-ifelse(household$a2==18,1,0)

household<-merge(household[,c("qno","addrdp_n","gac_name","i", "domestic")], matching, by = "qno")

#Collapse individual-level data on on household to get household members
household<-household %>% group_by(qno)%>%
  summarise(household.members=sum(i, na.rm=T), 
            domestic=max(domestic, na.rm=T), district=head(addrdp_n,1), 
            gac_name=head(gac_name.y,1), village=head(vill_nam,1))%>% 
  arrange(qno)
household<-as.data.frame(household)
household$id<-rep(1, nrow(household))

#-------------#
# Geolocation
#-------------#
#Collapse on village
village<-household %>% group_by(village)%>%
  summarise(district=head(district,1), 
            gac_name=head(gac_name,1))%>% 
  arrange(village)
village<-as.data.frame(village)
dim(village)
village$vill_id = rownames(village)

#Source 1 for geolocation
populated.places<-read.csv("Data/Exposure_data/Malawi/shapefiles/Geonames/mi/populated_places.csv")
populated.places$names<-populated.places$FULL_NAME_RO
pp.coords<-cbind(populated.places$LONG, populated.places$LAT)
populated.places2<-as.data.frame(populated.places[,"names"])
colnames(populated.places2)<-c("name")
#Turn into shapefile
populated.places.sp <- SpatialPointsDataFrame(coords = pp.coords, data = populated.places2, proj4string = CRS(proj))
#Add district data
populated.places.districts<-over(populated.places.sp, districts)
populated.places.sp@data$DISTRICT<-populated.places.districts$ADM2_EN

#Source 2 for geolocation
villagesgeo<-readOGR("Data/Exposure_data/Malawi/shapefiles/villagesgeo", layer="villagesgeo", p4s=proj)
#Harmonize columns
villagesgeo@data$name= villagesgeo@data$NAME

#Source 3 for geolocation
villagespoints<-readOGR("Data/Exposure_data/Malawi/shapefiles/mw_villagespoints", layer="mw_villagespoints", p4s=proj)
#Harmonize columns
villagespoints@data$name= villagespoints@data$NAME

#Bind all 3
bind2<-bind(populated.places.sp, villagesgeo, villagespoints)
#Create location (village_district) variable
bind2$location= ifelse(!is.na(bind2$DISTRICT), paste(bind2$name, bind2$DISTRICT, sep="_"), as.character(bind2$name))

#Create location variable for TARIS data
village$location<-paste(village$village, village$district, sep="_")

#Merge TARIS data with geolocation data
village.gps<-merge(village, bind2, by="location", all.x=TRUE)

#Collapse on location (village)
village.gps<-village.gps %>% group_by(location)%>%
  summarise(village=head(village,1),
            district=head(district,1),
            coords.x1=mean(coords.x1, na.rm=T),
            coords.x2=mean(coords.x2, na.rm=T))%>% 
  arrange(location)
village.gps<-as.data.frame(village.gps)

#Check how successful geolocation was
nomiss<-village.gps[which(!is.na(village.gps$coords.x1)),]
nlevels(as.factor(nomiss$location))
missing<-village.gps[which(is.na(village.gps$coords.x1)),]
nlevels(as.factor(missing$location))

#Delete villages with no geolocation
villages<-village.gps[which(is.na(village.gps$coords.x2)==FALSE),] 

#-----------------#
# Wealth variables
#-----------------#
#Land ownership
household<-merge(household, land[,c("c2_area","c2__unit","c3","qno")], by = "qno")
household$land<-ifelse(household$c3==6,NA, household$c3)
household$land<-as.factor(household$land)
levels(household$land)<-c(1,1,1,0,0)
household$landarea<-ifelse(household$land==1,household$c2_area, 0)
household$landarea<-ifelse(household$c2__unit==1, household$landarea, household$landarea*2.47105)
household$land<-as.integer(as.character(household$land))

#Collapse
household<-household %>% group_by(qno)%>%
  summarise(land=max(land, na.rm=T), 
            landarea=sum(landarea, na.rm=T),
            household.members=head(household.members,1),
            domestic=max(domestic, na.rm=T), gac_name=head(gac_name,1),
            village=head(village,1), id=head(id,1))%>% 
  arrange(qno)
household<-as.data.frame(household)
household$land<-ifelse(household$land==0|household$land==1,household$land, NA)

#Home ownership (b7), number of rooms (b6), roof (b5), floor (b4), wall (b3), water (b8), electricity (b9)
household<-merge(household, main[c("qno", "b7", "b6", "b5", "b3", "b4", "b8", "b9")], by="qno")
household$house<-ifelse(household$b7==1,1,0)
household$rooms<-ifelse(household$b6==0,1,household$b6)
household$electricity<-ifelse(household$b9==1,1,0)

household$walls<-as.factor(ifelse(household$b3==8, NA, household$b3))
levels(household$walls)<-c("brick_stone", "burnt_bricks", "mud_bricks", "iron", "wood", "wood_mud", "mud")
levels(household$walls)<-c(3,3,1,2,2,2,0)

household$floor<-as.factor(household$b4)
levels(household$floor)<-c("earth", "wood_planks", "polished_wood_vinyl_tiles", "cement")
levels(household$floor)<-c(0,1,2,2)

household$roof<-as.factor(ifelse(household$b5==4, NA, household$b5))
levels(household$roof)<-c("straw_thatch", "sheet_metal", "tile_slate")
levels(household$roof)<-c(0,1,2)

household$water<-as.factor(ifelse(household$b8==9, NA, household$b8))
levels(household$water)<-c("piped_residence", "public_tap", "private_well", 
                           "public_well", "river/stream", "pond/lake", "rainwater", "borehole")
levels(household$water)<-c(4,3,2,1,0,0,1,2)

#Bank account
household<-merge(household, bank[c("qno", "h3_s1")], by="qno")
household$bank<-ifelse(household$h3_s1==1,1,0)

#Add get urbanicity last
villages<-merge(villages, matching, by.x="village", by.y="vill_nam")
TARIS<-merge(household, villages, by="qno")

TARIS.coords<-data.frame(coords.x1=TARIS$coords.x1, coords.x2=TARIS$coords.x2)
urbanNum<-getUrbanicity(year=2006, loc=TARIS.coords, regMap=districts, country="Malawi")
urbanFac = factor(x = urbanNum,
                  levels = c(0, 1))
TARIS$urban<-urbanFac

TARIS$memroom<-TARIS$rooms/TARIS$household.members
TARIS$survey<-rep("TARIS", nrow(TARIS))
TARIS$year<-rep(2006, nrow(TARIS))
TARIS$cluster<-TARIS$village.x

missing.vars<-which(needed%ni%names(TARIS))
for(i in 1:length(missing.vars)){
  current_names<-names(TARIS)
  newvar<-rep(NA, nrow(TARIS))
  TARIS<-cbind(TARIS, newvar)
  colnames(TARIS)<-c(current_names, needed[missing.vars[i]])
}

###########
# IHPS 2013
###########
ihps13<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_a_filt_13.dta", convert.dates = T,convert.factors = F) 
ihps13$district<-ifelse(ihps13$district==101, "Chitipa", 
                        ifelse(ihps13$district==102, "Karonga", 
                               ifelse(ihps13$district==103, "Nkhatabay", 
                                      ifelse(ihps13$district==104, "Rumphi", 
                                             ifelse(ihps13$district==107, "Mzuzu City", 
                                                    ifelse(ihps13$district==105, "Mzimba", 
                                                           ifelse(ihps13$district==201, "Kasungu", 
                                                                  ifelse(ihps13$district==202, "Nkhota Kota", 
                                                                         ifelse(ihps13$district==203, "Ntchisi", 
                                                                                ifelse(ihps13$district==204, "Dowa",
                                                                                       ifelse(ihps13$district==205, "Salima", 
                                                                                              ifelse(ihps13$district==206, "Lilongwe", 
                                                                                                     ifelse(ihps13$district==207, "Mchinji", 
                                                                                                            ifelse(ihps13$district==208, "Dedza", 
                                                                                                                   ifelse(ihps13$district==209, "Ntcheu", 
                                                                                                                          ifelse(ihps13$district==210, "Lilongwe City", 
                                                                                                                                 ifelse(ihps13$district==301, "Mangochi", 
                                                                                                                                        ifelse(ihps13$district==302, "Machinga", 
                                                                                                                                               ifelse(ihps13$district==303, "Zomba", 
                                                                                                                                                      ifelse(ihps13$district==304, "Chiradzulu", 
                                                                                                                                                             ifelse(ihps13$district==305, "Blantyre", 
                                                                                                                                                                    ifelse(ihps13$district==306, "Mwanza", 
                                                                                                                                                                           ifelse(ihps13$district==307, "Thyolo", 
                                                                                                                                                                                  ifelse(ihps13$district==308, "Mulanje", 
                                                                                                                                                                                         ifelse(ihps13$district==309, "Phalombe", 
                                                                                                                                                                                                ifelse(ihps13$district==310, "Chikwawa", 
                                                                                                                                                                                                       ifelse(ihps13$district==311, "Nsanje", 
                                                                                                                                                                                                              ifelse(ihps13$district==312, "Balaka", 
                                                                                                                                                                                                                     ifelse(ihps13$district==313, "Neno", 
                                                                                                                                                                                                                            ifelse(ihps13$district==314, "Zomba City", "Blantyre City"))))))))))))))))))))))))))))))




household_b<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_b_13.dta", convert.dates = T,convert.factors = F)
household_b$domestic<-ifelse(household_b$hh_b04==13,1,0)

household_b<-merge(household_b, ihps13[,c("y2_hhid", "district")], by="y2_hhid")

land<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/ag_mod_d_13.dta", convert.dates = T,convert.factors = F)
land$land<-ifelse(land$ag_d03==1|land$ag_d03==2|land$ag_d03==3|
                    land$ag_d03==4|land$ag_d03==5,1,0)

household_b$i<-rep(1, nrow(household_b))
household<-household_b %>% group_by(y2_hhid)%>%
  summarise(domestic=max(domestic, na.rm=T), 
            members=sum(i, na.rm=T), 
            district=head(district,1))%>% 
  arrange(y2_hhid)
household<-as.data.frame(household)

land.hh<-land %>% group_by(y2_hhid)%>%
  summarise(land=max(land, na.rm=T))%>% 
  arrange(y2_hhid)
land.hh<-as.data.frame(land.hh)
land.hh$land<-ifelse(land.hh$land==1,1,
                     ifelse(land.hh$land==0,0,NA))

household<-merge(household, land.hh, by="y2_hhid", all.x=TRUE)

household_f<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_f_13.dta", convert.dates = T,convert.factors = F)
household_f$house<-ifelse(household_f$hh_f01==1,1,0)
household_f$rooms<-ifelse(household_f$hh_f10==0,1,household_f$hh_f10)
household_f$electricity<-ifelse(household_f$hh_f19==1,1,0)
household_f$mobile<-ifelse(household_f$hh_f34>=1,1,0)

household_f$floor<-as.factor(ifelse(household_f$hh_f09==6,NA, household_f$hh_f09))
levels(household_f$floor)<-c("sand", "smoothed_mud", "smooth_cement", "wood", "tile")
levels(household_f$floor)<-c(0,0,2,1,2)

household_f$fuel<-as.factor(ifelse(household_f$hh_f12==10,NA, household_f$hh_f12))
levels(household_f$fuel)<-c("collected_firewood", "purchased_firewood", "paraffin", "electricity", "gas", 
                            "charcoal", "crop_residue")
levels(household_f$fuel)<-c(1,1,3,4,3,2,0)

household_f$walls<-as.factor(ifelse(household_f$hh_f07==9,NA, household_f$hh_f07))
levels(household_f$walls)<-c("grass", "mud", "compacted_earth", "mud_brick", "burnt_bricks", 
                             "concrete", "wood")
levels(household_f$walls)<-c(0,0,0,1,3,3,2)

household_f$roof<-as.factor(ifelse(household_f$hh_f08==6,NA, household_f$hh_f08))
levels(household_f$roof)<-c("grass", "iron_sheets", "clay_tiles", "plastic sheeting")
levels(household_f$roof)<-c(0,2,3,1)

household_f$water<-as.factor(ifelse(household_f$hh_f36==16,NA, household_f$hh_f36))
levels(household_f$water)<-c("piped_dwelling", "piped_yard", "standpipe", "open_well_yard", "open_well_public",
                             "protected_well_yard", "protected_well_public", "borehole", "spring", "river/stream", 
                             "pond/lake", "tanker_truck/bowser")
levels(household_f$water)<-c(4,4,3,1,1,2,2,2,0,0,0,2)

household_f$toilet<-as.factor(ifelse(household_f$hh_f41==6,NA, household_f$hh_f41))
levels(household_f$toilet)<-c("flush", "VIP", "traditional_latrine_roof", "traditional_latrine_nroof","none")
levels(household_f$toilet)<-c(3,2,2,1,0)

household_f$toilet_share<-ifelse(household_f$hh_f42==1,1,0)
household_f$bank<-ifelse(household_f$hh_f48==1,1,0)

household<-merge(household, household_f[,c("y2_hhid", "house", "rooms", "electricity", "floor", "fuel", "walls", "roof", "water", "toilet", "toilet_share", "bank", "mobile")], by="y2_hhid")

goods<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_l_13.dta", convert.dates = T,convert.factors = F)

goods$goods<-as.factor(goods$hh_l02)
levels(goods$goods)<-c("mortar_pestle", "bed", "table", "chairs", "fan", "aircon", 
                       "radio","tape_cd_dvd_hifi", "tv", "vcr", "sewingmachine", "stove_kp", 
                       "stove_eg", "fridge", "clotheswash", "bicycle", "motorcycle_scooter", "car", 
                       "mini_bus", "lorry", "beer_drum", "sofa", "coffee_table", "bureau", "lamp", "desk", "clock", 
                       "iron", "computer", "satellite", "solar_panel", "generator")

newvars<-levels(goods$goods)

goods2<-create_var(newvars, goods, "hh_l01", "y2_hhid")

goods2<- goods2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))

household<-merge(household, goods2, by="y2_hhid")

household$vehicle<-ifelse(household$car==1|household$mini_bus==1|household$lorry==1,2,
                          ifelse(household$motorcycle_scooter==1,1,
                                 ifelse(household$bicycle==1,0,NA)))

possess<-levels(goods$goods)

household$possessions<-(rowSums(household[,c(possess)], na.rm=TRUE))/length(possess)

household$survey<-rep("IHPS13", nrow(household))
household$year<-rep(2013, nrow(household))

household$memroom<-household$rooms/household$members

geo<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/HouseholdGeovariables_IHPS_13.dta", convert.dates = T,convert.factors = F) 
geo.ea<-geo %>% group_by(y2_hhid)%>%
  summarise(lat= mean (LAT_DD_MOD, na.rm=T), long = mean (LON_DD_MOD, na.rm=T))%>% 
  arrange(y2_hhid)
geo.ea<-as.data.frame(geo.ea)

household<-merge(household, geo.ea[,c("y2_hhid", "lat", "long")], by="y2_hhid")

IHPS13.coords<-data.frame(coords.x1=household$long, coords.x2=household$lat)

household_a<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_a_filt_13.dta", convert.dates = T,convert.factors = F) 
household_a$urban<-ifelse(household_a$reside==1,1,0)

household<-merge(household, household_a[,c("y2_hhid", "urban")], by="y2_hhid")

missing.vars<-which(needed%ni%names(household))
for(i in 1:length(missing.vars)){
  current_names<-names(household)
  newvar<-rep(NA, nrow(household))
  household<-cbind(household, newvar)
  colnames(household)<-c(current_names, needed[missing.vars[i]])
}

IHPS13<-household

###########
# IHPS 2016
###########
ihps16<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/ag_mod_r1_16.dta", convert.dates = T,convert.factors = F)
household16<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_a_filt_16.dta", convert.dates = T,convert.factors = F) 
ihps16<-merge(household16[,c(1,3:14,18)], ihps16, by = "y3_hhid")
ihps16$district<-ifelse(ihps16$district==101, "Chitipa", 
                        ifelse(ihps16$district==102, "Karonga", 
                               ifelse(ihps16$district==103, "Nkhatabay", 
                                             ifelse(ihps16$district==107, "Mzuzu City", 
                                                    ifelse(ihps16$district==105, "Mzimba", 
                                                           ifelse(ihps16$district==201, "Kasungu", 
                                                                  ifelse(ihps16$district==202, "Nkhota Kota", 
                                                                         ifelse(ihps16$district==203, "Ntchisi", 
                                                                                ifelse(ihps16$district==204, "Dowa",
                                                                                       ifelse(ihps16$district==205, "Salima", 
                                                                                              ifelse(ihps16$district==206, "Lilongwe", 
                                                                                                     ifelse(ihps16$district==207, "Mchinji", 
                                                                                                            ifelse(ihps16$district==208, "Dedza", 
                                                                                                                   ifelse(ihps16$district==209, "Ntcheu", 
                                                                                                                          ifelse(ihps16$district==210, "Lilongwe City", 
                                                                                                                                 ifelse(ihps16$district==301, "Mangochi", 
                                                                                                                                        ifelse(ihps16$district==302, "Machinga", 
                                                                                                                                               ifelse(ihps16$district==303, "Zomba", 
                                                                                                                                                      ifelse(ihps16$district==304, "Chiradzulu", 
                                                                                                                                                             ifelse(ihps16$district==305, "Blantyre", 
                                                                                                                                                                    ifelse(ihps16$district==306, "Mwanza", 
                                                                                                                                                                           ifelse(ihps16$district==307, "Thyolo", 
                                                                                                                                                                                  ifelse(ihps16$district==308, "Mulanje", 
                                                                                                                                                                                         ifelse(ihps16$district==309, "Phalombe", 
                                                                                                                                                                                                ifelse(ihps16$district==310, "Chikwawa", 
                                                                                                                                                                                                       ifelse(ihps16$district==311, "Nsanje", 
                                                                                                                                                                                                              ifelse(ihps16$district==312, "Balaka", 
                                                                                                                                                                                                                     ifelse(ihps16$district==313, "Neno", 
                                                                                                                                                                                                                            ifelse(ihps16$district==314, "Zomba City", "Blantyre City")))))))))))))))))))))))))))))




ihps16<-ihps16 %>% group_by(y3_hhid)%>%
  summarise(district=head(district,1))%>% 
  arrange(y3_hhid)
ihps16<-as.data.frame(ihps16)

household_b<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_b_16.dta", convert.dates = T,convert.factors = F)
household_b$domestic<-ifelse(household_b$hh_b04==13,1,0)

household_b<-merge(household_b, ihps16[,c("y3_hhid", "district")], by="y3_hhid")


land<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/ag_mod_b2_16.dta", convert.dates = T,convert.factors = F)
land$land<-ifelse(land$ag_b203==1|land$ag_b203==2|land$ag_b203==3|
                    land$ag_b203==4|land$ag_b203==12|land$ag_b203==13,1,0)

household_b$i<-rep(1, nrow(household_b))
household<-household_b %>% group_by(y3_hhid)%>%
  summarise(domestic=max(domestic, na.rm=T), 
            members=sum(i, na.rm=T), 
            district=head(district,1))%>% 
  arrange(y3_hhid)
household<-as.data.frame(household)

land.hh<-land %>% group_by(y3_hhid)%>%
  summarise(land=max(land, na.rm=T))%>% 
  arrange(y3_hhid)
land.hh<-as.data.frame(land.hh)
land.hh$land<-ifelse(land.hh$land==1,1,
                     ifelse(land.hh$land==0,0,NA))

household<-merge(household, land.hh, by="y3_hhid", all.x=TRUE)

household_f<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_f_16.dta", convert.dates = T,convert.factors = F)
household_f$house<-ifelse(household_f$hh_f01==1,1,0)
household_f$rooms<-ifelse(household_f$hh_f10==0,1,household_f$hh_f10)
household_f$electricity<-ifelse(household_f$hh_f19==1,1,0)
household_f$landline<-ifelse(household_f$hh_f31==1,1,0)
household_f$mobile<-ifelse(household_f$hh_f34>=1,1,0)

household_f$floor<-as.factor(ifelse(household_f$hh_f09==6,NA, household_f$hh_f09))
levels(household_f$floor)<-c("sand", "smoothed_mud", "smooth_cement", "tile", "wood")
levels(household_f$floor)<-c(0,0,2,2,1)

household_f$fuel<-as.factor(ifelse(household_f$hh_f12==10,NA, household_f$hh_f12))
levels(household_f$fuel)<-c("collected_firewood", "purchased_firewood", "paraffin", "electricity", "gas", 
                            "charcoal", "crop_residue", "saw_dust")
levels(household_f$fuel)<-c(1,1,3,4,3,2,0,0)

household_f$walls<-as.factor(ifelse(household_f$hh_f07==9,NA, household_f$hh_f07))
levels(household_f$walls)<-c("grass", "mud", "compacted_earth", "mud_brick", "burnt_bricks", 
                             "concrete", "wood", "iron_sheets")
levels(household_f$walls)<-c(0,0,0,1,3,3,2,2)

household_f$roof<-as.factor(ifelse(household_f$hh_f08==6,NA, household_f$hh_f08))
levels(household_f$roof)<-c("grass", "iron_sheets", "clay_tiles", "plastic sheeting")
levels(household_f$roof)<-c(0,2,3,1)

household_f$water<-as.factor(ifelse(household_f$hh_f36==16,NA, household_f$hh_f36))
levels(household_f$water)<-c("piped_dwelling", "piped_yard", "standpipe", "open_well_yard", "open_well_public",
                             "protected_well_yard", "protected_well_public", "borehole", "spring", "river/stream", 
                             "pond/lake", "dam", "tanker_truck/bowser")
levels(household_f$water)<-c(4,4,3,1,1,2,2,2,0,0,0,0,2)

household_f$toilet<-as.factor(ifelse(household_f$hh_f41==6,NA, household_f$hh_f41))
levels(household_f$toilet)<-c("flush", "VIP", "traditional_latrine_roof", "traditional_latrine_nroof","none")
levels(household_f$toilet)<-c(3,2,2,1,0)

household_f$toilet_share<-ifelse(household_f$hh_f42==1,1,0)
household_f$bank<-ifelse(household_f$hh_f48==1,1,0)

household<-merge(household, household_f[,c("y3_hhid", "house", "rooms", "electricity", "floor", "fuel", "walls", "roof", "water", "toilet", "toilet_share", "bank", "mobile", "landline")], by="y3_hhid")

goods<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_l_16.dta", convert.dates = T,convert.factors = F)

goods$goods<-as.factor(goods$hh_l02)
levels(goods$goods)<-c("mortar_pestle", "bed", "table", "chairs", "fan", "aircon", 
                       "radio","tape_cd_dvd_hifi", "tv", "vcr", "sewingmachine", "stove_kp", 
                       "stove_eg", "fridge", "clotheswash", "bicycle", "motorcycle_scooter", "car", 
                       "mini_bus", "lorry", "beer_drum", "sofa", "coffee_table", "bureau", "lamp", "desk", "clock", 
                       "iron", "computer", "satellite", "solar_panel", "generator", "radio_fancy")

newvars<-levels(goods$goods)

goods2<-create_var(newvars, goods, "hh_l01", "y3_hhid")

goods2<- goods2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))

household<-merge(household, goods2, by="y3_hhid")

household$vehicle<-ifelse(household$car==1|household$mini_bus==1|household$lorry==1,2,
                          ifelse(household$motorcycle_scooter==1,1,
                                 ifelse(household$bicycle==1,0,NA)))

possess<-levels(goods$goods)

household$possessions<-(rowSums(household[,c(possess)], na.rm=TRUE))/length(possess)

household$survey<-rep("IHPS16", nrow(household))
household$year<-rep(2016, nrow(household))

household$memroom<-household$rooms/household$members

geo<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/HouseholdGeovariablesIHPSY3.dta", convert.dates = T,convert.factors = F) 
geo.ea<-geo %>% group_by(y3_hhid)%>%
  summarise(lat= mean (lat_modified, na.rm=T), long = mean (lon_modified, na.rm=T))%>% 
  arrange(y3_hhid)
geo.ea<-as.data.frame(geo.ea)

household<-merge(household, geo.ea[,c("y3_hhid", "lat", "long")], by="y3_hhid")

IHPS16.coords<-data.frame(coords.x1=household$long, coords.x2=household$lat)

household_a<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_a_filt_16.dta", convert.dates = T,convert.factors = F) 
household_a$urban<-ifelse(household_a$reside==1,1,0)
household<-merge(household, household_a[,c("y3_hhid", "urban")], by="y3_hhid")

missing.vars<-which(needed%ni%names(household))
for(i in 1:length(missing.vars)){
  current_names<-names(household)
  newvar<-rep(NA, nrow(household))
  household<-cbind(household, newvar)
  colnames(household)<-c(current_names, needed[missing.vars[i]])
}

IHPS16<-household

###########
# IHS 2004
###########
ihs04_a<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/sec_a.dta")
ihs04.hh<-ihs04_a %>% group_by(case_id)%>%
  summarise(district=head(dist,1), hhid=head(hhid,1))%>% 
  arrange(case_id)
ihs04.hh<-as.data.frame(ihs04.hh)

household_b<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/sec_b.dta")
household_b$domestic<-ifelse(household_b$b04==13,1,0)
household_b<-merge(ihs04.hh[,c("case_id", "district")], household_b, by="case_id")

land<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/sec_o.dta")
land$land<-ifelse(land$o09==1|land$o09==2|land$o09==3|
                    land$o09==4|land$o09==5,1,0)

household_b$i<-rep(1, nrow(household_b))
household<-household_b %>% group_by(hhid)%>%
  summarise(domestic=max(domestic, na.rm=T), 
            members=sum(i, na.rm=T), 
            EACODE=head(psu,1), 
            district=head(district,1), case_id=head(case_id, 1))%>% 
  arrange(hhid)
household<-as.data.frame(household)

household<-merge(ihs04.hh, household, by="hhid")

land.hh<-land %>% group_by(hhid)%>%
  summarise(land=max(land, na.rm=T))%>% 
  arrange(hhid)
land.hh<-as.data.frame(land.hh)
land.hh$land<-ifelse(land.hh$land==1,1,
                     ifelse(land.hh$land==0,0,NA))

household<-merge(household, land.hh, by="hhid", all.x=TRUE)

household_g<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/sec_g.dta")
household_g$house<-ifelse(household_g$g01=="Owned",1,0)
household_g$rooms<-ifelse(household_g$g11==0,1,household_g$g11)
household_g$electricity<-ifelse(household_g$g20=="Yes",1,0)
household_g$landline<-ifelse(household_g$g25=="Yes",1,0)
household_g$mobile<-ifelse(household_g$g28=="Yes",1,0)

household_g$floor<-household_g$g10
levels(household_g$floor)<-c(0,0,2,1,2,NA)

household_g$fuel<-household_g$g13
levels(household_g$fuel)<-c(1,1,3,4,3,2,0,0,0,NA)

household_g$walls<-household_g$g08
levels(household_g$walls)<-c(0,0,0,1,3,3,2,2,NA)

household_g$roof<-household_g$g09
levels(household_g$roof)<-c(0,2,3,3,1,NA)

household_g$water<-household_g$g30
levels(household_g$water)<-c(4,4,3,2,2,2,1,1,0,0,NA)

household_g$toilet<-household_g$g36
levels(household_g$toilet)<-c(3,2,2,1,0,NA)

household_g$toilet_share<-ifelse(household_g$g37==1,1,0)

household_g<-household_g %>% group_by(hhid)%>%
  summarise(house=head(house,1),
            rooms=head(rooms,1),
            electricity=head(electricity,1),
            floor=head(floor,1),
            walls=head(walls,1),
            fuel=head(fuel,1),
            roof=head(roof,1),
            water=head(water,1),
            toilet=head(toilet,1),
            toilet_share=head(toilet_share,1), 
            mobile=max(mobile, na.rm=T),
            landline=max(landline, na.rm=T))%>% 
  arrange(hhid)
household_g<-as.data.frame(household_g)

household<-merge(household, household_g[,c("hhid", "house", "rooms", "electricity", "floor", "fuel", "walls", "roof", "water", "toilet", "toilet_share", "mobile", "landline")], by="hhid")

goods<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/sec_m1.dta")

goods$goods<-as.factor(goods$m0a)
levels(goods$goods)<-c("aircon","bed","bicycle","car",                       
                       "chairs","stove_eg", "fan","stove_kp",
                       "lorry","mini_bus","mortar_pestle","motorcycle_scooter","radio",      
                       "fridge","sewingmachine","table",                     
                       "tape_cd_dvd_hifi","tv","clotheswash")

newvars<-levels(goods$goods)

goods2<-create_var(newvars, goods, "m01a", "hhid")

goods2<- goods2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))

goods$urban<-ifelse(goods$reside==0,1,0)

goods.hh<-goods %>% group_by(hhid)%>%
  summarise(urban=max(urban, na.rm=T))%>% 
  arrange(hhid)
goods.hh<-as.data.frame(goods.hh)

household<-merge(household, goods2, by="hhid")
household<-merge(household, goods.hh[,c("hhid", "urban")], by="hhid")

household$vehicle<-ifelse(household$car==1|household$mini_bus==1|household$lorry==1,2,
                          ifelse(household$motorcycle_scooter==1,1,
                                 ifelse(household$bicycle==1,0,NA)))

possess<-levels(goods$goods)

household$possessions<-(rowSums(household[,c(possess)], na.rm=TRUE))/length(possess)

household$survey<-rep("IHS04", nrow(household))
household$year<-rep(2004, nrow(household))

household$memroom<-household$rooms/household$members

geo<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2004/2004_IHS_EA.csv")
geo.ea<-geo %>% group_by(EACODE)%>%
  summarise(lat= mean (Lat, na.rm=T), long = mean (Long, na.rm=T))%>% 
  arrange(EACODE)
geo.ea<-as.data.frame(geo.ea)

household<-merge(household, geo.ea[,c("EACODE", "lat", "long")], by="EACODE")

IHS04.coords<-data.frame(coords.x1=household$long, coords.x2=household$lat)

household$cluster<-household$EACODE

missing.vars<-which(needed%ni%names(household))
for(i in 1:length(missing.vars)){
  current_names<-names(household)
  newvar<-rep(NA, nrow(household))
  household<-cbind(household, newvar)
  colnames(household)<-c(current_names, needed[missing.vars[i]])
}

IHS04<-household

###########
# IHS 2010
###########
household_b<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Household/hh_mod_b.csv")
household_b$domestic<-ifelse(household_b$hh_b04==13,1,0)

land<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Agriculture/ag_mod_d.csv")
land$land<-ifelse(land$ag_d03==1|land$ag_d03==2|land$ag_d03==3|
                    land$ag_d03==4|land$ag_d03==5,1,0)

household_b$i<-rep(1, nrow(household_b))
household<-household_b %>% group_by(case_id)%>%
  summarise(domestic=max(domestic, na.rm=T), 
            members=sum(i, na.rm=T), 
            ea_id=head(ea_id,1))%>% 
  arrange(case_id)
household<-as.data.frame(household)

land.hh<-land %>% group_by(case_id)%>%
  summarise(land=max(land, na.rm=T))%>% 
  arrange(case_id)
land.hh<-as.data.frame(land.hh)
land.hh$land<-ifelse(land.hh$land==1,1,
                     ifelse(land.hh$land==0,0,NA))

household<-merge(household, land.hh, by="case_id", all.x=TRUE)

household_f<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Household/hh_mod_f.csv")
household_f$house<-ifelse(household_f$hh_f01==1,1,0)
household_f$rooms<-ifelse(household_f$hh_f10==0,1,household_f$hh_f10)
household_f$electricity<-ifelse(household_f$hh_f19==1,1,0)
household_f$landline<-ifelse(household_f$hh_f31==1,1,0)
household_f$mobile<-ifelse(household_f$hh_f34>=1,1,0)

household_f$floor<-as.factor(ifelse(household_f$hh_f09==6,NA, household_f$hh_f09))
levels(household_f$floor)<-c("sand", "smoothed_mud", "smooth_cement", "wood", "tile")
levels(household_f$floor)<-c(0,0,2,1,2)

household_f$fuel<-as.factor(ifelse(household_f$hh_f12==10|household_f$hh_f12==14,NA, household_f$hh_f12))
levels(household_f$fuel)<-c("collected_firewood", "purchased_firewood", "paraffin", "electricity", "gas", 
                            "charcoal", "crop_residue", "saw_dust", "animal_waste")
levels(household_f$fuel)<-c(1,1,3,4,3,2,0,0,0)

household_f$walls<-as.factor(ifelse(household_f$hh_f07==9,NA, household_f$hh_f07))
levels(household_f$walls)<-c("grass", "mud", "compacted_earth", "mud_brick", "burnt_bricks", 
                             "concrete", "wood", "iron_sheets")
levels(household_f$walls)<-c(0,0,0,1,3,3,2,2)

household_f$roof<-as.factor(ifelse(household_f$hh_f08==6,NA, household_f$hh_f08))
levels(household_f$roof)<-c("grass", "iron_sheets", "clay_tiles", "concrete", "plastic_sheeting")
levels(household_f$roof)<-c(0,2,3,3,1)

household_f$water<-as.factor(ifelse(household_f$hh_f36==16,NA, household_f$hh_f36))
levels(household_f$water)<-c("piped_dwelling", "piped_yard", "standpipe", "open_well_yard", "open_well_public",
                             "protected_well_yard", "protected_well_public", "borehole", "spring", "river/stream", 
                             "pond/lake", "dam")
levels(household_f$water)<-c(4,4,3,1,1,2,2,2,0,0,0,0)

household_f$toilet<-as.factor(ifelse(household_f$hh_f41==6,NA, household_f$hh_f41))
levels(household_f$toilet)<-c("flush", "VIP", "traditional_latrine_roof", "traditional_latrine_nroof","none")
levels(household_f$toilet)<-c(3,2,2,1,0)

household_f$toilet_share<-ifelse(household_f$hh_f42==1,1,0)

household<-merge(household, household_f[,c("case_id", "house", "rooms", "electricity", "floor", "fuel", "walls", "roof", "water", "toilet", "toilet_share", "mobile", "landline")], by="case_id")

goods<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Household/hh_mod_l.csv")

goods$goods<-as.factor(goods$hh_l02)
levels(goods$goods)<-c("mortar_pestle", "bed", "table", "chairs", "fan", "aircon", 
                       "radio","tape_cd_dvd_hifi", "tv", "vcr", "sewingmachine", "stove_kp", 
                       "stove_eg", "fridge", "clotheswash", "bicycle", "motorcycle_scooter", "car", 
                       "mini_bus", "lorry", "beer_drum", "sofa", "coffee_table", "bureau", "lamp", "desk", "clock", 
                       "iron", "computer", "satellite", "solar_panel", "generator")

newvars<-levels(goods$goods)

goods2<-create_var(newvars, goods, "hh_l01", "case_id")

goods2<- goods2 %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))

household<-merge(household, goods2, by="case_id")

household$vehicle<-ifelse(household$car==1|household$mini_bus==1|household$lorry==1,2,
                          ifelse(household$motorcycle_scooter==1,1,
                                 ifelse(household$bicycle==1,0,NA)))

possess<-levels(goods$goods)

household$possessions<-(rowSums(household[,c(possess)], na.rm=TRUE))/length(possess)

household$survey<-rep("IHS10", nrow(household))
household$year<-rep(2010, nrow(household))

household$memroom<-household$rooms/household$members

geo<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Geovariables/HH_level/householdgeovariables.csv")
geo.ea<-geo %>% group_by(ea_id)%>%
  summarise(lat= mean (lat_modified, na.rm=T), long = mean (lon_modified, na.rm=T))%>% 
  arrange(ea_id)
geo.ea<-as.data.frame(geo.ea)

com10<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Community/com_ca.csv")
com10$district<-com10$com_ca01
com10$district<-ifelse(com10$district==101, "Chitipa", 
                       ifelse(com10$district==102, "Karonga", 
                              ifelse(com10$district==103, "Nkhatabay", 
                                     ifelse(com10$district==104, "Rumphi", 
                                            ifelse(com10$district==107, "Mzuzu City", 
                                                   ifelse(com10$district==105, "Mzimba", 
                                                          ifelse(com10$district==201, "Kasungu", 
                                                                 ifelse(com10$district==202, "Nkhota Kota", 
                                                                        ifelse(com10$district==203, "Ntchisi", 
                                                                               ifelse(com10$district==204, "Dowa",
                                                                                      ifelse(com10$district==205, "Salima", 
                                                                                             ifelse(com10$district==206, "Lilongwe", 
                                                                                                    ifelse(com10$district==207, "Mchinji", 
                                                                                                           ifelse(com10$district==208, "Dedza", 
                                                                                                                  ifelse(com10$district==209, "Ntcheu", 
                                                                                                                         ifelse(com10$district==210, "Lilongwe City", 
                                                                                                                                ifelse(com10$district==301, "Mangochi", 
                                                                                                                                       ifelse(com10$district==302, "Machinga", 
                                                                                                                                              ifelse(com10$district==303, "Zomba", 
                                                                                                                                                     ifelse(com10$district==304, "Chiradzulu", 
                                                                                                                                                            ifelse(com10$district==305, "Blantyre", 
                                                                                                                                                                   ifelse(com10$district==306, "Mwanza", 
                                                                                                                                                                          ifelse(com10$district==307, "Thyolo", 
                                                                                                                                                                                 ifelse(com10$district==308, "Mulanje", 
                                                                                                                                                                                        ifelse(com10$district==309, "Phalombe", 
                                                                                                                                                                                               ifelse(com10$district==310, "Chikwawa", 
                                                                                                                                                                                                      ifelse(com10$district==311, "Nsanje", 
                                                                                                                                                                                                             ifelse(com10$district==312, "Balaka", 
                                                                                                                                                                                                                    ifelse(com10$district==313, "Neno", 
                                                                                                                                                                                                                           ifelse(com10$district==314, "Zomba City", "Blantyre City"))))))))))))))))))))))))))))))



geo.ea<-merge(geo.ea, com10[,c("ea_id", "district")], by="ea_id")
household<-merge(household, geo.ea[,c("ea_id", "lat", "long", "district")], by="ea_id")

IHS10.coords<-data.frame(coords.x1=household$long, coords.x2=household$lat)

household_a<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Household/hh_mod_a_filt.csv")
household_a$urban<-ifelse(household_a$reside==1,1,0)

household<-merge(household, household_a[,c("urban", "case_id")], by="case_id")

missing.vars<-which(needed%ni%names(household))
for(i in 1:length(missing.vars)){
  current_names<-names(household)
  newvar<-rep(NA, nrow(household))
  household<-cbind(household, newvar)
  colnames(household)<-c(current_names, needed[missing.vars[i]])
}

IHS10<-household

###########
# IHS 2016
###########
household_b<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HH_MOD_B.dta")
household_b$domestic<-ifelse(household_b$hh_b04==13,1,0)

land<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/AG_MOD_B2.dta")
land$land<-ifelse(land$ag_b203==1|land$ag_b203==2|land$ag_b203==3|
                    land$ag_b203==4|land$ag_b203==13,1,0)

household_b$i<-rep(1, nrow(household_b))
household<-household_b %>% group_by(case_id)%>%
  summarise(domestic=max(domestic, na.rm=T), 
            members=sum(i, na.rm=T))%>% 
  arrange(case_id)
household<-as.data.frame(household)

land.hh<-land %>% group_by(case_id)%>%
  summarise(land=max(land, na.rm=T))%>% 
  arrange(case_id)
land.hh<-as.data.frame(land.hh)
land.hh$land<-ifelse(land.hh$land==1,1,
                     ifelse(land.hh$land==0,0,NA))

household<-merge(household, land.hh, by="case_id", all.x=TRUE)

household_f<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HH_MOD_F.dta")
household_f$house<-ifelse(household_f$hh_f01=="OWNED",1,0)
household_f$rooms<-household_f$hh_f10
household_f$electricity<-ifelse(household_f$hh_f19=="YES",1,0)
household_f$landline<-ifelse(household_f$hh_f31=="YES",1,0)
household_f$mobile<-ifelse(household_f$hh_f34>=1,1,0)
household_f$bank<-ifelse(household_f$hh_f48=="YES"|household_f$hh_f50=="YES",1,0)

household_f$floor<-household_f$hh_f09
levels(household_f$floor)<-c(0,0,2,1,2,NA)

household_f$fuel<-household_f$hh_f12
levels(household_f$fuel)<-c(1,1,3,4,3,2,0,0,0,NA)

household_f$walls<-household_f$hh_f07
levels(household_f$walls)<-c(0,0,0,1,3,3,2,2,NA)

household_f$roof<-household_f$hh_f08
levels(household_f$roof)<-c(0,2,3,3,1,NA)

household_f$water<-household_f$hh_f36
levels(household_f$water)<-c(4,4,3,1,1,2,2,2,0,0,0,0,2,2,2,NA)

household_f$toilet<-household_f$hh_f41
levels(household_f$toilet)<-c(3,2,2,1,0,NA)

household_f$toilet_share<-ifelse(household_f$hh_f42==1,1,0)

household_f<-household_f %>% group_by(case_id)%>%
  summarise(house=head(house,1),
            rooms=head(rooms,1),
            electricity=head(electricity,1),
            floor=head(floor,1),
            walls=head(walls,1),
            fuel=head(fuel,1),
            roof=head(roof,1),
            water=head(water,1),
            toilet=head(toilet,1),
            toilet_share=head(toilet_share,1), 
            mobile=max(mobile, na.rm=T),
            landline=max(landline, na.rm=T), 
            bank=max(bank, na.rm=T))%>% 
  arrange(case_id)
household_f<-as.data.frame(household_f)

household_f<- household_f %>% 
  mutate_all(.funs = function(x) replace(x, which(x == -Inf | x == "N/A"), NA))

household<-merge(household, household_f[,c("case_id", "house", "rooms", "electricity", "floor", "fuel", "walls", "roof", "water", "toilet", "toilet_share", "mobile", "landline", "bank")], by="case_id")

goods<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HH_MOD_L.dta")

goods$goods<-as.factor(goods$hh_l02)
levels(goods$goods)<-c("mortar_pestle", "bed", "table", "chairs", "fan", "aircon", 
                      "radio","tape_cd_dvd_hifi", "tv", "vcr", "sewingmachine", "stove_kp", 
                      "stove_eg", "fridge", "clotheswash", "bicycle", "motorcycle_scooter", "car", 
                      "mini_bus", "lorry", "beer_drum", "sofa", "coffee_table", "bureau", "lamp", "desk", "clock", 
                      "iron", "computer", "satellite", "solar_panel", "generator", "radio_fancy")

newvars<-levels(goods$goods)

goods2<-create_var(newvars, goods, "hh_l01", "case_id")

household<-merge(household, goods2, by="case_id")

household$vehicle<-ifelse(household$car==1|household$mini_bus==1|household$lorry==1,2,
                          ifelse(household$motorcycle_scooter==1,1,
                                 ifelse(household$bicycle==1,0,NA)))

possess<-levels(goods$goods)

household$possessions<-(rowSums(household[,c(possess)], na.rm=TRUE))/length(possess)

household$survey<-rep("IHS16", nrow(household))
household$year<-rep(2016, nrow(household))

household$memroom<-household$rooms/household$members

urb<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/IHS4 Consumption Aggregate.dta")
urb$urban<-ifelse(urb$urban=="URBAN",1,0)
household<-merge(household, urb[,c("urban", "case_id")], by="case_id")

geo<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HouseholdGeovariables_stata11/HouseholdGeovariablesIHS4.dta")
geo.ea<-geo %>% group_by(case_id)%>%
  summarise(lat= mean (lat_modified, na.rm=T), long = mean (lon_modified, na.rm=T))%>% 
  arrange(case_id)
geo.ea<-as.data.frame(geo.ea)

ihs16.indv<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HH_MOD_A_FILT.dta", convert.dates = T,convert.factors = F)
ihs16.hh<-ihs16.indv %>% group_by(case_id)%>%
  summarise(district=head(district,1))%>% 
  arrange(case_id)
ihs16.hh<-as.data.frame(ihs16.hh)

household<-merge(household, geo.ea[,c("case_id", "lat", "long")], by="case_id")
household<-merge(household, ihs16.hh[,c("case_id", "district")], by="case_id")

household$district<-ifelse(household$district==101, "Chitipa", 
                             ifelse(household$district==102, "Karonga", 
                                    ifelse(household$district==103, "Nkhatabay", 
                                           ifelse(household$district==104, "Rumphi", 
                                                  ifelse(household$district==107, "Mzuzu City", 
                                                         ifelse(household$district==105, "Mzimba", 
                                                                ifelse(household$district==201, "Kasungu", 
                                                                       ifelse(household$district==202, "Nkhota Kota", 
                                                                              ifelse(household$district==203, "Ntchisi", 
                                                                                     ifelse(household$district==204, "Dowa",
                                                                                            ifelse(household$district==205, "Salima", 
                                                                                                   ifelse(household$district==206, "Lilongwe", 
                                                                                                          ifelse(household$district==207, "Mchinji", 
                                                                                                                 ifelse(household$district==208, "Dedza", 
                                                                                                                        ifelse(household$district==209, "Ntcheu", 
                                                                                                                               ifelse(household$district==210, "Lilongwe City", 
                                                                                                                                      ifelse(household$district==301, "Mangochi", 
                                                                                                                                             ifelse(household$district==302, "Machinga", 
                                                                                                                                                    ifelse(household$district==303, "Zomba", 
                                                                                                                                                           ifelse(household$district==304, "Chiradzulu", 
                                                                                                                                                                  ifelse(household$district==305, "Blantyre", 
                                                                                                                                                                         ifelse(household$district==306, "Mwanza", 
                                                                                                                                                                                ifelse(household$district==307, "Thyolo", 
                                                                                                                                                                                       ifelse(household$district==308, "Mulanje", 
                                                                                                                                                                                              ifelse(household$district==309, "Phalombe", 
                                                                                                                                                                                                     ifelse(household$district==310, "Chikwawa", 
                                                                                                                                                                                                            ifelse(household$district==311, "Nsanje", 
                                                                                                                                                                                                                   ifelse(household$district==312, "Balaka", 
                                                                                                                                                                                                                          ifelse(household$district==313, "Neno", 
                                                                                                                                                                                                                                 ifelse(household$district==314, "Zomba City", "Blantyre City"))))))))))))))))))))))))))))))






IHS16.coords<-data.frame(coords.x1=household$long, coords.x2=household$lat)

missing.vars<-which(needed%ni%names(household))
for(i in 1:length(missing.vars)){
  current_names<-names(household)
  newvar<-rep(NA, nrow(household))
  household<-cbind(household, newvar)
  colnames(household)<-c(current_names, needed[missing.vars[i]])
}

IHS16<-household

IHS16$urban<-ifelse(IHS16$urban=="URBAN",1,0)

###########
# DHS
###########
dhs<-read.csv("Data/Exposure_data/Malawi/inputs/DHS/idhs_00017.csv")
dhs$district<-as.factor(dhs$GEOALT_MW2010_2016)

levels(dhs$district)<-c("Chitipa", "Karonga", 
                             "Nkhatabay and Likoma", "Rumphi", 
                             "Mzimba and Mzuzu City", "Kasungu", 
                             "Nkhotakota", "Ntchisi", "Dowa", 
                             "Salima", "Lilongwe", "Mchinji", 
                             "Dedza", "Ntcheu", "Mangochi", "Machinga", 
                             "Zomba", "Chradzulu", "Blantyre", 
                             "Mwanza", "Thyolo", "Mulange", 
                             "Phalombe", "Chikwawa", "Ndanjge", 
                             "Balaka", "Neno")

dhs$hhid<-paste(dhs$DHSID, dhs$HHNUMALL)
dhs$urban<-ifelse(dhs$URBANHH==2,0,1)
dhs$members<-dhs$DEJURENO
dhs$rooms<-ifelse(dhs$SLEEPROOMS==0,1,dhs$SLEEPROOMS)
dhs$bank<-ifelse(dhs$BANKACC==1,1,
                 ifelse(dhs$BANKACC==8,NA, 0))
dhs$survey<-paste("DHS", dhs$YEAR, sep="")
dhs$year<-dhs$YEAR
dhs$domestic<-ifelse(dhs$HHRELATE==41,1,0)
dhs$land<-ifelse(dhs$AGLANDYN==1,1,
                 ifelse(dhs$AGLANDYN==8,NA, 0))
dhs$toilet_share<-ifelse(dhs$TOILETSHAREYN==1,1,
                 ifelse(dhs$TOILETSHAREYN==8|dhs$TOILETSHAREYN==9,NA, 0))
dhs$electricity<-ifelse(dhs$ELECTRCHH==1,1,
                         ifelse(dhs$ELECTRCHH==8,NA, 0))
dhs$cluster<-dhs$DHSID

dhs$fuel<-as.factor(ifelse(dhs$COOKFUEL==800|dhs$COOKFUEL==995|dhs$COOKFUEL==998,NA, dhs$COOKFUEL))
levels(dhs$fuel)<-c("Electricity","LPG_natural gas","LPG", "Natural_gas", "Kerosene", "Biogas",
                    "Coal_lignite", "Wood", "Charcoal", "Firewood_straw", "Straw_shrub_grass",
                    "Dung","Agricultural_crop_based")
levels(dhs$fuel)<-c(4,3,3,3,3,3,2,1,2,1,0,0,0)

dhs$toilet<-as.factor(ifelse(dhs$TOILETTYPE==5000|dhs$TOILETTYPE==9998,NA, dhs$TOILETTYPE))
levels(dhs$toilet)<-c("NO FACILITY","Unspecified type of flush toilet",
                      "Own flush toilet (unspecified type)","Shared flush toilet (unspecified type)",
                      "Flush to piped sewer system",
                      "Flush to septic tank","Flush to pit latrine","Flush to somewhere else",
                      "Flush, dont know where","Composting toilet","Traditional pit toilet or latrine",
                      "Pit latrine without slab or open pit","Pit latrine with slab","Ventilated improved pit latrine",
                      "Bucket toilet","Hanging latrine over water source")
levels(dhs$toilet)<-c(0,3,3,3,3,3,3,3,3,2,1,1,2,2,1,1)

dhs$water<-as.factor(ifelse(dhs$DRINKWTR==6000|dhs$DRINKWTR==9998,NA, dhs$DRINKWTR))

levels(dhs$water)<-c("Piped into own dwelling","Piped into own yard/plot","Public tap/standpipe",
                     "Piped into neighbors dwelling/yard", "Unprotected/open well",
                     "Open well in own dwelling/yard/plot","Open public well",
                     "Protected well","Protected well in own dwelling/yard/plot","Protected public well",
                     "Tube well or borehole","Unspecified public well","Spring","Protected spring/surface water",
                     "Unprotected spring/surface water","River/dam/lake/ponds/streams/canal/irrigation channel",
                     "River, stream","Pond, lake","Dam",
                     "RAINWATER","Tanker truck","Cart with small tank","Bottled water")
levels(dhs$water)<-c(3,3,2,2,0,0,0,2,2,2,2,0,0,2,0,0,0,0,0,2,1,1,1)

dhs$roof<-as.factor(ifelse(dhs$ROOF==400|dhs$ROOF==998,NA, dhs$ROOF))
levels(dhs$roof)<-c("No roof","Thatch/palm leaf",
                    "Grass, thatch","Sod","Rustic mat","Palm/bamboo","Palm, bamboo, grass",
                    "Wood planks","Cardboard","Iron sheets","Metal","Asbetos","Wood",
                    "Cement","Ceramic tiles","Iron and tiles","Calamine/cement fiber","Roofing shingles")
levels(dhs$roof)<-c(0,0,0,0,1,1,1,1,1,2,2,3,3,3,3,3,3,3)

dhs$walls<-as.factor(ifelse(dhs$WALL==400|dhs$WALL==998,NA, dhs$WALL))
levels(dhs$walls)<-c("No walls","Cane/palm/trunks","Dirt",
                     "Bamboo with mud","Stone with mud","Plywood",
                     "Reused wood","Poles and mud","Cardboard",
                     "Uncovered adobe","Unburnt bricks","Cement/concrete",
                     "Bricks","Finished/burnt bricks","Cement blocks",
                     "Wood planks/shingles","Stone with lime/cement","Covered adobe")
levels(dhs$walls)<-c(0,0,0,0,3,2,2,0,0,1,1,3,3,3,3,3,3,3)

dhs$floor<-as.factor(ifelse(dhs$FLOOR==400|dhs$FLOOR==998,NA, dhs$FLOOR))
levels(dhs$floor)<-c("Earth, sand","Dung","Wood planks","Palm/bamboo",
                     "Broken bricks","Parquet/polished wood","Vinyl/asphalt strips/linoleum",
                     "Tiles/mosaic","Ceramic tiles","Cement/concrete","Carpet","Bricks")
levels(dhs$floor)<-c(0,0,1,1,1,2,2,2,2,2,2,2)

dhs2<- dhs %>% rename(landline= HHPHONEHH,
                     mobile = MOBPHONE,
                     computer=PC,
                     bicycle=BIKEHH,
                     car=CARHH,
                     motorcycle_scooter=MOTORCYCLHH,
                     boat=BOATWMOTOR,
                     fridge=FRIDGEHH,
                     bed=BED,
                     radio=RADIOHH,
                     animal_cart=DRAWNCART,
                     tv=TVHH,
                     clock=WATCHCLOCK,
                     sofa=SOFA,
                     table=TABLECHARYN,
                     lamp1=PARALAMPYN,
                     lamp2=OPARALAMPYN)

dhs2$chairs<-dhs2$table

ridmiss.dat<-as.data.frame(apply(dhs2[,c("landline", "mobile", "computer", 
                                    "bicycle", "car", "motorcycle_scooter", "bed",
                                    "boat", "fridge", "radio", "tv", "animal_cart",
                                    "clock", "sofa", "table", "chairs", "lamp1", "lamp2")], MARGIN=2, ridmiss, yIdx="1", mIdx="8"))

dhs<-cbind(dhs, ridmiss.dat)

dhs$lamp<-max(dhs$lamp1, dhs$lamp2, na.rm=T)
dhs$vehicle<-ifelse(dhs$car==1,2,
                          ifelse(dhs$motorcycle_scooter==1,1,
                                 ifelse(dhs$bicycle==1|dhs$animal_cart==1,0,NA)))

possess<-c("landline", "mobile", "computer", 
           "bicycle", "car", "motorcycle_scooter", "bed",
           "boat", "fridge", "radio", "tv", "animal_cart",
           "clock", "sofa", "table", "chairs", "lamp1", "lamp2")

dhs$possessions<-(rowSums(dhs[,c(possess)], na.rm=TRUE))/length(possess)

dhs$memroom<-dhs$rooms/dhs$members

missing.vars<-which(needed%ni%names(dhs))
for(i in 1:length(missing.vars)){
  current_names<-names(dhs)
  newvar<-rep(NA, nrow(dhs))
  dhs<-cbind(dhs, newvar)
  colnames(dhs)<-c(current_names, needed[missing.vars[i]])
}

DHS00<-subset(dhs, year==2000)
DHS04<-subset(dhs, year==2004)
DHS10<-subset(dhs, year==2010)
DHS16<-subset(dhs, year==2016)

dhs$fuel<-ifelse(dhs$year==2000, NA, dhs$fuel)

mshapedhs00<-readOGR("Data/Exposure_data/Malawi/inputs/DHS/shapefiles/MWGE43FL", layer="MWGE43FL") 
mshapedhs04<-readOGR("Data/Exposure_data/Malawi/inputs/DHS/shapefiles/MWGE4BFL", layer="MWGE4BFL") 
mshapedhs10<-readOGR("Data/Exposure_data/Malawi/inputs/DHS/shapefiles/MWGE62FL", layer="MWGE62FL") 
mshapedhs16<-readOGR("Data/Exposure_data/Malawi/inputs/DHS/shapefiles/MWGE7AFL", layer="MWGE7AFL") 

dhs00<-merge(dhs, mshapedhs00,  by="DHSID")
dhs04<-merge(dhs, mshapedhs04, by="DHSID")
dhs10<-merge(dhs, mshapedhs10,  by="DHSID")
dhs16<-merge(dhs, mshapedhs16,  by="DHSID")

DHS00.coords<-data.frame(coords.x1=dhs00$coords.x1, coords.x2=dhs00$coords.x2)
DHS04.coords<-data.frame(coords.x1=dhs04$coords.x1, coords.x2=dhs04$coords.x2)
DHS10.coords<-data.frame(coords.x1=dhs10$coords.x1, coords.x2=dhs10$coords.x2)
DHS16.coords<-data.frame(coords.x1=dhs16$coords.x1, coords.x2=dhs16$coords.x2)

############
#MICS 2013
############
mics5<-read.spss("Data/Exposure_data/Malawi/inputs/MICS5/Malawi MICS 2013-14 SPSS Datasets/hh.sav", to.data.frame=TRUE)
mics5<-dplyr::rename(mics5, cluster_num= HH1, hh_num= HH2, members=HH11, electricity=HC8A, 
                     solar_panel=HC8A1, radio=HC8B, tv=HC8C, landline=HC8D, fridge=HC8E, 
                     lamp=HC8F, bed=HC8G, table_chairs=HC8H, lamp2=HC8I, lamp3=HC8J, 
                     computer=HC8K, watch=HC9A, mobile=HC9B, bicycle=HC9C, motorcycle_scooter=HC9D,
                     animal_cart=HC9E, car_truck=HC9F, boat=HC9G, canoe=HC9H, fishing_net=HC9I,
                     house=HC10, land=HC11, landarea_units=HC12A, landarea=HC12B, fuel=HC6,
                     bank=HC15, water=WS1, toilet=WS8, toilet_share=WS9, rooms=HC2, roof=HC4, 
                     floor=HC3, walls=HC5)
mics5[mics5 %in% c("<NA>")] <- NA
mics5<-mics5[which(mics5$HH9=="Completed"),]

mics5geo<-read.csv("Data/Exposure_data/Malawi/inputs/MICS5/Malawi MICS 2013-14 SPSS Datasets/2013_14_MICS_EA.csv")
mics5geo$district<-mics5geo$DISTNAME

mics5<-merge(mics5, mics5geo, by.x="cluster_num", by.y="CLUSTERNO", all.x=TRUE)
conversion_hectares<-as.numeric(as.character(mics5$landarea))*2.47105
conversion_football<-as.numeric(as.character(mics5$landarea))*1.7643308639291
mics5$landarea<-ifelse(mics5$landarea_units=="Acres", mics5$landarea, 
                       ifelse(mics5$landarea_units=="Hectares", conversion_hectares,
                              ifelse(mics5$landarea_units=="Football pitches", conversion_football, NA)))
mics5$land<-ifelse(mics5$land=="Yes",1, 
                          ifelse(mics5$land=="No",0,NA))
mics5$landarea<-ifelse(mics5$land==1,mics5$landarea, 
                       ifelse(mics5$land==0, 0,NA))

mics5$urban<-ifelse(mics5$HH6=="Rural", 0,1)
mics5$survey=rep("MICS2013", nrow(mics5))
mics5$year<-rep(2013, nrow(mics5))

levels(mics5$toilet)<-c(3,3,3,3,3,2,2,1,2,2,1,1,0,NA,NA)
mics5$toilet_share<-ifelse(mics5$toilet_share=="Yes",1,
                           ifelse(mics5$toilet_share=="No", 0, NA))

levels(mics5$water)<-c(4,4,3,3,2,2,1,2,0,2,2,2,0,2,NA,NA)

levels(mics5$fuel)<-c(4,3,3,3,3,2,2,1,0,0,0,NA,NA,NA)

levels(mics5$roof)<-c(0,0,1,1,1,2,3,3,3,NA,NA)

levels(mics5$walls)<-c(0,0,0,0,3,1,2,0,2,1,3,3,3,3,3,2,NA,NA)

levels(mics5$floor)<-c(0,0,1,1,2,2,2,2,2,NA,NA)

mics5$electricity<-ifelse(mics5$electricity=="Yes",1, 
                          ifelse(mics5$electricity=="No",0,NA))

mics5$house<-ifelse(mics5$house=="Own",1,
                    ifelse(mics5$house=="Missing", NA,0))

mics5$bank<-ifelse(mics5$bank=="Yes",1, 
                          ifelse(mics5$bank=="No",0,NA))

goods<-c("solar_panel","radio", "tv", "landline", "fridge", 
         "lamp", "bed", "table_chairs", "lamp2", "lamp3", 
         "computer", "watch", "mobile", "bicycle", "motorcycle_scooter",
         "animal_cart", "car_truck", "boat", "canoe", "fishing_net")

goods<-mics5[,c("hh_num", "solar_panel","radio", "tv", "landline", "fridge", 
                "lamp", "bed", "table_chairs", "lamp2", "lamp3", 
                "computer", "watch", "mobile", "bicycle", "motorcycle_scooter",
                "animal_cart", "car_truck", "boat", "canoe", "fishing_net")]
  
mics5<-mics5[,c("hh_num", "cluster_num", "members", "electricity", "floor", "fuel", "walls", "roof", "water", "toilet", "toilet_share", "rooms", "house", "land", "landarea", "district", "bank", "urban", "Lat", "Long", "EACODE", "survey","year")]

ridmiss.dat<-as.data.frame(apply(goods[,c("solar_panel","radio", "tv", "landline", "fridge", 
                                          "lamp", "bed", "table_chairs", "lamp2", "lamp3", 
                                          "computer", "watch", "mobile", "bicycle", "motorcycle_scooter",
                                          "animal_cart", "car_truck", "boat", "canoe", "fishing_net")], MARGIN=2, ridmiss, yIdx="Yes", mIdx="Missing"))

mics5<-cbind(mics5, ridmiss.dat)

mics5$vehicle<-ifelse(mics5$car_truck==1,2,
                    ifelse(mics5$motorcycle_scooter==1,1,
                           ifelse(mics5$bicycle==1|mics5$animal_cart==1,0,NA)))

mics5$lamp<-max(mics5$lamp1, mics5$lamp2, mics5$lamp3, na.rm=T)

mics5$table<-mics5$chairs<-mics5$table_chairs

possess<-c("solar_panel","radio", "tv", "landline", "fridge", 
           "lamp", "bed", "table_chairs", "lamp2", "lamp3", 
           "computer", "watch", "mobile", "bicycle", "motorcycle_scooter",
           "animal_cart", "car_truck", "boat", "canoe", "fishing_net")

mics5$possessions<-(rowSums(mics5[,c(possess)], na.rm=TRUE))/length(possess)

mics5$rooms<-ifelse(mics5$rooms=="Missing", NA, as.numeric(as.character(mics5$rooms)))

mics5$memroom<-as.numeric(as.character(mics5$rooms))/as.numeric(as.character(mics5$members))

keepvars<-names(mics5)[names(mics5)%ni%c("table_chairs", "lamp2", "lamp3", "members", "rooms")]

codeline=paste0(keepvars[1]," =head(",keepvars[1],", 1)")
for(i in 2:length(keepvars)){
  codeline_new=paste0(keepvars[i]," =head(",keepvars[i],", 1)")
  codeline=paste(codeline,codeline_new, sep=", ")
}

mics5$HHID<-paste(mics5$cluster_num, mics5$hh_num, sep="_")

mics5$clock<-mics5$watch

mics5.l<-read.spss("Data/Exposure_data/Malawi/inputs/MICS5/Malawi MICS 2013-14 SPSS Datasets/hl.sav", to.data.frame=TRUE)
mics5.l$domestic<-ifelse(mics5.l$HL3==14,1,0)

mics5.l$HHID<-paste(mics5.l$HH1, mics5.l$HH2, sep="_")

mics5.l.hh<-mics5.l %>% group_by(HHID)%>%
  summarise(domestic=max(domestic,na.rm=T))%>% 
  arrange(HHID)
mics5.l.hh<-as.data.frame(mics5.l.hh)

mics5<-merge(mics5, mics5.l.hh, by="HHID")

missing.vars<-which(needed%ni%names(mics5))
for(i in 1:length(missing.vars)){
  current_names<-names(mics5)
  newvar<-rep(NA, nrow(mics5))
  mics5<-cbind(mics5, newvar)
  colnames(mics5)<-c(current_names, needed[missing.vars[i]])
}

MICS13.coords<-data.frame(coords.x1=mics5$Long, coords.x2=mics5$Lat)

MICS13<-mics5

###########
# MIS 2012
###########
household12<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2012/MWHR6HDT/MWHR6HFL.DTA", convert.dates = T,convert.factors = F)
individual12<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2012/MWIR6HDT/MWIR6HFL.DTA", convert.dates = T,convert.factors = F)
geospatial12<-read.csv("Data/Exposure_data/Malawi/inputs/MIS_2012/MWGC6AFL/MWGC6AFL.csv")
geographic12<-readOGR("Data/Exposure_data/Malawi/inputs/MIS_2012/MWGE6AFL", layer="MWGE6AFL") 

household12<-dplyr::rename(household12, DHSCLUST= hv001, hh_id=hv002, 
                           members=hv012, rooms=hv216, land=hv244, landarea=hv245, 
                           bank=hv247,electricity=hv206, radio=hv207, tv=hv208,
                           fridge=hv209, bicycle=hv210, motorcycle_scooter=hv211,
                           car_truck=hv212, floor=hv213, walls=hv214, roof=hv215,
                           landline=hv221, fuel=hv226, toilet_share=hv225,
                           mobile=hv243a, clock=hv243b, animal_cart=hv243c, 
                           boat=hv243d, toilet=hv205, water=hv201)

household12$HHID<-paste(household12$DHSCLUST, household12$hh_id, sep="_")
household12$district<-as.factor(household12$shdistr)
levels(household12$district)<-c("Chitipa", "Karonga", "Nkhatabay", "Rumphi", "Mzimba", "Mzuzu City", 
                                "Kasungu", "Nkhotakota", "Ntchisi", "Dowa", "Salima", "Lilongwe", 
                                "Mchinji", "Dedza", "Ntcheu", "Lilongwe city", "Mangochi", "Manchinga", 
                                "Zomba", "Chiradzulu", "Blantyre", "Thyolo", "Mulanje", "Phalombe", 
                                "Chikwawa", "Nsanje", "Balaka", "Zomba city", "Blantyre city")

household12<-household12[,c("DHSCLUST", "hh_id", 
         "members", "rooms", "land", "landarea", 
         "bank","electricity", "radio", "tv",
         "fridge", "bicycle", "motorcycle_scooter",
         "car_truck", "floor", "walls", "roof",
         "landline", "fuel", "toilet_share",
         "mobile", "clock", "animal_cart", 
         "boat", "toilet", "water", "district")]

household12<-merge(household12, geographic12@data[,c(1,4,15)], by = "DHSCLUST")
household12$urban<-ifelse(household12$URBAN_RURA=="R",0,1)

conversion_hectares<-as.numeric(as.character(household12$landarea))*2.47105
household12$landarea<-ifelse(household12$land==1, conversion_hectares, NA)

household12$survey=rep("MIS2012", nrow(household12))
household12$year<-rep(2012, nrow(household12))

household12$toilet<-as.factor(ifelse(household12$toilet==96, NA, household12$toilet))
levels(household12$toilet)<-c("Flush to piped sewer system", "VIP", "Pit latrine with slab", 
                              "Pit latrine without slab/open pit","No facility/bush/field",
                              "Composting toilet","Hanging toilet/latrine")
levels(household12$toilet)<-c(3,2,2,1,0,2,1)

household12$water<-as.factor(ifelse(household12$water==96, NA, household12$water))
levels(household12$water)<-c("Piped into dwelling", "Piped to yard/plot", "Public tap/standpipe", 
                             "Tube well or borehole", "Protected well", "Unprotected well", 
                             "Protected spring", "Unprotected spring", "River/dam/lake/ponds/stream/canal/irrigation channel",
                             "Tanker truck")
levels(household12$water)<-c(4,4,3,2,2,1,2,0,0,2)

household12$fuel<-as.factor(ifelse(household12$fuel==96|household12$fuel==95, NA, household12$fuel))
levels(household12$fuel)<-c("Electricity","LPG/ natural gas","Coal, lignite","Charcoal","Wood","Straw/shrubs/grass")
levels(household12$fuel)<-c(4,3,2,2,1,0)

household12$roof<-as.factor(ifelse(household12$roof==96, NA, household12$roof))
levels(household12$roof)<-c("No roof","Thatch / palm leaf","Palm / bamboo / grass","Wood planks","Cardboard","Iron sheets","Wood","Calamine / cement fiber","Ceramic tiles","Cement","Roofing shingles")
levels(household12$roof)<-c(0,0,1,1,1,2,3,3,3,3,3)

household12$walls<-as.factor(ifelse(household12$walls==96, NA, household12$walls))
levels(household12$walls)<-c("No walls","Cane / palm / trunks","Dirt","Bamboo/tree trunks with mud",
                             "Stone with mud","Plywood","Cardboard","Reused wood","Cement","Stone with lime / cement",
                             "Burnt bricks","Unburnt bricks","Cement blocks")
levels(household12$walls)<-c(0,0,0,0,3,2,0,2,3,3,3,1,3)

household12$floor<-as.factor(ifelse(household12$floor==96, NA, household12$floor))
levels(household12$floor)<-c("Earth, sand","Dung","Broken bricks","Parquet, polished wood","Vinyl, asphalt strips",
                             "Ceramic tiles","Cement","Carpet")
levels(household12$floor)<-c(0,0,1,2,2,2,2,2)

household12$rooms<-ifelse(household12$rooms==0,1,household12$rooms)

household12$vehicle<-ifelse(household12$car_truck==1,2,
                      ifelse(household12$motorcycle_scooter==1,1,
                             ifelse(household12$bicycle==1|household12$animal_cart==1,0,NA)))

possess<-c("radio", "tv",
           "fridge", "bicycle", "motorcycle_scooter",
           "car_truck", 
           "landline", 
           "mobile", "clock", "animal_cart", 
           "boat")

household12$possessions<-(rowSums(household12[,c(possess)], na.rm=TRUE))/length(possess)

household12$memroom<-as.numeric(as.character(household12$rooms))/as.numeric(as.character(household12$members))

missing.vars<-which(needed%ni%names(household12))
for(i in 1:length(missing.vars)){
  current_names<-names(household12)
  newvar<-rep(NA, nrow(household12))
  household12<-cbind(household12, newvar)
  colnames(household12)<-c(current_names, needed[missing.vars[i]])
}

household<-merge(household12, geographic12[,c("DHSCLUST", "LATNUM", "LONGNUM")], by="DHSCLUST")

MIS12.coords<-data.frame(coords.x1=household$LONGNUM, coords.x2=household$LATNUM)

MIS12<-household

###########
# MIS 2014
###########
household14<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2014/MWHR71DT/MWHR71FL.DTA", convert.dates = T,convert.factors = F)
household.members14<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2014/MWPR71DT/MWPR71FL.DTA", convert.dates = T,convert.factors = F)
individual14<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2014/MWIR71DT/MWIR71FL.DTA", convert.dates = T,convert.factors = F)
geospatial14<-read.csv("Data/Exposure_data/Malawi/inputs/MIS_2014/MWGC71FL/MWGC71FL.csv")
geographic14<-readOGR("Data/Exposure_data/Malawi/inputs/MIS_2014/MWGE71FL", layer="MWGE71FL") 

household14<-dplyr::rename(household14, DHSCLUST= hv001, hh_id=hv002, 
                           members=hv012, rooms=hv216, land=hv244, landarea=hv245, 
                           bank=hv247,electricity=hv206, radio=hv207, tv=hv208,
                           fridge=hv209, bicycle=hv210, motorcycle_scooter=hv211,
                           car_truck=hv212, floor=hv213, walls=hv214, roof=hv215,
                           landline=hv221, fuel=hv226, toilet_share=hv225,
                           mobile=hv243a, clock=hv243b, animal_cart=hv243c, 
                           boat=hv243d, toilet=hv205, water=hv201)

household14$HHID<-paste(household14$DHSCLUST, household14$hh_id, sep="_")

household14$district<-as.factor(household14$shdistr)
levels(household14$district)<-c("Chitipa", "Karonga", "Nkhatabay", "Rumphi", "Mzimba", "Mzuzu City", 
                                "Kasungu", "Nkhotakota", "Ntchisi", "Dowa", "Salima", "Lilongwe", 
                                "Mchinji", "Dedza", "Ntcheu", "Lilongwe city", "Mangochi", "Manchinga", 
                                "Zomba", "Chiradzulu", "Blantyre", "Thyolo", "Mulanje", "Phalombe", 
                                "Chikwawa", "Nsanje", "Balaka", "Zomba city", "Blantyre city")

household14<-household14[,c("DHSCLUST", "hh_id", 
                            "members", "rooms", "land", "landarea", 
                            "bank","electricity", "radio", "tv",
                            "fridge", "bicycle", "motorcycle_scooter",
                            "car_truck", "floor", "walls", "roof",
                            "landline", "fuel", "toilet_share",
                            "mobile", "clock", "animal_cart", 
                            "boat", "toilet", "water", "district")]

household14<-merge(household14, geographic14@data[,c(1,4,15)], by = "DHSCLUST")
household14$urban<-ifelse(household14$URBAN_RURA=="R",0,1)

conversion_hectares<-as.numeric(as.character(household14$landarea))*2.47105
household14$landarea<-ifelse(household14$land==1, conversion_hectares, NA)

household14$survey=rep("MIS2014", nrow(household14))
household14$year<-rep(2014, nrow(household14))

household14$toilet<-as.factor(ifelse(household14$toilet==96, NA, household14$toilet))
levels(household14$toilet)<-c("Flush to piped sewer system", "VIP", "Pit latrine with slab", 
                              "Pit latrine without slab/open pit","No facility/bush/field",
                              "Composting toilet","Hanging toilet/latrine")
levels(household14$toilet)<-c(3,2,2,1,0,2,1)

household14$water<-as.factor(ifelse(household14$water==96, NA, household14$water))
levels(household14$water)<-c("Piped into dwelling", "Piped to yard/plot", "Public tap/standpipe", 
                             "Tube well or borehole", "Protected well", "Unprotected well", 
                             "Unprotected spring", "River/dam/lake/ponds/stream/canal/irrigation channel",
                             "Bottled water")
levels(household14$water)<-c(4,4,3,2,2,1,0,0,2)

household14$fuel<-as.factor(ifelse(household14$fuel==95, NA, household14$fuel))
levels(household14$fuel)<-c("Electricity","LPG/ natural gas","Coal, lignite","Charcoal","Wood","Straw/shrubs/grass")
levels(household14$fuel)<-c(4,3,2,2,1,0)

household14$roof<-as.factor(ifelse(household14$roof==96, NA, household14$roof))
levels(household14$roof)<-c("No roof","Thatch / palm leaf","Rustic mat", "Palm / bamboo / grass","Wood planks","Iron sheets","Wood","Calamine / cement fiber","Ceramic tiles","Cement","Roofing shingles")
levels(household14$roof)<-c(0,0,1,1,1,2,3,3,3,3,3)

household14$walls<-as.factor(ifelse(household14$walls==96, NA, household14$walls))
levels(household14$walls)<-c("No walls","Cane / palm / trunks","Dirt","Bamboo/tree trunks with mud",
                             "Stone with mud","Plywood","Cardboard","Cement","Stone with lime / cement",
                             "Burnt bricks","Unburnt bricks","Cement blocks", "Wood planks / shingles")
levels(household14$walls)<-c(0,0,0,0,2,2,0,2,2,2,1,2,2)

household14$floor<-as.factor(ifelse(household14$floor==96, NA, household14$floor))
levels(household14$floor)<-c("Earth, sand","Dung","Wood planks", "Broken bricks","Parquet, polished wood","Vinyl, asphalt strips",
                             "Ceramic tiles","Cement","Carpet")
levels(household14$floor)<-c(0,0,1,1,2,2,2,2,2)

household14$rooms<-ifelse(household14$rooms==0,1,household14$rooms)

household14$vehicle<-ifelse(household14$car_truck==1,2,
                            ifelse(household14$motorcycle_scooter==1,1,
                                   ifelse(household14$bicycle==1|household14$animal_cart==1,0,NA)))

possess<-c("radio", "tv",
           "fridge", "bicycle", "motorcycle_scooter",
           "car_truck", 
           "landline", 
           "mobile", "clock", "animal_cart", 
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
# MIS 2017
###########
household17<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWHR7IDT/MWHR7IFL.DTA", convert.dates = T,convert.factors = F)
household.members17<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWPR7IDT/MWPR7IFL.DTA", convert.dates = T,convert.factors = F)
individual17<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWIR7IDT/MWIR7IFL.DTA", convert.dates = T,convert.factors = F)
geographic17<-readOGR("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWGE7IFL", layer="MWGE7IFL") 
geospatial17<-read.csv("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWGC7JFL/MWGC7JFL.csv")

household17<-dplyr::rename(household17, DHSCLUST= hv001, hh_id=hv002, 
                           members=hv012, rooms=hv216, land=hv244, landarea=hv245, 
                           bank=hv247,electricity=hv206, radio=hv207, tv=hv208,
                           fridge=hv209, bicycle=hv210, motorcycle_scooter=hv211,
                           car_truck=hv212, floor=hv213, walls=hv214, roof=hv215,
                           landline=hv221, fuel=hv226, toilet_share=hv225,
                           mobile=hv243a, clock=hv243b, animal_cart=hv243c, 
                           boat=hv243d, toilet=hv205, water=hv201)

household17$HHID<-paste(household17$DHSCLUST, household17$hh_id, sep="_")

household17<-household17[,c("DHSCLUST", "hh_id", 
                            "members", "rooms", "land", "landarea", 
                            "bank","electricity", "radio", "tv",
                            "fridge", "bicycle", "motorcycle_scooter",
                            "car_truck", "floor", "walls", "roof",
                            "landline", "fuel", "toilet_share",
                            "mobile", "clock", "animal_cart", 
                            "boat", "toilet", "water")]

household17<-merge(household17, geographic17@data[,c(1,4,15)], by = "DHSCLUST")
household17$urban<-ifelse(household17$URBAN_RURA=="R",0,1)

conversion_hectares<-as.numeric(as.character(household17$landarea))*2.47105
household17$landarea<-ifelse(household17$land==1, conversion_hectares, NA)

household17$survey=rep("MIS2017", nrow(household17))
household17$year<-rep(2017, nrow(household17))

household17$toilet<-as.factor(ifelse(household17$toilet==96, NA, household17$toilet))
levels(household17$toilet)<-c("Flush to piped sewer system", "Flush to septic tank","Flush to pit latrine",
                              "Flush, dont know where","VIP", "Pit latrine with slab", 
                              "Pit latrine without slab/open pit","No facility/bush/field",
                              "Composting toilet","Hanging toilet/latrine")
levels(household17$toilet)<-c(3,3,3,3,2,2,1,0,2,1)
household17$toilet<-as.integer(as.character(household17$toilet))

household17$water<-as.factor(ifelse(household17$water==96, NA, household17$water))
levels(household17$water)<-c("Piped into dwelling", "Piped to yard/plot", "Piped to neighbor", "Public tap/standpipe", 
                             "Tube well or borehole", "Protected well", "Unprotected well", "Protected spring",
                             "Unprotected spring", "River/dam/lake/ponds/stream/canal/irrigation channel",
                             "Tanker truck", "Cart with small tank")
levels(household17$water)<-c(3,3,2,2,1,1,0,1,0,0,1,1)
household17$water<-as.integer(as.character(household17$water))

household17$fuel<-as.factor(ifelse(household17$fuel==96|household17$fuel==95, NA, household17$fuel))
levels(household17$fuel)<-c("Electricity","Natural gas","Biogas", "Coal, lignite","Charcoal","Wood","Straw/shrubs/grass", "Agricultural crop")
levels(household17$fuel)<-c(4,3,3,2,2,1,0,0)
household17$fuel<-as.integer(as.character(household17$fuel))

household17$roof<-as.factor(ifelse(household17$roof==96, NA, household17$roof))
levels(household17$roof)<-c("No roof","Thatch / palm leaf","Sod", "Rustic mat", "Palm / bamboo","Wood planks",
                            "Cardboard", "Metal","Wood","Calamine / cement fiber","Ceramic tiles","Cement","Roofing shingles")
levels(household17$roof)<-c(0,0,0,1,1,1,1,1,2,2,2,2,2)
household17$roof<-as.integer(as.character(household17$roof))

household17$walls<-as.factor(ifelse(household17$walls==96, NA, household17$walls))
levels(household17$walls)<-c("Cane / palm / trunks","Dirt","Bamboo with mud",
                             "Stone with mud","Uncovered adobe", "Plywood","Cardboard","Cement",
                             "Stone with lime / cement",
                             "Bricks","Cement blocks", "Covered adobe", "Wood planks / shingles")
levels(household17$walls)<-c(0,0,0,2,1,2,0,2,2,2,2,2,2)
household17$walls<-as.integer(as.character(household17$walls))

household17$floor<-as.factor(ifelse(household17$floor==96, NA, household17$floor))
levels(household17$floor)<-c("Earth, sand","Dung","Wood planks", "Palm/bamboo","Vinyl, asphalt strips",
                             "Ceramic tiles","Cement","Carpet")
levels(household17$floor)<-c(0,0,1,1,2,2,2,2,2)
household17$floor<-as.integer(as.character(household17$floor))

household17$rooms<-ifelse(household17$rooms==0,1,household17$rooms)

household17$vehicle<-ifelse(household17$car_truck==1,2,
                            ifelse(household17$motorcycle_scooter==1,1,
                                   ifelse(household17$bicycle==1|household17$animal_cart==1,0,NA)))

possess<-c("radio", "tv",
           "fridge", "bicycle", "motorcycle_scooter",
           "car_truck", 
           "landline", 
           "mobile", "clock", "animal_cart", 
           "boat")

household17$possessions<-(rowSums(household17[,c(possess)], na.rm=TRUE))/length(possess)

household17$memroom<-as.numeric(as.character(household17$rooms))/as.numeric(as.character(household17$members))

missing.vars<-which(needed%ni%names(household17))
for(i in 1:length(missing.vars)){
  current_names<-names(household17)
  newvar<-rep(NA, nrow(household17))
  household17<-cbind(household17, newvar)
  colnames(household17)<-c(current_names, needed[missing.vars[i]])
}

household<-merge(household17, geographic17[,c("DHSCLUST", "LATNUM", "LONGNUM")], by="DHSCLUST")

MIS17.coords<-data.frame(coords.x1=household$LONGNUM, coords.x2=household$LATNUM)

MIS17<-household

########################
# Bind into one dataset
########################
merged_wealth<-rbind(WHS[,c(needed)], TARIS[,c(needed)], IHPS13[,c(needed)], IHPS16[,c(needed)],
                     IHS04[,c(needed)], IHS10[,c(needed)], IHS16[,c(needed)], DHS00[,c(needed)], DHS04[,c(needed)], 
                     DHS10[,c(needed)], DHS16[,c(needed)], MICS13[,c(needed)], MIS12[,c(needed)], MIS14[,c(needed)], MIS17[,c(needed)])

merged_wealth$toilet_share<-ifelse(merged_wealth$toilet_share==1,0,1)

merged_coords<-rbind(WHS.coords, TARIS.coords, IHPS13.coords, IHPS16.coords,
                    IHS04.coords, IHS10.coords, IHS16.coords, DHS00.coords, DHS04.coords, 
                    DHS10.coords, DHS16.coords, MICS13.coords, MIS12.coords, MIS14.coords, MIS17.coords)
  
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

save(merged_wealth, file="Data/Confounder_data/Malawi/merged_wealth.csv")

#################
#Factor analysis 
#################
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

#(2) Create correlation matrix
matrices<-list()
for(i in 1:length(newdat_keep)){
  data=newdat_keep[[i]]
  matrix<-mixedCor(data= data, ncat=2)
  matrices<-c(matrices, list(matrix))
}

data<-newdat_keep[[1]]
idx<-which(names(data)=="domestic")
data<-data[,-idx]
newdat_keep[[1]]<-data

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
  model<-lm(wealth_check$comb.score~wealth_check[,idx])
  res<-model$coefficients[2]
  if(res<0){
    miss<-(as.numeric((summary(is.na(wealth_check[,idx]))[3]))/nrow(wealth_check))*100
    message(paste0(variable, ":(miss ",round(miss),"%)"))
    print(res)
    
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
merge<-over(wealth.sp, districts)
merge$districts<-merge$ADM2_EN
wealth.sp@data$district<-ifelse(is.na(wealth.sp@data$district), as.character(merge$districts), as.character(wealth.sp@data$district))
idxNA = which (is.na(wealth.sp@data$district))
missing_district = wealth.sp[idxNA,]
plot(districts)
plot(missing_district, add=TRUE, col="red")

dists <- dist2Line(missing_district,districts)
missing_district@data$district<-as.factor(districts@data$ADM2_EN[dists[,4]])
wealth.sp@data$district<-ifelse(is.na(wealth.sp@data$district), as.character(missing_district@data$district), as.character(wealth.sp@data$district))

wealth.sp@data$district<-as.factor(wealth.sp@data$district)
levels(wealth.sp@data$district)<-c("Balaka","Blantyre", "Blantyre and Blantyre City", "Blantyre City", 
                                "Blantyre City", "Blantyre","Blantyre and Blantyre City", "Chikwawa", "Chiradzulu",
                                "Chitipa", "Chiradzulu", "Dedza","Dowa",
                                "Karonga","Kasungu","Likoma", "Lilongwe",  "Lilongwe and Lilongwe City",   
                                "Lilongwe City", "Lilongwe City", "Lilongwe", "Lilongwe and Lilongwe City", "Machinga", 
                                "Machinga", "Mangochi","Mchinji","Mulanje", 
                                "Mulanje", "Mwanza", "Mzimba", "Mzimba and Mzuzu City", "Mzimba and Mzuzu City",
                                "Mzuzu City", "Nsanje","Neno","Nkhata Bay", 
                                "Nkhata Bay", "Nkhata Bay and Likkoma", "Nkhotakota", "Nkhotakota", 
                                "Nkhotakota","Nsanje","Ntcheu","Ntchisi","Phalombe","Rumphi","Salima","Thyolo",       
                                "Zomba", "Zomba and Zomba City", "Zomba City", "Zomba City", "Zomba and Zomba City")


poverty2<-raster("Data/Confounder_data/Malawi/mwi11povcons200.tif")
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

writeOGR(obj=wealth, dsn="Data/Confounder_data/Malawi/wealth", layer="wealth", driver="ESRI Shapefile", overwrite_layer = TRUE)