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
libs <- c("rgdal", "maptools", "gridExtra")
lapply(libs, require, character.only = TRUE)

source('codes/my_functions.R')
geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
states<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/gadm36_SSD_shp", layer="gadm36_SSD_1")
counties<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/gadm36_SSD_shp", layer="gadm36_SSD_2")
subcounties<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/gadm36_SSD_shp", layer="gadm36_SSD_3")
proj<-proj4string(states)

all<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/ssd_admbndl_admall_ssnbs_itos_20160114", layer="ssd_admbndl_admALL_ssnbs_itos_20160114")

humdata<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/ssd_admbnda_adm2_imwg_nbs_20180817", layer="ssd_admbnda_adm2_imwg_nbs_20180817")

#----------------------#
#   Read in data
#----------------------#
###########
# Census
###########
if (!require("ipumsr")) stop("Reading IPUMS data into R requires the ipumsr package. It can be installed using the following command: install.packages('ipumsr')")
ddi <- read_ipums_ddi("Data/Exposure_data/South_Sudan/inputs/Census_2008/ipumsi_00004.xml")
data <- read_ipums_micro(ddi)
census<-as.data.frame(data)

census<-dplyr::rename(census, adm1_08=GEO1_SS2008, adm1=GEO1_SS, adm2_08=GEO2_SS2008, adm2=GEO2_SS2008,
                      region=REGNSS, members=PERSONS, stratum=STRATA, own=OWNERSHIP, own_detail=OWNERSHIPD, 
                      urban=SS2008A_URBAN, rooms=SS2008A_BEDRMS, own2=SS2008A_OWNRSHP, water=SS2008A_WATSRC,
                      fuel=SS2008A_FUELCOOK, toilet=SS2008A_TOILET, car=SS2008A_VEHICLE, motorcycle=SS2008A_MOTOCYC,
                      bicycle=SS2008A_BICYCLE, boat=SS2008A_BOAT, cart=SS2008A_ANIMTRAN, tractor=SS2008A_TRACTOR, 
                      tv=SS2008A_TV, radio=SS2008A_RADIO, mobile=SS2008A_MOBILE, landline=SS2008A_PHONE, 
                      computer=SS2008A_COMPUTER, fridge=SS2008A_REFRIDGE, satellite=SS2008A_SATELL, fan=SS2008A_FAN, 
                      aircon=SS2008A_AIRCON, land=SS2008A_LANDOWN, landarea=SS2008A_AREACULT,
                      cattle=SS2008A_CATTLE, pigs=SS2008A_PIGS)

summary(is.na(census$cattle))
summary(is.na(census$pigs))
census$cattle.dens<-census$cattle/census$members
census$pig.dens<-census$pigs/census$members

census.des<-svydesign(ids=~1, weights=~HHWT, data=census) 

density_c1<-svyby(~cattle, by= ~adm1, denominator = ~members, census.des, svyratio)
density_c1<-dplyr::rename(density_c1, c_dens.se="se.cattle/members", c_dens="cattle/members")
density_c2<-svyby(~cattle, by= ~adm2, denominator = ~members, census.des, svyratio)
density_c2<-dplyr::rename(density_c2, c_dens.se="se.cattle/members", c_dens="cattle/members")

density_p1<-svyby(~pigs, by= ~adm1, denominator = ~members, census.des, svyratio)
density_p1<-dplyr::rename(density_p1, p_dens.se="se.pigs/members", p_dens="pigs/members")
density_p2<-svyby(~pigs, by= ~adm2, denominator = ~members, census.des, svyratio)
density_p2<-dplyr::rename(density_p2, p_dens.se="se.pigs/members", p_dens="pigs/members")

county.census<-merge(density_c2, density_p2, by="adm2")
state.census<-merge(density_c1, density_p1, by="adm1")

#Make state shapefile
state.census$survey<-rep("Census", nrow(state.census))
state.census$year<-rep("2008", nrow(state.census))
state.census$state<-as.factor(state.census$adm1)
levels(state.census$state)<-c("Upper Nile", "Jonglei", "Unity", "Warrap", 
                              "Northern Bahr el Ghazal", "Western Bahr el Ghazal", 
                              "Lakes", "Western Equatoria", "Central Equatoria", 
                              "Eastern Equatoria")

data.state<-levels(state.census$state)
shape.state<-humdata@data$ADM1_EN
'%ni%' <- Negate('%in%')

data.state[which(data.state%ni%shape.state)]

humdata.df <- as(humdata, "data.frame")

humdata.df.agg<-humdata.df%>%group_by(ADM1_EN)%>%
  summarise(ADM1_PCODE=head(ADM1_PCODE,1))%>% 
  arrange(ADM1_EN)
humdata.df.agg<-as.data.frame(humdata.df.agg)
rownames(humdata.df.agg)<-humdata.df.agg$ADM1_EN

humdata.union <- unionSpatialPolygons(humdata, humdata@data$ADM1_EN)

humdata.sp <- SpatialPolygonsDataFrame(humdata.union, humdata.df.agg)

state.sp<-merge(humdata.sp, state.census, by.y="state", by.x="ADM1_EN")

state.sp@data$struct<-1:nrow(state.sp@data)
state.sp@data$unstruct<-1:nrow(state.sp@data)
state.sp@data$ld_c<-log(state.sp@data$c_dens)
state.sp@data$ld_p<-log(state.sp@data$p_dens)

#Delta method: var(p)/(p^2)
var_c<-(state.sp@data$c_dens.se^2)/(state.sp@data$c_dens^2)
state.sp@data$prec_c<-1/var_c
var_p<-(state.sp@data$p_dens.se^2)/(state.sp@data$p_dens^2)
state.sp@data$prec_p<-1/var_p

remove<-c("adm1", "adm1.1", "adm1.2", "adm1.3")
rIdx<-which(names(state.sp@data)%in%remove)
state.sp<-state.sp[,-rIdx]

writeOGR(obj=state.sp, dsn="Data/Exposure_data/South_Sudan/Created_datasets/census_data_state", layer="census_data_state", driver="ESRI Shapefile", overwrite_layer = TRUE)

#Make county shapefile
county.census$county<-as.factor(county.census$adm2)
levels(county.census$county)<-c("Renk", "Manyo, Melut", "Fashoda", 
                                "Maban", "Maiwut", "Luakpiny/Nasir", 
                                "Longochuk", "Ulang", "Baliet", "Malakal", 
                                "Panyikang", "Fangak", "Khorflus", "Ayod", 
                                "Duk", "Uror", "Nyirol", "Akobo", "Pibor, Pochalla", 
                                "Twic East", "Bor South", "Pariang", "Mayom, Abiemnhom",
                                "Rubkona", "Guit", "Koch", "Leer", "Mayendit", "Panyijiar", 
                                "Twic, Abyei", "Gogrial West", "Gogrial East", "Tonj North", 
                                "Tonj East", "Tonj South", "Aweil North", "Aweil East", "Aweil South", 
                                "Aweil West", "Aweil Centre", "Raja", "Jur River", "Wau", 
                                "Cueibet, Rumbek North", "Rumbek Centre", "Wulu", "Rumbek East", 
                                "Yirol West", "Yirol East", "Awerial", "Tambura, Nagero", "Nzara", 
                                "Ezo", "Yambio", "Ibba", "Maridi", "Mvolo, Mundri West", "Mundri East", 
                                "Terekeka", "Juba", "Lainya", "Yei", "Morobo", "Kajo-keji", "Torit", 
                                "Lopa/Lafon", "Kapoeta North", "Kapoeta East", "Kapoeta South", "Budi", "Ikotos", "Magwi")

data.county<-levels(county.census$county)
shape.county<-humdata@data$ADM2_EN

data.county[which(data.county%ni%shape.county)]

humdata@data$new_adm2<-rep(NA, nrow(humdata@data))
for(i in 1:nrow(humdata@data)){
  original_adm<-humdata@data$ADM2_EN[i]
  if(original_adm=="Twic"|original_adm=="Abyei Region"){
    id<-which(grepl("Twic, Abyei", data.county))
  }else{
  id<-which(grepl(original_adm, data.county))
  }
  if(length(id)>=1){
  new_county<-as.character(county.census$county[id])
  humdata@data$new_adm2[i]<-new_county
  }
}

humdata@data[which(is.na(humdata@data$new_adm2)),]

humdata.union <- unionSpatialPolygons(humdata, humdata@data$new_adm2)

humdata.df <- as(humdata, "data.frame")

humdata.df.agg<-humdata.df%>%group_by(new_adm2)%>%
  summarise(Shape_area = sum(Shape_Area, na.rm=T), 
            ADM1_EN=head(ADM1_EN,1), ADM1_PCODE=head(ADM1_PCODE,1))%>% 
  arrange(new_adm2)
humdata.df.agg<-as.data.frame(humdata.df.agg)
humdata.df.agg<-humdata.df.agg[!is.na(humdata.df.agg$new_adm2),]
rownames(humdata.df.agg)<-humdata.df.agg$new_adm2

humdata.sp <- SpatialPolygonsDataFrame(humdata.union, humdata.df.agg)

county.sp<-merge(humdata.sp, county.census, by.y="county", by.x="new_adm2")
remove<-c("adm2", "adm2.1", "adm2.2", "adm2.3")
rIdx<-which(names(county.sp@data)%in%remove)
county.sp<-county.sp[,-rIdx]

county.sp@data$struct<-1:nrow(county.sp@data)
county.sp@data$unstruct<-1:nrow(county.sp@data)
county.sp@data$ld_c<-log(county.sp@data$c_dens)
p_dens_o<-ifelse(county.sp@data$p_dens==0,1e-10, county.sp@data$p_dens)
county.sp@data$ld_p<-log(p_dens_o)

p_dens_se_o<-ifelse(county.sp@data$p_dens.se==0,1e-10, county.sp@data$p_dens.se)
var_c<-(county.sp@data$c_dens.se^2)/(county.sp@data$c_dens^2)
county.sp@data$prec_c<-1/var_c
var_p<-(p_dens_se_o^2)/(p_dens_o^2)
county.sp@data$prec_p<-1/var_p

writeOGR(obj=county.sp, dsn="Data/Exposure_data/South_Sudan/Created_datasets/census_data", layer="census_data", driver="ESRI Shapefile", overwrite_layer = TRUE)
