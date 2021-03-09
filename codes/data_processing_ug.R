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
UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
districts<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_4")
proj<-proj4string(districts)

#----------------------#
#   Read in data
#----------------------#
#######
# DHS
#######
dhs.data <- read.csv("Data/Exposure_data/Uganda/inputs/DHS/idhs_00016.csv")
dhs.data<-subset(dhs.data, YEAR>2001)
dhs.data<-subset(dhs.data, COUNTRY==800)
dhs.data$district<-as.factor(dhs.data$GEO_UG2006_2016)

dhs.data$HHID<-paste(dhs.data$DHSID, dhs.data$HHNUMALL)
dhs.data$URBANHH<-ifelse(dhs.data$URBANHH==2,0,1)

dhs.data$CATTLENUM<-ifelse(dhs.data$CATTLENUM>95,NA, dhs.data$CATTLENUM)
dhs.data$CBNUM<-ifelse(dhs.data$COWBULLNUM>95,NA, dhs.data$COWBULLNUM)
dhs.data$CATTLENUM<-rowSums(dhs.data[,c("CATTLENUM", "CBNUM")], na.rm=TRUE)
dhs.data$PIGNUM<-ifelse(dhs.data$PIGNUM>95,NA, dhs.data$PIGNUM)

#Collapse on household
dhs.hh<-dhs.data %>% group_by(HHID, YEAR)%>%
  summarise(num.id=n(), DHSCLUST= head(DHSID,1), HHMEMBERS= head(HHMEMBERS,1), urban=head(URBANHH,1),  
            district= head (district,1), cattle.num = mean(CATTLENUM, na.rm=T), 
            pig.num = mean(PIGNUM, na.rm=T), cattle.bin= head(CATTLEYN, 1), pig.bin=head(PIGYN, 1),
            livestock = mean (LIVESTOCKYN, na.rm=T))%>% 
  arrange(HHID)

dhs.hh<-as.data.frame(dhs.hh)

dhs.hh$hh.indicator<-rep(1, nrow(dhs.hh))

#Collapse on cluster
dhs.cluster<-dhs.hh %>% group_by(DHSCLUST, YEAR)%>%
  summarise(cluster.members=sum(HHMEMBERS, na.rm=T), hh.size= mean(HHMEMBERS, na.rm=T), households=n(),
            urban=head(urban,1), district=head(district,1), 
            num.cattle = sum(cattle.num, na.rm=T), 
            cattle.hhs= sum(cattle.bin, na.rm=T), pig.hhs = sum(pig.bin, na.rm=T), 
            num.pigs = sum(pig.num, na.rm=T), 
            livestock.hhs = sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)
dhs.cluster<-as.data.frame(dhs.cluster)

dhs.cluster$clustermembers.cattle<-dhs.cluster$clustermembers.pigs<-dhs.cluster$cluster.members

#Density = # cattle/ # in in house (not #s surveyed!), summed across households
dhs.cluster$cattle.dens<-dhs.cluster$num.cattle/dhs.cluster$clustermembers.cattle
dhs.cluster$pig.dens<-dhs.cluster$num.pigs/dhs.cluster$clustermembers.pigs

DHS06.cluster<-subset(dhs.cluster, YEAR==2006)
DHS11.cluster<-subset(dhs.cluster, YEAR==2011)
DHS16.cluster<-subset(dhs.cluster, YEAR==2016)

DHS06.cluster$survey <-rep("DHS2006", nrow(DHS06.cluster))
DHS11.cluster$survey <-rep("DHS2011", nrow(DHS11.cluster))
DHS16.cluster$survey <-rep("DHS2016", nrow(DHS16.cluster))

###########
# MIS 2009
###########
household09<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2009/UGHR5HDT/UGHR5HFL.DTA", convert.dates = T,convert.factors = F)
individual09<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2009/UGIR5HDT/UGIR5HFL.DTA", convert.dates = T,convert.factors = F)
geographic09<-readOGR("Data/Exposure_data/Uganda/inputs/MIS_2009/UGGE5AFL", layer="UGGE5AFL")
geographic09.mis<-geographic09 #re-use this filename in NPS09

household09$livestock<-ifelse(household09$hv246==1,1,0)
household09$hv246a<-ifelse(household09$hv246a>95, NA, household09$hv246a)
household09$hv246b<-ifelse(household09$hv246b>95, NA, household09$hv246b)
household09$num.cattle<-rowSums(household09[,c("hv246a", "hv246b")], na.rm=TRUE)
household09$num.pigs<-ifelse(household09$hv246g>95, NA, household09$hv246g)

household09$household.members<-household09$hv009
household09$householdmembers.pigs<-ifelse(is.na(household09$num.pigs), NA, household09$household.members)
household09$householdmembers.cattle<-ifelse(is.na(household09$num.cattle), NA, household09$household.members)

household09$district<-as.factor(household09$hv024)
levels(household09$district)<-c("Central 1", "Central 2", "Kampala", "East Central", "Mid Eastern", "North East", 
                                "Mid Northern", "West Nile", "Mid Western", "South")

household09$DHSCLUST<-household09$hv001
household09<-merge(household09, geographic09@data[,c(1,4,15, 16, 17)], by = "DHSCLUST")
household09$urban<-ifelse(household09$URBAN_RURA=="R",0,1)

#Collapse on household
MIS09.hh<-household09 %>% group_by(hhid)%>%
  summarise(num.id=n(), numid.pigs=mean(householdmembers.pigs, na.rm=T),
            numid.cattle=mean(householdmembers.cattle, na.rm=T), DHSCLUST= head(DHSID,1), urban=head(urban,1), num.cattle = mean(num.cattle, na.rm=T), 
            num.pigs = mean(num.pigs, na.rm=T), 
            livestock=mean(livestock, na.rm=T), district=head(district,1))%>% 
  arrange(hhid)
MIS09.hh<-as.data.frame(MIS09.hh)
MIS09.hh$hh.indicator<-rep(1, nrow(MIS09.hh))

#Collapse on cluster
MIS09.cluster<-MIS09.hh %>% group_by(DHSCLUST)%>%
  summarise(ea.members=sum(num.id, na.rm=T), hh.size= mean(num.id, na.rm=T), 
            eamembers.cattle=sum(numid.cattle, na.rm=T), hhsize.cattle= mean(numid.cattle, na.rm=T), 
            eamembers.pigs=sum(numid.pigs, na.rm=T), hhsize.pigs= mean(numid.pigs, na.rm=T), households=sum(hh.indicator, na.rm=T),
            urban=head(urban,1), 
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), district=head(district,1),
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)

MIS09.cluster<-as.data.frame(MIS09.cluster)

#Density = # cattle/ # in household, summed across clusters
MIS09.cluster$cattle.dens<-MIS09.cluster$num.cattle/MIS09.cluster$eamembers.cattle
MIS09.cluster$pig.dens<-MIS09.cluster$num.pigs/MIS09.cluster$eamembers.pigs

MIS09.cluster$survey <-rep("MIS2009", nrow(MIS09.cluster))

###########
# MIS 2014
###########
household14<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2014/UGHR72DT/UGHR72FL.DTA", convert.dates = T,convert.factors = F)
individual14<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2014/UGIR72DT/UGIR72FL.DTA", convert.dates = T,convert.factors = F)
geographic14<-readOGR("Data/Exposure_data/Uganda/inputs/MIS_2014/UGGE71FL", layer="UGGE71FL")

household14$livestock<-ifelse(household14$hv246==1,1,0)
household14$hv246a<-ifelse(household14$hv246a>95, NA, household14$hv246a)
household14$hv246b<-ifelse(household14$hv246b>95, NA, household14$hv246b)
household14$num.cattle<-rowSums(household14[,c("hv246a", "hv246b")], na.rm=TRUE)
household14$num.pigs<-ifelse(household14$hv246g>95, NA, household14$hv246g)

household14$household.members<-household14$hv009
household14$householdmembers.pigs<-ifelse(is.na(household14$num.pigs), NA, household14$household.members)
household14$householdmembers.cattle<-ifelse(is.na(household14$num.cattle), NA, household14$household.members)

household14$district<-as.factor(household14$hv024)
levels(household14$district)<-c("Central 1", "Central 2", "East Central", "Kampala", "Mid Northern","Mid Western", "Mid Eastern", "North East", 
                                "South", "West Nile")

household14$DHSCLUST<-household14$hv001
household14<-merge(household14, geographic14@data[,c(1,4,15, 16, 17)], by = "DHSCLUST")
household14$urban<-ifelse(household14$URBAN_RURA=="R",0,1)

#Collapse on household
MIS14.hh<-household14 %>% group_by(hhid)%>%
  summarise(num.id=mean(household.members, na.rm=T), numid.cattle=mean(householdmembers.cattle, na.rm=T), numid.pigs=mean(householdmembers.pigs, na.rm=T), 
            DHSCLUST= head(DHSID,1), urban=head(urban,1), num.cattle = mean(num.cattle, na.rm=T), 
            num.pigs = mean(num.pigs, na.rm=T), district=head(district,1),
            livestock=mean(livestock, na.rm=T))%>% 
  arrange(hhid)
MIS14.hh<-as.data.frame(MIS14.hh)
MIS14.hh$hh.indicator<-rep(1, nrow(MIS14.hh))

#Collapse on cluster
MIS14.cluster<-MIS14.hh %>% group_by(DHSCLUST)%>%
  summarise(ea.members=sum(num.id, na.rm=T), hh.size= mean(num.id, na.rm=T), households=sum(hh.indicator, na.rm=T),
            eamembers.cattle=sum(numid.cattle, na.rm=T), hhsize.cattle= mean(numid.cattle, na.rm=T), 
            eamembers.pigs=sum(numid.pigs, na.rm=T), hhsize.pigs= mean(numid.pigs, na.rm=T), 
            urban=head(urban,1), district=head(district,1),
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), 
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)

MIS14.cluster<-as.data.frame(MIS14.cluster)

#Density = # cattle/ # in household, summed across clusters
MIS14.cluster$cattle.dens<-MIS14.cluster$num.cattle/MIS14.cluster$eamembers.cattle
MIS14.cluster$pig.dens<-MIS14.cluster$num.pigs/MIS14.cluster$eamembers.pigs

MIS14.cluster$survey <-rep("MIS2014", nrow(MIS14.cluster))

###########
# MIS 2018
###########
household18<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2018/UGHR7IDT/UGHR7IFL.DTA", convert.dates = T,convert.factors = F)
individual18<-read.dta("Data/Exposure_data/Uganda/inputs/MIS_2018/UGIR7IDT/UGIR7IFL.DTA", convert.dates = T,convert.factors = F)
geographic18<-readOGR("Data/Exposure_data/Uganda/inputs/MIS_2018/shapefiles/UGGE7IFL", layer="UGGE7IFL")

household18$livestock<-ifelse(household18$hv246==1,1,0)
household18$hv246a<-ifelse(household18$hv246a>95, NA, household18$hv246a)
household18$hv246b<-ifelse(household18$hv246b>95, NA, household18$hv246b)
household18$num.cattle<-rowSums(household18[,c("hv246a", "hv246b")], na.rm=TRUE)
household18$num.pigs<-ifelse(household18$hv246g>95, NA, household18$hv246g)
household18$household.members<-household18$hv009
household18$householdmembers.pigs<-household18$householdmembers.cattle<-household18$household.members

household18$DHSCLUST<-household18$hv001
household18<-merge(household18, geographic18@data[,c(1,4,13,15)], by = "DHSCLUST")
household18$urban<-ifelse(household18$URBAN_RURA=="R",0,1)
household18$district<-household18$DHSREGNA

#Collapse on household
MIS18.hh<-household18 %>% group_by(hhid)%>%
  summarise(num.id=mean(household.members, na.rm=T), numid.cattle=mean(householdmembers.cattle, na.rm=T), numid.pigs=mean(householdmembers.pigs, na.rm=T), 
            DHSCLUST= head(DHSID,1), urban=head(urban,1), num.cattle = mean(num.cattle, na.rm=T), 
            num.pigs = mean(num.pigs, na.rm=T), district=head(district, 1), 
            livestock=mean(livestock, na.rm=T))%>% 
  arrange(hhid)
MIS18.hh<-as.data.frame(MIS18.hh)
MIS18.hh$hh.indicator<-rep(1, nrow(MIS18.hh))

#Collapse on cluster
MIS18.cluster<-MIS18.hh %>% group_by(DHSCLUST)%>%
  summarise(ea.members=sum(num.id, na.rm=T), hh.size= mean(num.id, na.rm=T), households=sum(hh.indicator, na.rm=T),
            eamembers.cattle=sum(numid.cattle, na.rm=T), hhsize.cattle= mean(numid.cattle, na.rm=T), 
            eamembers.pigs=sum(numid.pigs, na.rm=T), hhsize.pigs= mean(numid.pigs, na.rm=T), 
            urban=head(urban,1), district=head(district,1),
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), 
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)
MIS18.cluster<-as.data.frame(MIS18.cluster)

#Density = # cattle/ # in household, summed across clusters
MIS18.cluster$cattle.dens<-MIS18.cluster$num.cattle/MIS18.cluster$eamembers.cattle
MIS18.cluster$pig.dens<-MIS18.cluster$num.pigs/MIS18.cluster$eamembers.pigs

MIS18.cluster$survey <-rep("MIS2018", nrow(MIS18.cluster))

###########
# AIS 2011
###########
household11<-read.dta("Data/Exposure_data/Uganda/inputs/AIS_2011/UGHR6ADT/UGHR6AFL.DTA", convert.dates = T,convert.factors = F)
individual11<-read.dta("Data/Exposure_data/Uganda/inputs/AIS_2011/UGIR6ADT/UGIR6AFL.DTA", convert.dates = T,convert.factors = F)
geographic11<-readOGR("Data/Exposure_data/Uganda/inputs/AIS_2011/shapefiles/UGGE6AFL", layer="UGGE6AFL")

household11$livestock<-ifelse(household11$hv246==1,1,0)
household11$hv246a<-ifelse(household11$hv246a>95, NA, household11$hv246a)
household11$hv246b<-ifelse(household11$hv246b>95, NA, household11$hv246b)
household11$hv246g<-ifelse(household11$hv246g>95, NA, household11$hv246g)
household11$num.cattle<-rowSums(household11[,c("hv246a", "hv246b", "hv246g")], na.rm=TRUE)
household11$num.pigs<-ifelse(household11$hv246h>95, NA, household11$hv246h)
household11$household.members<-household11$hv009
household11$householdmembers.pigs<-ifelse(is.na(household11$num.pigs), NA, household11$household.members)
household11$householdmembers.cattle<-ifelse(is.na(household11$num.cattle), NA, household11$household.members)

household11$DHSCLUST<-household11$hv001
household11<-merge(household11, geographic11@data[,c(1,4,13,15)], by = "DHSCLUST")
household11$urban<-ifelse(household11$URBAN_RURA=="R",0,1)
household11$district<-household11$DHSREGNA

#Collapse on household
AIS11.hh<-household11 %>% group_by(hhid)%>%
  summarise(num.id=mean(household.members, na.rm=T), numid.cattle=mean(householdmembers.cattle, na.rm=T), numid.pigs=mean(householdmembers.pigs, na.rm=T), 
            DHSCLUST= head(DHSID,1), urban=head(urban,1), num.cattle = mean(num.cattle, na.rm=T), 
            num.pigs = mean(num.pigs, na.rm=T), district=head(district,1),
            livestock=mean(livestock, na.rm=T))%>% 
  arrange(hhid)
AIS11.hh<-as.data.frame(AIS11.hh)
AIS11.hh$hh.indicator<-rep(1, nrow(AIS11.hh))

#Collapse on cluster
AIS11.cluster<-AIS11.hh %>% group_by(DHSCLUST)%>%
  summarise(ea.members=sum(num.id, na.rm=T), hh.size= mean(num.id, na.rm=T), households=sum(hh.indicator, na.rm=T),
            eamembers.cattle=sum(numid.cattle, na.rm=T), hhsize.cattle= mean(numid.cattle, na.rm=T), 
            eamembers.pigs=sum(numid.pigs, na.rm=T), hhsize.pigs= mean(numid.pigs, na.rm=T), 
            urban=head(urban,1), district=head(district,1),
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), 
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)
AIS11.cluster<-as.data.frame(AIS11.cluster)

#Density = # cattle/ # in household, summed across clusters
AIS11.cluster$cattle.dens<-AIS11.cluster$num.cattle/AIS11.cluster$eamembers.cattle
AIS11.cluster$pig.dens<-AIS11.cluster$num.pigs/AIS11.cluster$eamembers.pigs

AIS11.cluster$survey <-rep("AIS2011", nrow(AIS11.cluster))
geographic11.ais<-geographic11

#############################
# National Panel Survey 2009
#############################
household09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC2.dta", convert.dates = T,convert.factors = F)
household09b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/GSEC1.dta", convert.dates = T,convert.factors = F)
livestock09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/AGSEC6A.dta", convert.dates = T,convert.factors = F)
livestock09b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/AGSEC6B.dta", convert.dates = T,convert.factors = F)
geographic09<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2009-2010/UNPS_Geovars_0910.dta", convert.dates = T,convert.factors = F)

livestock09$livestock<-livestock09$a6aq3
idx<-which(is.na(livestock09$livestock)==FALSE)
livestock09<-livestock09[idx,]
livestock09$species<-ifelse(
  livestock09$livestock==1|
    livestock09$livestock==2|
    livestock09$livestock==3|
    livestock09$livestock==4|
    livestock09$livestock==5|
    livestock09$livestock==6,"Cattle", "Other")

livestock09b$small<-livestock09b$a6bq3
idx<-which(is.na(livestock09b$small)==FALSE)
livestock09b<-livestock09b[idx,]
livestock09b$species<-ifelse(
  livestock09b$small==21,"Pigs", "Other")

cattle<-subset(livestock09, species=="Cattle")
cattle.hh<-cattle %>% group_by(HHID)%>%
  summarise(bin.cattle=mean(a6aq4, na.rm=T), num.cattle= sum(a6aq5, na.rm=T))%>% #this is sum over cattle types, not household members
  arrange(HHID)
cattle.hh<-as.data.frame(cattle.hh)
cattle.hh<-cattle.hh[which(!is.na(cattle.hh$num.cattle)),]

pigs<-subset(livestock09b, species=="Pigs")
pigs.hh<-pigs %>% group_by(HHID)%>%
  summarise(bin.pigs=mean(a6bq4), num.pigs= sum(a6bq5, na.rm=T))%>% 
  arrange(HHID)
pigs.hh<-as.data.frame(pigs.hh)
pigs.hh<-pigs.hh[which(!is.na(pigs.hh$num.pigs)),]

nps09.hhv2<-merge(cattle.hh, pigs.hh, by="HHID")

indv.hh<-household09 %>% group_by(HHID)%>%
  summarise(household.members=n())%>% 
  arrange(HHID)%>% as.data.frame()

household09b$district<-as.factor(household09b$region)
levels(household09b$district)<-c("Kampala", "Central", "Eastern", "Northern", "Western")
household09b$CLUSTID<-household09b$comm

NPS09.hh<-merge(nps09.hhv2, indv.hh, by = "HHID")
NPS09.hh<-merge(NPS09.hh, geographic09[,c("HHID", "lat_mod", "lon_mod", "urban")], by = "HHID")
NPS09.hh<-merge(NPS09.hh, household09b[,c("HHID", "district", "CLUSTID")], by = "HHID")

NPS09.hh$hh.indicator<-rep(1, nrow(NPS09.hh))

NPS09.hh$householdmembers.pigs<-ifelse(is.na(NPS09.hh$num.pigs), NA, NPS09.hh$household.members)
NPS09.hh$householdmembers.cattle<-ifelse(is.na(NPS09.hh$num.cattle), NA, NPS09.hh$household.members)

#Collapse on cluster
NPS09.cluster<-NPS09.hh %>% group_by(CLUSTID)%>%
  summarise(ea.members=sum(household.members, na.rm=T), hh.size= mean(household.members, na.rm=T), households=sum(hh.indicator, na.rm=T),
            eamembers.cattle=sum(householdmembers.cattle, na.rm=T), 
            eamembers.pigs=sum(householdmembers.pigs, na.rm=T), 
            urban=head(urban,1), 
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), 
            district=head(district,1), 
            latitude=mean(lat_mod, na.rm=T),
            longitude=mean(lon_mod, na.rm=T))%>% 
  arrange(CLUSTID)%>%as.data.frame()

#Density = # cattle/ # in household, summed across clusters
NPS09.cluster$cattle.dens<-NPS09.cluster$num.cattle/NPS09.cluster$eamembers.cattle
NPS09.cluster$pig.dens<-NPS09.cluster$num.pigs/NPS09.cluster$eamembers.pigs

NPS09.cluster$survey <-rep("NPS2009", nrow(NPS09.cluster))

#############################
# National Panel Survey 2010
#############################
household10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC2.dta", convert.dates = T,convert.factors = F)
household10b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/GSEC1.dta", convert.dates = T,convert.factors = F)
livestock10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/AGSEC6A.dta", convert.dates = T,convert.factors = F)
livestock10b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/AGSEC6B.dta", convert.dates = T,convert.factors = F)
geographic10<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2010-2011/UNPS_Geovars_1011.dta", convert.dates = T,convert.factors = F)

livestock10$livestock<-livestock10$a6aq3
idx<-which(is.na(livestock10$livestock)==FALSE)
livestock10<-livestock10[idx,]
livestock10$species<-ifelse(
  livestock10$livestock==1|
    livestock10$livestock==2|
    livestock10$livestock==3|
    livestock10$livestock==4|
    livestock10$livestock==5|
    livestock10$livestock==6|
    livestock10$livestock==7|
    livestock10$livestock==8|
    livestock10$livestock==9|
    livestock10$livestock==10,"Cattle", "Other")

livestock10b$small<-livestock10b$a6bq3
idx<-which(is.na(livestock10b$small)==FALSE)
livestock10b<-livestock10b[idx,]
livestock10b$species<-ifelse(
  livestock10b$small==21,"Pigs", "Other")

cattle<-subset(livestock10, species=="Cattle")
cattle.hh<-cattle %>% group_by(HHID)%>%
  summarise(bin.cattle=mean(a6aq4, na.rm=T), num.cattle= sum(a6aq5a, na.rm=T))%>% #this is sum over cattle types, not household members
  arrange(HHID)
cattle.hh<-as.data.frame(cattle.hh)

pigs<-subset(livestock10b, species=="Pigs")
pigs.hh<-pigs %>% group_by(HHID)%>%
  summarise(bin.pigs=mean(a6bq4), num.pigs= sum(a6bq5a, na.rm=T))%>% 
  arrange(HHID)
pigs.hh<-as.data.frame(pigs.hh)

nps10.hhv2<-merge(cattle.hh, pigs.hh, by="HHID")

indv.hh<-household10 %>% group_by(HHID)%>%
  summarise(household.members=n())%>% 
  arrange(HHID)%>% as.data.frame()

household10b$district<-as.factor(household10b$region)
levels(household10b$district)<-c("Kampala", "Central", "Eastern", "Northern", "Western")
household10b$CLUSTID<-household10b$comm

NPS10.hh<-merge(nps10.hhv2, indv.hh, by = "HHID")
NPS10.hh<-merge(NPS10.hh, geographic10[,c("HHID", "lat_mod", "lon_mod", "urban")], by = "HHID")
NPS10.hh<-merge(NPS10.hh, household10b[,c("HHID", "district", "CLUSTID")], by = "HHID")

NPS10.hh$hh.indicator<-rep(1, nrow(NPS10.hh))

NPS10.hh$householdmembers.pigs<-ifelse(is.na(NPS10.hh$num.pigs), NA, NPS10.hh$household.members)
NPS10.hh$householdmembers.cattle<-ifelse(is.na(NPS10.hh$num.cattle), NA, NPS10.hh$household.members)

#Collapse on cluster
NPS10.cluster<-NPS10.hh %>% group_by(CLUSTID)%>%
  summarise(ea.members=sum(household.members, na.rm=T), hh.size= mean(household.members, na.rm=T), households=sum(hh.indicator, na.rm=T),
            eamembers.cattle=sum(householdmembers.cattle, na.rm=T), 
            eamembers.pigs=sum(householdmembers.pigs, na.rm=T), 
            urban=head(urban,1), 
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), 
            district=head(district,1), 
            latitude=mean(lat_mod, na.rm=T),
            longitude=mean(lon_mod, na.rm=T))%>% 
  arrange(CLUSTID)%>%as.data.frame()

#Density = # cattle/ # in household, summed across clusters
NPS10.cluster$cattle.dens<-NPS10.cluster$num.cattle/NPS10.cluster$eamembers.cattle
NPS10.cluster$pig.dens<-NPS10.cluster$num.pigs/NPS10.cluster$eamembers.pigs

NPS10.cluster$survey <-rep("NPS2010", nrow(NPS10.cluster))

#############################
# National Panel Survey 2011
#############################
household11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC2.dta", convert.dates = T,convert.factors = F)
household11b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/GSEC1.dta", convert.dates = T,convert.factors = F)
livestock11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/AGSEC6A.dta", convert.dates = T,convert.factors = F)
livestock11b<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/AGSEC6B.dta", convert.dates = T,convert.factors = F)
geographic11<-read.dta("Data/Exposure_data/Uganda/inputs/LSMS_2011-2012/UNPS_Geovars_1112.dta", convert.dates = T,convert.factors = F)

livestock11$livestock<-livestock11$lvstid
idx<-which(is.na(livestock11$livestock)==FALSE)
livestock11<-livestock11[idx,]
livestock11$species<-ifelse(
  livestock11$livestock==1|
    livestock11$livestock==2|
    livestock11$livestock==3|
    livestock11$livestock==4|
    livestock11$livestock==5|
    livestock11$livestock==6|
    livestock11$livestock==7|
    livestock11$livestock==8|
    livestock11$livestock==9|
    livestock11$livestock==10,"Cattle", "Other")

livestock11b$small<-livestock11b$lvstid
idx<-which(is.na(livestock11b$small)==FALSE)
livestock11b<-livestock11b[idx,]
livestock11b$species<-ifelse(
  livestock11b$small==22,"Pigs", "Other")

cattle<-subset(livestock11, species=="Cattle")
cattle.hh<-cattle %>% group_by(HHID)%>%
  summarise(bin.cattle=mean(a6aq2, na.rm=T), num.cattle= sum(a6aq3a, na.rm=T))%>% #this is sum over cattle types, not household members
  arrange(HHID)
cattle.hh<-as.data.frame(cattle.hh)
cattle.hh$num.cattle<-ifelse(is.na(cattle.hh$bin.cattle), NA, cattle.hh$num.cattle)

pigs<-subset(livestock11b, species=="Pigs")
pigs.hh<-pigs %>% group_by(HHID)%>%
  summarise(bin.pigs=mean(a6bq2), num.pigs= sum(a6bq3a, na.rm=T))%>% 
  arrange(HHID)
pigs.hh<-as.data.frame(pigs.hh)
pigs.hh$num.pigs<-ifelse(is.na(pigs.hh$bin.pigs), NA, pigs.hh$num.pigs)

nps11.hhv2<-merge(cattle.hh, pigs.hh, by="HHID")

indv.hh<-household11 %>% group_by(HHID)%>%
  summarise(household.members=n())%>% 
  arrange(HHID)%>% as.data.frame()

household11b$district<-as.factor(household11b$region)
levels(household11b$district)<-c("Central", "Eastern", "Northern", "Western")
household11b$CLUSTID<-household11b$comm

NPS11.hh<-merge(nps11.hhv2, indv.hh, by = "HHID")
NPS11.hh<-merge(NPS11.hh, geographic11[,c("HHID", "lat_mod", "lon_mod", "urban")], by = "HHID")
NPS11.hh<-merge(NPS11.hh, household11b[,c("HHID", "district", "CLUSTID")], by = "HHID")

NPS11.hh$hh.indicator<-rep(1, nrow(NPS11.hh))

NPS11.hh$householdmembers.pigs<-ifelse(is.na(NPS11.hh$num.pigs), NA, NPS11.hh$household.members)
NPS11.hh$householdmembers.cattle<-ifelse(is.na(NPS11.hh$num.cattle), NA, NPS11.hh$household.members)

#Collapse on cluster
NPS11.cluster<-NPS11.hh %>% group_by(CLUSTID)%>%
  summarise(ea.members=sum(household.members, na.rm=T), hh.size= mean(household.members, na.rm=T), households=sum(hh.indicator, na.rm=T),
            eamembers.cattle=sum(householdmembers.cattle, na.rm=T), 
            eamembers.pigs=sum(householdmembers.pigs, na.rm=T), 
            urban=head(urban,1), 
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), 
            district=head(district,1), 
            latitude=mean(lat_mod, na.rm=T),
            longitude=mean(lon_mod, na.rm=T))%>% 
  arrange(CLUSTID)%>%as.data.frame()

#Density = # cattle/ # in household, summed across clusters
NPS11.cluster$cattle.dens<-NPS11.cluster$num.cattle/NPS11.cluster$eamembers.cattle
NPS11.cluster$pig.dens<-NPS11.cluster$num.pigs/NPS11.cluster$eamembers.pigs

NPS11.cluster$survey <-rep("NPS2011", nrow(NPS11.cluster))
NPS11.cluster<-NPS11.cluster[which(!is.na(NPS11.cluster$latitude)),]

############################
# Make column names match up
############################
DHS06.cluster$geotype<-rep("cluster", nrow(DHS06.cluster))
DHS06.cluster$year = 2006
DHS06.cluster$ea.members<-DHS06.cluster$cluster.members
DHS06.cluster$eamembers.cattle<-DHS06.cluster$clustermembers.cattle
DHS06.cluster$eamembers.pigs<-DHS06.cluster$clustermembers.pigs
DHS06<-DHS06.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
DHS06$baseline_urban<-rep(NA, length(nrow(DHS06)))

DHS11.cluster$geotype<-rep("cluster", nrow(DHS11.cluster))
DHS11.cluster$year = 2011
DHS11.cluster$ea.members<-DHS11.cluster$cluster.members
DHS11.cluster$eamembers.cattle<-DHS11.cluster$clustermembers.cattle
DHS11.cluster$eamembers.pigs<-DHS11.cluster$clustermembers.pigs
DHS11<-DHS11.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
DHS11$baseline_urban<-rep(NA, length(nrow(DHS11)))

DHS16.cluster$geotype<-rep("cluster", nrow(DHS16.cluster))
DHS16.cluster$year= 2016
DHS16.cluster$ea.members<-DHS16.cluster$cluster.members
DHS16.cluster$eamembers.cattle<-DHS16.cluster$clustermembers.cattle
DHS16.cluster$eamembers.pigs<-DHS16.cluster$clustermembers.pigs
DHS16<-DHS16.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
DHS16$baseline_urban<-rep(NA, length(nrow(DHS16)))

MIS09.cluster$geotype<-rep("cluster", nrow(MIS09.cluster))
MIS09.cluster$year = 2009
MIS09<-MIS09.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
MIS09$baseline_urban<-rep(NA, length(nrow(MIS09)))

MIS14.cluster$geotype<-rep("cluster", nrow(MIS14.cluster))
MIS14.cluster$year = 2014
MIS14<-MIS14.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
MIS14$baseline_urban<-rep(NA, length(nrow(MIS14)))

MIS18.cluster$geotype<-rep("cluster", nrow(MIS18.cluster))
MIS18.cluster$year = 2018
MIS18<-MIS18.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
MIS18$baseline_urban<-rep(NA, length(nrow(MIS18)))

AIS11.cluster$geotype<-rep("cluster", nrow(AIS11.cluster))
AIS11.cluster$year = 2018
AIS11<-AIS11.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]
AIS11$baseline_urban<-rep(NA, length(nrow(AIS11)))

NPS09.coords <- cbind(NPS09.cluster$longitude, NPS09.cluster$latitude)
NPS09.cluster$geotype<-rep("ea", nrow(NPS09.cluster))
NPS09.cluster$year = 2009
NPS09.cluster<-NPS09.cluster[,c("CLUSTID", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]

NPS10.coords <- cbind(NPS10.cluster$longitude, NPS10.cluster$latitude)
NPS10.cluster$geotype<-rep("ea", nrow(NPS10.cluster))
NPS10.cluster$year = 2010
NPS10.cluster<-NPS10.cluster[,c("CLUSTID", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]

NPS11.coords <- cbind(NPS11.cluster$longitude, NPS11.cluster$latitude)
NPS11.coords<-na.omit(NPS11.coords)
NPS11.cluster$geotype<-rep("ea", nrow(NPS11.cluster))
NPS11.cluster$year = 2011
NPS11.cluster<-NPS11.cluster[,c("CLUSTID", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households")]

##########################
# Shapefiles
##########################
mshapedhs06<-readOGR("Data/Exposure_data/Uganda/inputs/DHS/shapefiles/UGGE53FL_06", layer="UGGE53FL") 
mshapedhs11<-readOGR("Data/Exposure_data/Uganda/inputs/DHS/shapefiles/UGGE61FL_11", layer="UGGE61FL") 
mshapedhs16<-readOGR("Data/Exposure_data/Uganda/inputs/DHS/shapefiles/UGGE7AFL_16", layer="UGGE7AFL") 

geographic<-proj4string(districts)
polygon<- polygon_from_extent(districts,proj4string=geographic)
mcrop06 <- crop(mshapedhs06, polygon)
mcrop11 <- crop(mshapedhs11, polygon)
mcrop16 <- crop(mshapedhs16, polygon)

dhs06<-merge(mcrop06, DHS06, by.x="DHSID", by.y="DHSCLUST")
dhs11<-merge(mcrop11, DHS11, by.x="DHSID", by.y="DHSCLUST")
dhs16<-merge(mcrop16, DHS16, by.x="DHSID", by.y="DHSCLUST")

mis09<-merge(geographic09.mis, MIS09, by.x="DHSID", by.y="DHSCLUST")
mis14<-merge(geographic14, MIS14, by.x="DHSID", by.y="DHSCLUST")
mis18<-merge(geographic18, MIS18, by.x="DHSID", by.y="DHSCLUST")

ais11<-merge(geographic11.ais, AIS11, by.x="DHSID", by.y="DHSCLUST")

proj<-proj4string(dhs06)
nps09.sp<-SpatialPointsDataFrame(coords = NPS09.coords, data = NPS09.cluster,
                                 proj4string = CRS(proj))

nps10.sp<-SpatialPointsDataFrame(coords = NPS10.coords, data = NPS10.cluster,
                                 proj4string = CRS(proj))

nps11.sp<-SpatialPointsDataFrame(coords = NPS11.coords, data = NPS11.cluster,
                                 proj4string = CRS(proj))

library(raster)

#Combine surveys
merged<-bind(dhs06, dhs11, dhs16, mis09, mis14, mis18, ais11, nps09.sp, nps10.sp, nps11.sp)
merged@data$cluster = rownames(merged@data)

##########################
# Add predictors
##########################
proj<-proj4string(merged)

#Districts
merge<-over(merged, districts) 
merge$parish<-merge$NAME_4
merge$subcounty<-merge$NAME_3
merge$district<-merge$NAME_2
merged@data$district<-merge$parish
idxNA = which (is.na(merged@data$district))
missing_district = merged[idxNA,]
missing_district = missing_district[which(missing_district@data$survey!="WHS2003"),]
plot(districts)
plot(missing_district, add=TRUE, col="red")

dists <- dist2Line(missing_district,districts)
missing_district@data$district<-as.factor(districts@data$NAME_4[dists[,4]])
merged@data$district<-ifelse(is.na(merged@data$district), as.character(missing_district@data$district), as.character(merged@data$district))

#Protected areas
protected<-readOGR("Data/Predictor_data/Protected_areas/WDPA_Jan2019-shapefile", layer="WDPA_Jan2019-shapefile-polygons", p4s=proj)
protected.crop<-crop(protected, polygon)
merge<-over(merged, protected.crop)
merge$protected<-ifelse(is.na(merge$WDPAID)==FALSE, 1,0)
merged@data$protected<-merge$protected

#Elevation
gmted2 <- raster("Data/Predictor_data/Elevation/10s000e_20101117_gmted_med075.tif") 
gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
gmted2.crop<-crop(gmted2, polygon) 
gmted6.crop<-crop(gmted6, polygon)
merged.trans <- spTransform(merged, crs(gmted2.crop))
elevation.1<-raster::extract(gmted2.crop, merged.trans)
elevation.2<-raster::extract(gmted6.crop, merged.trans)
merged@data$elevation.1<-elevation.1
merged@data$elevation.2<-elevation.2
merged@data$elevation<-ifelse(is.na(merged@data$elevation.1)==TRUE, merged@data$elevation.2, merged@data$elevation.1)

#Lakes and wetlands
lakes1<-readOGR("Data/Predictor_data/Lakes_wetlands/GLWD-level1", layer="glwd_1", p4s=proj)
lakes2<-readOGR("Data/Predictor_data/Lakes_wetlands/GLWD-level2", layer="glwd_2", p4s=proj)
lakes1.crop<-crop(lakes1, polygon)
lakes2.crop<-crop(lakes2, polygon)
bind.lakes<-bind(lakes1.crop, lakes2.crop)
victoria<-bind.lakes[which(grepl("Victoria", bind.lakes@data$LAKE_NAME)==TRUE),]
merge<-over(merged, bind.lakes)
merge$water.body<-ifelse(is.na(merge$TYPE)==FALSE, 1,0)
merged@data$water.body<-merge$water.body

merged@data$baseline_urban<-as.factor(merged@data$baseline_urban)
merged@data<-merged@data[,c("ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", 
                            "num.pigs", "survey", "year", "baseline_urban", "district", "cluster", "water.body",
                            "protected", "elevation")]
writeOGR(obj=merged, dsn="Data/Exposure_data/Uganda/Created_datasets/merged_data", layer="merged_map", driver="ESRI Shapefile", overwrite_layer = TRUE)

#----------------------#
# Outcome data
#----------------------#
HAT.tg<-read.csv("Data/Outcome_data/Uganda/HAT_ug_Tbg.csv", na.strings="")
HAT.tg<-HAT.tg[which(!is.na(HAT.tg$Latitude)),]#Remove missing location
HAT.tg.sub<-HAT.tg[,c("District", "Parish", "Focus", "Location_name", "Year", "Parasite", "Surveillance_type", "Census", "People_screened", "New_HAT_cases", "P1", "P2", "P_na", "Latitude", "Longitude")]

HAT.tr<-read.csv("Data/Outcome_data/Uganda/HAT_ug_Tbr.csv", na.strings="")
HAT.tr<-HAT.tr[which(!is.na(HAT.tr$Latitude)),]#Remove missing location
HAT.tr.sub<-HAT.tr[,c("District", "Parish", "Focus", "Location_name", "Year", "Parasite", "Surveillance_type", "Census", "People_screened", "New_HAT_cases", "P1", "P2", "P_na", "Latitude", "Longitude")]

HAT<-rbind(HAT.tg.sub, HAT.tr.sub)

HAT_coords<-cbind(HAT$Longitude, HAT$Latitude)

HAT$District<-as.factor(HAT$District)
HAT$Parish<-as.factor(HAT$Parish)

HAT.sp <- SpatialPointsDataFrame(coords = HAT_coords, data = HAT,
                                 proj4string = crs(geographic))

#(1) Protected areas
merge.protected<-over(HAT.sp, protected)
protectedv2<-ifelse(is.na(merge.protected$WDPAID)==FALSE, 1,0)

#(2) Elevation
merge.trans <- spTransform(HAT.sp, crs(gmted2.crop))
gmted2.crop<-crop(gmted2, districts)
gmted6.crop<-crop(gmted6, districts)
elevation.1<-raster::extract(gmted2, merge.trans)
elevation.2<-raster::extract(gmted6, merge.trans)
elevation<-ifelse(is.na(elevation.1)==TRUE, elevation.2, elevation.1)

#(3) Lakes and wetlands
bind.lakes<-bind(lakes1, lakes2)
lakes.crop<-crop(bind.lakes, districts)
merge.proj<-over(HAT.sp, bind.lakes)
water.body<-ifelse(is.na(merge.proj$TYPE)==FALSE, 1,0)

HATsp.data<-HAT.sp@data
HATsp.data$protected<-protectedv2
HATsp.data$elevation<-elevation
HATsp.data$water.body<-water.body

HAT.sp <- SpatialPointsDataFrame(coords = coordinates(HAT.sp), data = HATsp.data,
                                   proj4string = CRS(proj))


writeOGR(obj=HAT.sp, dsn="Data/Outcome_data/Uganda/HAT_Data", layer="HAT_data", driver="ESRI Shapefile", overwrite_layer=TRUE)

#----------------------#
#Random field prediction
#----------------------#
mapExt=extent(districts)
nx=330 
ny=330
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
merge.trans <- spTransform(gridP, crs(gmted3))
elevation.1<-raster::extract(gmted2, merge.trans)
elevation.2<-raster::extract(gmted6, merge.trans)
elevation<-ifelse(is.na(elevation.1)==TRUE, elevation.2, elevation.1)

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

writeOGR(obj=gridP.sp, dsn="Data/Exposure_data/Uganda/Created_datasets/prediction_data", layer="prediction_data", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=protected.crop, dsn="Data/Exposure_data/Uganda/Created_datasets/prediction_data_protected", layer="prediction_data_protected", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeRaster(gmted2.crop, "Data/Exposure_data/Uganda/Created_datasets/prediction_data_elevation1", format="GTiff", overwrite=TRUE)
writeRaster(gmted6.crop, "Data/Exposure_data/Uganda/Created_datasets/prediction_data_elevation2", format="GTiff", overwrite=TRUE)
writeOGR(obj=lakes.crop, dsn="Data/Exposure_data/Uganda/Created_datasets/prediction_data_lakes", layer="prediction_data_lakes", driver="ESRI Shapefile", overwrite_layer=TRUE)

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
  filename<-paste0("Data/denominator data/LandScan/Uganda/LandScan_", year)
  writeRaster(landScancrop2, filename, format="GTiff", overwrite=TRUE)
}