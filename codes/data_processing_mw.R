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
library(dplyr)
source('codes/my_functions.R')
geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
districts<-readOGR("Data/Exposure_data/Malawi/shapefiles/mwi_admbnda_adm2_nso_20181016", layer="mwi_admbnda_adm2_nso_20181016")
proj<-proj4string(districts)

#----------------------#
#   Read in data
#----------------------#
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

WHS$livestock.bin<-ifelse(WHS$q0718==1,1,0) #17 missing

WHS3<-read.dta13("Data/Exposure_data/Malawi/inputs/WHS_2003/MWI_2003_WHS_v01_M_Stata8/WHS-Malawi_F3.dta") 

WHS3$i<-rep(1, nrow(WHS3))
WHS3.hh<-WHS3 %>% group_by(id)%>%
  summarise(household.members=sum(i, na.rm=T))%>% 
  arrange(id)
WHS3.hh<-as.data.frame(WHS3.hh)

WHS<-merge(WHS, WHS3.hh, by="id")

WHS$location<-paste(WHS$lat, WHS$long, sep="_")
WHS$i<-rep(1, nrow(WHS))

WHS.cluster<-WHS %>% group_by(location)%>%
  summarise(cluster.members=sum(household.members, na.rm=T), 
            households = sum(i, na.rm=T),
            livestock.hhs= sum(livestock.bin, na.rm=T), 
            lat = head (lat, 1), 
            lon = head (lon, 1))%>% 
  arrange(location)
WHS.cluster<-as.data.frame(WHS.cluster)
WHS.cluster<-WHS.cluster[,-1]

################################################
# Technology Adoption and Risk Initiation Survey
################################################
matching<-read.dta13("Data/Exposure_data/Malawi/inputs/TARIHS_2006/MTARI_matching.dta", convert.dates = T,convert.factors = F)
qno<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/_qno.csv")
livestock<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/blivesto_puf.csv")
livestock$district<-as.factor(livestock$addrdp_n)
levels(livestock$district)<-c(NA, "Kasungu", "Lilongwe", "Mchinji", "Nkhotakota", "Salima")
household<-read.csv("Data/Exposure_data/Malawi/inputs/TARIHS_2006/bhroster_puf.csv")
household$i<-rep(1, length(nrow(household)))

livestock$bin.livestock<-ifelse(livestock$i1_donke==1|livestock$i1_cows==1|
                                  livestock$i1_oxenb==1|livestock$i1_goats==1|
                                  livestock$i1_pigs==1|livestock$i1_poult==1|
                                  livestock$i1_sheep==1,1,0)
livestock$bin.cattle<-ifelse(livestock$i1_cows==1 | livestock$i1_oxenb==1,1,0) #No missing values
livestock$bin.pigs<-ifelse(livestock$i1_pigs==1,1,0) #No missing values
livestock$num.pigs<-livestock$i2_pigs
livestock$num.pigs<-ifelse(livestock$bin.pigs==0 & is.na(livestock$num.pigs)==TRUE,0,livestock$num.pigs) #there is some missingness here, however
livestock$num.cattle<-rowSums(livestock[,c("i2_cows", "i2_oxen")], na.rm=TRUE)

dim(livestock)#1087
dim(matching)#954

household<-merge(household[,c("qno","addrdp_n","gac_name","i")], matching, by = "qno")

#Collapse on household
household<-household %>% group_by(qno)%>%
  summarise(household.members=sum(i, na.rm=T), district=head(addrdp_n,1), 
            gac_name=head(gac_name.y,1), village=head(vill_nam,1))%>% 
  arrange(qno)
household<-as.data.frame(household)
household$id<-rep(1, nrow(household))
household<-merge(household, livestock[,c("qno","num.cattle","num.pigs", "bin.livestock", "district")], by="qno")

household$hhmembers.cattle<-ifelse(!is.na(household$num.cattle), household$household.members, NA)
household$hhmembers.pigs<-ifelse(!is.na(household$num.pigs), household$household.members, NA)

#Collapse on village
village<-household %>% group_by(village)%>%
  summarise(village.members=sum(household.members, na.rm=T), 
            village.members.cattle= sum(hhmembers.cattle, na.rm=T),
            village.members.pigs= sum(hhmembers.pigs, na.rm=T),
            household.size=mean(household.members, na.rm=T),
            household.size.cattle=mean (hhmembers.cattle, na.rm=T),
            household.size.pigs=mean (hhmembers.pigs, na.rm=T),
            households=sum(id, na.rm=T),
            district=head(district.x,1), 
            gac_name=head(gac_name,1), 
            num.pigs=sum(num.pigs, na.rm=T), 
            num.cattle=sum(num.cattle, na.rm=T), 
            livestock.hhs = sum(bin.livestock, na.rm=T))%>% 
  arrange(village)
village<-as.data.frame(village)
village$pig.dens<-village$num.pigs/village$village.members.pigs
summary(village$pig.dens)
village$cattle.dens<-village$num.cattle/village$village.members.cattle
summary(village$cattle.dens)
dim(village)#340 villages
village$vill_id = rownames(village)

populated.places<-read.csv("Data/Exposure_data/shapefiles/Geonames/mi/populated_places.csv")
populated.places$names<-populated.places$FULL_NAME_RO
pp.coords<-cbind(populated.places$LONG, populated.places$LAT)
populated.places2<-as.data.frame(populated.places[,37])
colnames(populated.places2)<-c("name")
populated.places.sp <- SpatialPointsDataFrame(coords = pp.coords, data = populated.places2,
                                              proj4string = CRS(proj))

populated.places.districts<-over(populated.places.sp, districts)
populated.places.sp@data$DISTRICT<-populated.places.districts$ADM2_EN

villagesgeo<-readOGR("Data/Exposure_data/Malawi/shapefiles/villagesgeo", layer="villagesgeo", p4s=proj)
villagespoints<-readOGR("Data/Exposure_data/Malawi/shapefiles/mw_villagespoints", layer="mw_villagespoints", p4s=proj)

'%ni%' <- Negate('%in%')
match<-which(villagesgeo@data$NAME %ni% populated.places.sp@data$name)
villagesgeo@data$name= villagesgeo@data$NAME
villagespoints@data$name= villagespoints@data$NAME
bind2<-bind(populated.places.sp, villagesgeo, villagespoints)
bind2$location= ifelse(!is.na(bind2$DISTRICT), paste(bind2$name, bind2$DISTRICT, sep="_"), as.character(bind2$name))

village$location<-paste(village$village, village$district, sep="_")
village.gps<-merge(village, bind2, by="location", all.x=TRUE)
village.gps<-village.gps %>% group_by(location)%>%
  summarise(village=head(village,1),
            ea.members=mean(village.members,na.rm=T),
            eamembers.cattle=mean(village.members.cattle,na.rm=T),
            eamembers.pigs=mean(village.members.pigs,na.rm=T),
            district=head(district,1),
            households=mean(households, na.rm=T), 
            livestock.hhs=mean(livestock.hhs, na.rm=T),
            num.pigs=mean(num.pigs, na.rm=T),
            num.cattle=mean(num.cattle, na.rm=T),
            coords.x1=mean(coords.x1, na.rm=T),
            coords.x2=mean(coords.x2, na.rm=T),
            district=head(district,1))%>% 
  arrange(location)
village.gps<-as.data.frame(village.gps)

nomiss<-village.gps[which(!is.na(village.gps$coords.x1)),]
nlevels(as.factor(nomiss$location))#Matched 143 out of 340
missing<-village.gps[which(is.na(village.gps$coords.x1)),]
nlevels(as.factor(missing$location))#Failed to match 197 out of 340

TARIS2006v2<-village.gps[which(is.na(village.gps$coords.x2)==FALSE),] 
TARIS06.coords<-data.frame(coords.x1=TARIS2006v2$coords.x1, coords.x2=TARIS2006v2$coords.x2)

TARIS2006v2$year<-rep(2006, nrow(TARIS2006v2))

urbanNum = getUrbanicity(year=unique(TARIS2006v2$year), loc=TARIS06.coords, regMap=districts)
urbanFac = factor(x = urbanNum,
                  levels = c(0, 1))
#2 urban clusters, rest rural
TARIS2006v2$urban<-urbanFac

###########
# IHPS 2013
###########
ihps13<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/ag_mod_r1_13.dta", convert.dates = T,convert.factors = F)
household13<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_a_filt_13.dta", convert.dates = T,convert.factors = F) #had error regarding factor levels
ihps13<-merge(household13[,c(2:13)], ihps13, by = "y2_hhid")

ihps13$livestock<-ifelse(ihps13$ag_r00==2,0,1) #did you or anyone in your household own any livestock in the last 12 months?
idx<-which(is.na(ihps13$livestock)==FALSE)
ihps13<-ihps13[idx,]
ihps13$species<-ifelse(
  ihps13$ag_r0a==301|
    ihps13$ag_r0a==302|
    ihps13$ag_r0a==303|
    ihps13$ag_r0a==304, "Cattle", 
  ifelse(ihps13$ag_r0a==309, "Pigs", 
         ifelse(ihps13$livestock==0, "None", "Other")))

ihps13$baseline_urban<-ifelse(ihps13$baseline_rural== 1,1,0) #use this for model-based estimation
ihps13$urban<-ifelse(ihps13$reside==1,1,0)

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


#Collapse on household
cattle<-subset(ihps13, species=="Cattle")
cattle$ag_r01<-ifelse(cattle$ag_r01==2,1,0)
cattle.hh<-cattle %>% group_by(y2_hhid)%>%
  summarise(cattle= mean(ag_r01, na.rm=T), 
            num.cattle= sum(ag_r02, na.rm=T))%>% 
  arrange(y2_hhid)
cattle.hh<-as.data.frame(cattle.hh)# no missingness in cattle.num
cattle.hh$cattle<-ifelse(cattle.hh$cattle>0, 1,0)

pigs<-subset(ihps13, species=="Pigs")
pigs$ag_r01<-ifelse(pigs$ag_r01==2,1,0)
pigs.hh<-pigs %>% group_by(y2_hhid)%>%
  summarise(pigs= mean(ag_r01, na.rm=T), 
            num.pigs= sum(ag_r02, na.rm=T))%>% 
  arrange(y2_hhid)
pigs.hh<-as.data.frame(pigs.hh)
pigs.hh$pigs<-ifelse(pigs.hh$pigs>0, 1, 0)# no missingness in pigs.num

ihps13.hhv2<-merge(cattle.hh, pigs.hh, by="y2_hhid")
ihps13.hhv3<-merge(ihps13.hhv2, ihps13, by="y2_hhid", all.y=TRUE)
ihps13.hhv3$num.pigs<-ifelse(is.na(ihps13.hhv3$num.pigs==TRUE), 0, ihps13.hhv3$num.pigs)
ihps13.hhv3$pigs<-ifelse(is.na(ihps13.hhv3$pigs==TRUE), 0, ihps13.hhv3$pigs)
ihps13.hhv3$num.cattle<-ifelse(is.na(ihps13.hhv3$num.cattle==TRUE), 0, ihps13.hhv3$num.cattle)
ihps13.hhv3$cattle<-ifelse(is.na(ihps13.hhv3$cattle==TRUE), 0, ihps13.hhv3$cattle)

ihps13.hh<-ihps13.hhv3 %>% group_by(y2_hhid)%>%
  summarise(case_id=head(case_id, 1), panel_type = head(qx_type,1), ea_id = head (ea_id,1), cattle= mean(cattle, na.rm=T), 
            num.cattle= mean(num.cattle, na.rm=T), pigs= mean(pigs, na.rm=T), 
            num.pigs= mean(num.pigs, na.rm=T), livestock = mean(livestock, na.rm=T), distance = head(dist_to_IHS3location,1), 
            baseline_urban = mean(baseline_urban, na.rm=T), urban=mean(urban,na.rm=T), district=head(district,1))%>% 
  arrange(y2_hhid)
ihps13.hh<-as.data.frame(ihps13.hh)
ihps13.hh$cattle<-ifelse(is.na(ihps13.hh$cattle)==TRUE, 0, ihps13.hh$cattle)
ihps13.hh$pigs<-ifelse(is.na(ihps13.hh$pigs)==TRUE, 0, ihps13.hh$pigs)
ihps13.hh$num.cattle<-ifelse(is.na(ihps13.hh$num.cattle)==TRUE, 0, ihps13.hh$num.cattle)
ihps13.hh$num.pigs<-ifelse(is.na(ihps13.hh$num.pigs)==TRUE, 0, ihps13.hh$num.pigs)

#Number of people per household
ihps13.indv<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_b_13.dta", convert.dates = T,convert.factors = F)
ihps13.indv$i<-rep(1, length(nrow(ihps13.indv)))
indv.hh<-ihps13.indv %>% group_by(y2_hhid)%>%
  summarise(panel_type = head(qx_type,1), hh.members=sum(i, na.rm=T))%>% 
  arrange(y2_hhid)
indv.hh<-as.data.frame(indv.hh)

ihps13.hh<-merge(ihps13.hh, indv.hh[,c(1,3)], by = "y2_hhid")
geo13<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/HouseholdGeovariables_IHPS_13.dta", convert.dates = T,convert.factors = F) #had error regarding factor levels
ihps13.hh<-merge(ihps13.hh, geo13[,c(1, 49, 50)], by = "y2_hhid")

#Number of people per EA
ihps13.hh$i<-rep(1, length(nrow(ihps13.hh)))
moved<-ifelse(household13$hh_a06==1,0,1)#744 households moved, 1245 didn't

ihps13.hh$hhmembers.cattle<-ihps13.hh$hhmembers.pigs<-ihps13.hh$hh.members

ihps13.ea<-ihps13.hh %>% group_by(LON_DD_MOD)%>%
  summarise(panel_type = head(panel_type,1), hh.size=mean(hh.members, na.rm=T), hh.size.cattle=mean(hhmembers.cattle, na.rm=T),
            hh.size.pigs=mean(hhmembers.pigs, na.rm=T), households=sum(i, na.rm=T), 
            ea.members = sum(hh.members, na.rm=T), eamembers.cattle=sum(hhmembers.cattle, na.rm=T),
            eamembers.pigs=sum(hhmembers.pigs, na.rm=T), ea_id = head(ea_id, 1), baseline_urban = mean (baseline_urban, na.rm=T),
            urban = mean (urban, na.rm=T), district= head(district,1), mean_distance = mean (distance, na.rm=T),
            cattle=mean(cattle, na.rm=T), pigs=mean(pigs, na.rm=T), num.cattle=sum(num.cattle, na.rm=T), 
            num.pigs = sum (num.pigs, na.rm=T), Latitude=head(LAT_DD_MOD, 1), Longitude=head(LON_DD_MOD,1), 
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(LON_DD_MOD)
ihps13.ea<-as.data.frame(ihps13.ea)
ihps13.ea$urban<-ifelse(ihps13.ea$urban>0,1,0)#only 1 ea had value 0.875

ihps13.ea$cattle<-ifelse(ihps13.ea$cattle>0, 1,0)
ihps13.ea$pigs<-ifelse(ihps13.ea$pigs>0, 1,0)
ihps13.ea$ea_id_2<-rownames(ihps13.ea)
ihps13.ea<-ihps13.ea[,-1]

ihps13.ea$cattle.dens<-ihps13.ea$num.cattle/ihps13.ea$eamembers.cattle
ihps13.ea$pig.dens<-ihps13.ea$num.pigs/ihps13.ea$eamembers.pigs

ihps13<-ihps13.ea
ihps13$survey <-rep("IHPS2013", nrow(ihps13))

###########
# IHPS 2016
###########
ihps16<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/ag_mod_r1_16.dta", convert.dates = T,convert.factors = F)
household16<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_a_filt_16.dta", convert.dates = T,convert.factors = F) #had error regarding factor levels
scratch<-merge(household16, household13[,c(2,7)], by="y2_hhid")
ihps16<-merge(household16[,c(1,3:14,18)], ihps16, by = "y3_hhid")
ihps16<-merge(scratch[,c(1:2,34)], ihps16, by="y3_hhid")

ihps16$livestock<-ifelse(ihps16$ag_r00==2,0,1) #did you or anyone in your household own any livestock in the last 12 months?
idx<-which(is.na(ihps16$livestock)==FALSE)
ihps16<-ihps16[idx,]
ihps16$species<-ifelse(
  ihps16$ag_r0a==301|
    ihps16$ag_r0a==302|
    ihps16$ag_r0a==303|
    ihps16$ag_r0a==304, "Cattle", 
  ifelse(ihps16$ag_r0a==309, "Pigs", 
         ifelse(ihps16$livestock==0, "None", "Other")))

ihps16$baseline_urban<-ifelse(ihps16$baseline_rural== 1,1,0) #use this for model-based estimation
ihps16$urban<-ifelse(ihps16$reside==1,1,0)

ihps16$district<-ifelse(ihps16$district==101, "Chitipa", 
                        ifelse(ihps16$district==102, "Karonga", 
                               ifelse(ihps16$district==103, "Nkhatabay", 
                                      ifelse(ihps16$district==104, "Rumphi", 
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
                                                                                                                                                                                                                            ifelse(ihps16$district==314, "Zomba City", "Blantyre City"))))))))))))))))))))))))))))))


#Collapse on household
cattle<-subset(ihps16, species=="Cattle")
cattle$ag_r01<-ifelse(cattle$ag_r01==2,1,0)
cattle.hh<-cattle %>% group_by(y3_hhid)%>%
  summarise(cattle= mean(ag_r01, na.rm=T), 
            num.cattle= sum(ag_r02, na.rm=T))%>% 
  arrange(y3_hhid)
cattle.hh<-as.data.frame(cattle.hh)#No missingness in cattle.num
cattle.hh$cattle<-ifelse(cattle.hh$cattle>0, 1,0)

pigs<-subset(ihps16, species=="Pigs")
pigs$ag_r01<-ifelse(pigs$ag_r01==2,1,0)
pigs.hh<-pigs %>% group_by(y3_hhid)%>%
  summarise(pigs= mean(ag_r01, na.rm=T), 
            num.pigs= sum(ag_r02, na.rm=T))%>% 
  arrange(y3_hhid)
pigs.hh<-as.data.frame(pigs.hh)#No missingness in pigs.num
pigs.hh$pigs<-ifelse(pigs.hh$pigs>0, 1, 0)

ihps16.hhv2<-merge(cattle.hh, pigs.hh, by="y3_hhid")
ihps16.hhv3<-merge(ihps16.hhv2, ihps16, by="y3_hhid", all.y=TRUE)
ihps16.hhv3$num.pigs<-ifelse(is.na(ihps16.hhv3$num.pigs==TRUE), 0, ihps16.hhv3$num.pigs)
ihps16.hhv3$pigs<-ifelse(is.na(ihps16.hhv3$pigs==TRUE), 0, ihps16.hhv3$pigs)
ihps16.hhv3$num.cattle<-ifelse(is.na(ihps16.hhv3$num.cattle==TRUE), 0, ihps16.hhv3$num.cattle)
ihps16.hhv3$cattle<-ifelse(is.na(ihps16.hhv3$cattle==TRUE), 0, ihps16.hhv3$cattle)

ihps16.hh<-ihps16.hhv3 %>% group_by(y3_hhid)%>%
  summarise(case_id=head(case_id, 1), HHID=head(HHID,1), ea_id = head (ea_id,1), cattle= mean(cattle, na.rm=T), 
            num.cattle= mean(num.cattle, na.rm=T), pigs= mean(pigs, na.rm=T), 
            num.pigs= mean(num.pigs, na.rm=T),livestock = mean(livestock, na.rm=T), distance_2010IHS3 = head(dist_to_IHS3location,1), 
            distance_2013IHSP3 = head(dist_to_IHPSlocation,1),  moved = mean(mover, na.rm=T),
            baseline_urban = mean(baseline_urban, na.rm=T), urban=mean(urban,na.rm=T), district=head(district,1), 
            y2_hhid = head(y2_hhid,1))%>% 
  arrange(y3_hhid)
ihps16.hh<-as.data.frame(ihps16.hh)
ihps16.hh$cattle<-ifelse(is.na(ihps16.hh$cattle)==TRUE, 0, ihps16.hh$cattle)
ihps16.hh$pigs<-ifelse(is.na(ihps16.hh$pigs)==TRUE, 0, ihps16.hh$pigs)
ihps16.hh$num.cattle<-ifelse(is.na(ihps16.hh$num.cattle)==TRUE, 0, ihps16.hh$num.cattle)
ihps16.hh$num.pigs<-ifelse(is.na(ihps16.hh$num.pigs)==TRUE, 0, ihps16.hh$num.pigs)

#Number of people per household
ihps16.indv<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/hh_mod_b_16.dta", convert.dates = T,convert.factors = F)
ihps16.indv$i<-rep(1, length(nrow(ihps16.indv)))
indv.hh<-ihps16.indv %>% group_by(y3_hhid)%>%
  summarise(panel_type = head(qx_type,1), hh.members=sum(i, na.rm=T))%>% 
  arrange(y3_hhid)
indv.hh<-as.data.frame(indv.hh)

indv.hh$hhmembers.cattle<-indv.hh$hhmembers.pigs<-indv.hh$hh.members

ihps16.hh<-merge(ihps16.hh, indv.hh[,c("y3_hhid", "hh.members", "hhmembers.cattle", "hhmembers.pigs")], by = "y3_hhid")
geo16<-read.dta("Data/Exposure_data/Malawi/inputs/IHPS_2010-2013-2016_long/MWI_2010-2013-2016_IHPS_v02_M_Stata/HouseholdGeovariablesIHPSY3.dta", convert.dates = T,convert.factors = F) #had error regarding factor levels
ihps16.hh<-merge(ihps16.hh, geo16[,c("y3_hhid", "lat_modified", "lon_modified")], by = "y3_hhid")

nlevels(as.factor(ihps16.hh$ea_id))
nlevels(as.factor(geo16$lat_modified))
table(ihps16.hh$moved)

#Number of people per EA
ihps16.hh$i<-rep(1, length(nrow(ihps16.hh)))
nlevels(as.factor(ihps16.hh$lon_modified))#348 after merge

ihps16.ea<-ihps16.hh %>% group_by(lat_modified)%>%
  summarise(hh.size=mean(hh.members, na.rm=T), hhsize.cattle=mean(hhmembers.cattle, na.rm=T),
            hhsize.pigs=mean(hhmembers.pigs, na.rm=T), households=sum(i, na.rm=T), 
            ea.members = sum(hh.members, na.rm=T), eamembers.cattle=sum(hhmembers.cattle, na.rm=T),
            eamembers.pigs=sum(hhmembers.pigs, na.rm=T), ea_id = head(ea_id, 1), baseline_urban = mean (baseline_urban, na.rm=T),
            urban = mean (urban, na.rm=T), district= head(district,1), mean_distance_2010 = mean (distance_2010IHS3, na.rm=T),
            mean_distance_2013 = mean (distance_2013IHSP3, na.rm=T), cattle=mean(cattle, na.rm=T), pigs=mean(pigs, na.rm=T), num.cattle=sum(num.cattle, na.rm=T), 
            num.pigs = sum (num.pigs, na.rm=T), Latitude=head(lat_modified, 1), Longitude=head(lon_modified,1), livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(lat_modified)
ihps16.ea<-as.data.frame(ihps16.ea)
ihps16.ea$urban<-ifelse(ihps16.ea$urban>0,1,0)#only 1 ea had value 0.875

ihps16.ea$cattle<-ifelse(ihps16.ea$cattle>0, 1,0)
ihps16.ea$pigs<-ifelse(ihps16.ea$pigs>0, 1,0)
ihps16.ea$ea_id_2<-rownames(ihps16.ea)
ihps16.ea<-ihps16.ea[,-1]

ihps16.ea$cattle.dens<-ihps16.ea$num.cattle/ihps16.ea$eamembers.cattle
ihps16.ea$pig.dens<-ihps16.ea$num.pigs/ihps16.ea$eamembers.pigs

ihps16<-ihps16.ea
ihps16$survey <-rep("IHPS2016", nrow(ihps16))

ihps16<-subset(ihps16, is.na(ihps16$Latitude)==FALSE)

###########
# IHS 2004
###########
ihs04_u<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/sec_u.dta")
ihs04_a<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/sec_a.dta")

ihs04_u$livestock<-ifelse(ihs04_u$u02==2,0,1) #did you or anyone in your household own any livestock in the last 12 months?
idx<-which(is.na(ihs04_u$livestock)==FALSE)
ihs04_u<-ihs04_u[idx,]
ihs04_u$species<-ifelse(
  ihs04_u$u03==51|
    ihs04_u$u03==52, "Cattle", 
  ifelse(ihs04_u$u03==55, "Pigs",
         ifelse(ihs04_u$livestock==0, "None", "Other")))
cattle<-subset(ihs04_u, species=="Cattle")
cattle$number<-cattle$u04

pigs<-subset(ihs04_u, species=="Pigs")
pigs$number<-pigs$u04

#Collapse on household
cattle.hh<-cattle %>% group_by(case_id)%>%
  summarise(CLUSTERNO = head (psu,1), cattle= mean(number, na.rm=T), 
            num.cattle= sum(number, na.rm=T))%>% 
  arrange(case_id)
cattle.hh<-as.data.frame(cattle.hh)
cattle.hh$num.cattle<-ifelse(is.na(cattle.hh$cattle), NA, cattle.hh$num.cattle)
cattle.hh$cattle.miss <- ifelse(is.na(cattle.hh$num.cattle),1,0)

pigs.hh<-pigs %>% group_by(case_id)%>%
  summarise(pigs= mean(number, na.rm=T), 
            num.pigs= sum(number, na.rm=T), CLUSTERNO=head(psu,1))%>% 
  arrange(case_id)
pigs.hh<-as.data.frame(pigs.hh)
pigs.hh$num.pigs<-ifelse(is.na(pigs.hh$pigs), NA, pigs.hh$num.pigs)
pigs.hh$pigs.miss<-ifelse(is.na(pigs.hh$num.pigs), 1,0)

ihs04.hhv2<-merge(cattle.hh, pigs.hh, by="case_id", all=TRUE)
ihs04.hhv3<-merge(ihs04.hhv2, ihs04_u, by="case_id", all=TRUE)
ihs04.hhv3$pigs.miss<-ifelse(is.na(ihs04.hhv3$pigs.miss),0,ihs04.hhv3$pigs.miss)
ihs04.hhv3$num.pigs<-ifelse(ihs04.hhv3$pigs.miss==1, NA,
                            ifelse(is.na(ihs04.hhv3$num.pigs),0, ihs04.hhv3$num.pigs))

ihs04.hhv3$cattle.miss<-ifelse(is.na(ihs04.hhv3$cattle.miss),0,ihs04.hhv3$cattle.miss)
ihs04.hhv3$num.cattle<-ifelse(ihs04.hhv3$cattle.miss==1, NA,
                              ifelse(is.na(ihs04.hhv3$num.cattle),0, ihs04.hhv3$num.cattle))

ihs04.hh<-ihs04.hhv3 %>% group_by(case_id)%>%
  summarise(CLUSTERNO = head (CLUSTERNO.x,1), num.cattle= mean(num.cattle, na.rm=T), num.pigs=mean(num.pigs, na.rm=T), 
            livestock=mean(livestock, na.rm=T))%>% 
  arrange(case_id)
ihs04.hh<-as.data.frame(ihs04.hh)

#Number of people per household
ihs04.indv<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2004/ihs2_individ.dta")
ihs04.indv$i<-rep(1, length(nrow(ihs04.indv)))
indv.hh<-ihs04.indv %>% group_by(case_id)%>%
  summarise( hh.members=sum(i, na.rm=T))%>% 
  arrange(case_id)
indv.hh<-as.data.frame(indv.hh)

ihs04.hh<-merge(ihs04.hh, indv.hh, by = "case_id")

#Number of people per EA
ihs04.hh$i<-rep(1, length(nrow(ihs04.hh)))
nomiss.p<-which(is.na(ihs04.hh$num.pigs)==FALSE)
nomiss.c<-which(is.na(ihs04.hh$num.cattle)==FALSE)
ihs04.hh.cattle<-ihs04.hh[nomiss.c,]
ihs04.hh.pigs<-ihs04.hh[nomiss.p,]

ihs04.ea<-ihs04.hh %>% group_by(CLUSTERNO)%>%
  summarise(hh.size=mean(hh.members, na.rm=T), households=sum(i, na.rm=T), 
            ea.members = sum(hh.members, na.rm=T), 
            livestock.hhs = sum(livestock, na.rm=T))%>% 
  arrange(CLUSTERNO)

ihs04.ea.cattle<-ihs04.hh.cattle %>% group_by(CLUSTERNO)%>%
  summarise(hh.size.cattle=mean(hh.members, na.rm=T), households.cattle=sum(i, na.rm=T), 
            eamembers.cattle = sum(hh.members, na.rm=T), num.cattle=sum(num.cattle, na.rm=T))%>% 
  arrange(CLUSTERNO)
ihs04.ea.pigs<-ihs04.hh.pigs %>% group_by(CLUSTERNO)%>%
  summarise(hh.size.pigs=mean(hh.members, na.rm=T), households.pigs=sum(i, na.rm=T), 
            eamembers.pigs = sum(hh.members, na.rm=T), num.pigs=sum(num.pigs, na.rm=T))%>% 
  arrange(CLUSTERNO)
ihs04.ea.cattle<-as.data.frame(ihs04.ea.cattle)
ihs04.ea.pigs<-as.data.frame(ihs04.ea.pigs)

ihs04.eav1<-merge(ihs04.ea.pigs, ihs04.ea, by="CLUSTERNO", all.y=TRUE)
ihs04.ea<-merge(ihs04.ea.cattle, ihs04.eav1, by="CLUSTERNO", all.y=TRUE)

ihs04.ea$cattle.dens<-ihs04.ea$num.cattle/ihs04.ea$eamembers.cattle#different denominator for these two b/c different # of households answered #pigs vs. #cattle question, and I'm only using those households for the offset
ihs04.ea$pig.dens<-ihs04.ea$num.pigs/ihs04.ea$eamembers.pigs

#Design variables
ihs04_a$CLUSTERNO<-ihs04_a$psu
ihs04.eav2<-ihs04_a %>% group_by(CLUSTERNO)%>%
  summarise(reside=head(reside,1), district=head(dist,1))%>% 
  arrange(CLUSTERNO)
ihs04.eav2<-as.data.frame(ihs04.eav2)
ihs04.ea<-merge(ihs04.eav2, ihs04.ea, by="CLUSTERNO")
ihs04.ea$urban<-ifelse(ihs04.ea$reside=="Urban",1,0)

ihs04_coords<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2004/2004_IHS_EA.csv")

ihs04.ea$EACODE<-ihs04.ea$CLUSTERNO

ihs04_merge<-merge(ihs04.ea, ihs04_coords, by="EACODE")

ihs04<-ihs04_merge
ihs04$survey <-rep("IHS2004", nrow(ihs04))

###########
# IHS 2010
###########
ihs10<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Agriculture/ag_mod_r1.csv", header=TRUE)

ihs10$livestock<-ifelse(ihs10$ag_r00==2,0,1) 
idx<-which(is.na(ihs10$livestock)==FALSE)
ihs10<-ihs10[idx,]
ihs10$species<-ifelse(
  ihs10$ag_r0a==301|
    ihs10$ag_r0a==302|
    ihs10$ag_r0a==303|
    ihs10$ag_r0a==304, "Cattle", 
  ifelse(ihs10$ag_r0a==309, "Pigs",
         ifelse(ihs10$livestock==0, "None", "Other")))
cattle<-subset(ihs10, species=="Cattle")
cattle$ag_r01<-ifelse(cattle$ag_r01==2,1,0)

pigs<-subset(ihs10, species=="Pigs")
pigs$ag_r01<-ifelse(pigs$ag_r01==2,1,0)

#Collapse on household
cattle.hh<-cattle %>% group_by(case_id)%>%
  summarise(ea_id = head (ea_id,1), cattle= mean(ag_r01, na.rm=T), 
            num.cattle= sum(ag_r02, na.rm=T))%>% 
  arrange(case_id)
cattle.hh<-as.data.frame(cattle.hh)#No missingness
cattle.hh$cattle<-ifelse(cattle.hh$cattle>0,1,0)

pigs.hh<-pigs %>% group_by(case_id)%>%
  summarise(pigs= mean(ag_r01, na.rm=T), 
            num.pigs= sum(ag_r02, na.rm=T), ea_id=head(ea_id,1))%>% 
  arrange(case_id)
pigs.hh<-as.data.frame(pigs.hh)#No missingness again

ihs10.hhv2<-merge(cattle.hh, pigs.hh, by="case_id", all=TRUE)
ihs10.hhv3<-merge(ihs10.hhv2, ihs10, by="case_id", all=TRUE)
ihs10.hhv3$num.pigs<-ifelse(is.na(ihs10.hhv3$num.pigs==TRUE), 0, ihs10.hhv3$num.pigs)
ihs10.hhv3$pigs<-ifelse(is.na(ihs10.hhv3$pigs==TRUE), 0, ihs10.hhv3$pigs)
ihs10.hhv3$num.cattle<-ifelse(is.na(ihs10.hhv3$num.cattle==TRUE), 0, ihs10.hhv3$num.cattle)
ihs10.hhv3$cattle<-ifelse(is.na(ihs10.hhv3$cattle==TRUE), 0, ihs10.hhv3$cattle)

ihs10.hh<-ihs10.hhv3 %>% group_by(case_id)%>%
  summarise(ea_id = head (ea_id,1), cattle= mean(cattle, na.rm=T), 
            num.cattle= mean(num.cattle, na.rm=T),
            pigs= mean (pigs, na.rm=T), num.pigs=mean(num.pigs, na.rm=T), 
            livestock = mean(livestock, na.rm=T))%>% 
  arrange(case_id)
ihs10.hh<-as.data.frame(ihs10.hh)
ihs10.hh$cattle<-ifelse(is.na(ihs10.hh$cattle)==TRUE, 0, ihs10.hh$cattle)
ihs10.hh$pigs<-ifelse(is.na(ihs10.hh$pigs)==TRUE, 0, ihs10.hh$pigs)
ihs10.hh$num.cattle<-ifelse(is.na(ihs10.hh$num.cattle)==TRUE, 0, ihs10.hh$num.cattle)
ihs10.hh$num.pigs<-ifelse(is.na(ihs10.hh$num.pigs)==TRUE, 0, ihs10.hh$num.pigs)

#Number of people per household. 
ihs10.indv<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Household/hh_mod_b.csv")
ihs10.indv$i<-rep(1, length(nrow(ihs10.indv)))
indv.hh<-ihs10.indv  %>% group_by(case_id)%>%
  summarise(panel_type = head(qx_type,1), hh.members=sum(i, na.rm=T))%>% 
  arrange(case_id)
indv.hh<-as.data.frame(indv.hh)

ihs10.hh<-merge(ihs10.hh, indv.hh, by = "case_id")

ihs10.hh$hhmembers.cattle<-ihs10.hh$hhmembers.pigs<-ihs10.hh$hh.members

#Number of people per EA
ihs10.hh$i<-rep(1, length(nrow(ihs10.hh)))
ihs10.ea<-ihs10.hh %>% group_by(ea_id)%>%
  summarise(hh.size=mean(hh.members, na.rm=T), hhsize.cattle=mean(hhmembers.cattle, na.rm=T),
            hhsize.pigs=mean(hhmembers.pigs, na.rm=T), households=sum(i, na.rm=T), 
            ea.members = sum(hh.members, na.rm=T), eamembers.cattle=sum(hhmembers.cattle, na.rm=T),
            eamembers.pigs=sum(hhmembers.pigs, na.rm=T),
            cattle=max(cattle, na.rm=T), pigs=max(pigs, na.rm=T), num.cattle=sum(num.cattle, na.rm=T), 
            num.pigs = sum (num.pigs, na.rm=T), 
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(ea_id)
ihs10.ea<-as.data.frame(ihs10.ea)
ihs10.ea$cattle.dens<-ihs10.ea$num.cattle/ihs10.ea$eamembers.cattle
ihs10.ea$pig.dens<-ihs10.ea$num.pigs/ihs10.ea$eamembers.pigs

#Design variables
ihs10.indv2<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Household/hh_mod_a_filt.csv")
ihs10.eav2<-ihs10.indv2 %>% group_by(ea_id)%>%
  summarise(reside=head(reside,1))%>% 
  arrange(ea_id)
ihs10.eav2<-as.data.frame(ihs10.eav2)
ihs10.ea<-merge(ihs10.eav2, ihs10.ea, by="ea_id")
ihs10.ea$urban<-ifelse(ihs10.ea$reside==1,1,0)

geo10<-read.csv("Data/Exposure_data/Malawi/inputs/IHS_2010/Full_Sample/Geovariables/HH_level/householdgeovariables.csv")
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


design10<-merge(geo10[,c(2:4)], com10[,c(1,22)])

design10.ea<- design10 %>% group_by(ea_id)%>%
  summarise(district=head(district,1), Latitude= head(lat_modified,1), Longitude=head(lon_modified,1))%>% 
  arrange(ea_id)
ihs10.ea<-merge(ihs10.ea, design10.ea, by="ea_id")
ihs10<-ihs10.ea
ihs10$survey <-rep("IHS2010", nrow(ihs10))

###########
# IHS 2016
###########
ihs16<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/AG_MOD_R1.dta", convert.dates = T,convert.factors = F)

ihs16$livestock<-ifelse(ihs16$ag_r00==2,0,1) #did you or anyone in your household own any livestock in the last 12 months?
idx<-which(is.na(ihs16$livestock)==FALSE)
ihs16<-ihs16[idx,]
ihs16$species<-ifelse(
  ihs16$ag_r0a==301|
    ihs16$ag_r0a==302|
    ihs16$ag_r0a==303|
    ihs16$ag_r0a==304, "Cattle", 
  ifelse(ihs16$ag_r0a==309, "Pigs",
         ifelse(ihs16$livestock==0, "None", "Other")))
cattle<-subset(ihs16, species=="Cattle")
cattle$ag_r01<-ifelse(cattle$ag_r01==2,1,0)

pigs<-subset(ihs16, species=="Pigs")
pigs$ag_r01<-ifelse(pigs$ag_r01==2,1,0)

#Collapse on household
cattle.hh<-cattle %>% group_by(case_id)%>%
  summarise(cattle= mean(ag_r01, na.rm=T), 
            num.cattle= sum(ag_r02, na.rm=T))%>% 
  arrange(case_id)
cattle.hh<-as.data.frame(cattle.hh)
cattle.hh$cattle<-ifelse(cattle.hh$cattle>0,1,0)#No missingness

pigs.hh<-pigs %>% group_by(case_id)%>%
  summarise(pigs= mean(ag_r01, na.rm=T), 
            num.pigs= sum(ag_r02, na.rm=T))%>% 
  arrange(case_id)
pigs.hh<-as.data.frame(pigs.hh)#No missingness

ihs16.hhv2<-merge(cattle.hh, pigs.hh, by="case_id")
ihs16.hhv3<-merge(ihs16.hhv2, ihs16, by="case_id", all.y=TRUE)
ihs16.hhv3$num.pigs<-ifelse(is.na(ihs16.hhv3$num.pigs==TRUE), 0, ihs16.hhv3$num.pigs)
ihs16.hhv3$pigs<-ifelse(is.na(ihs16.hhv3$pigs==TRUE), 0, ihs16.hhv3$pigs)
ihs16.hhv3$num.cattle<-ifelse(is.na(ihs16.hhv3$num.cattle==TRUE), 0, ihs16.hhv3$num.cattle)
ihs16.hhv3$cattle<-ifelse(is.na(ihs16.hhv3$cattle==TRUE), 0, ihs16.hhv3$cattle)

ihs16.hh<-ihs16.hhv3 %>% group_by(case_id)%>%
  summarise(cattle= mean(cattle, na.rm=T), 
            num.cattle= mean(num.cattle, na.rm=T),
            pigs= mean (pigs, na.rm=T), num.pigs=mean(num.pigs, na.rm=T), 
            livestock=mean(livestock, na.rm=T))%>% 
  arrange(case_id)
ihs16.hh<-as.data.frame(ihs16.hh)
ihs16.hh$cattle<-ifelse(is.na(ihs16.hh$cattle)==TRUE, 0, ihs16.hh$cattle)
ihs16.hh$pigs<-ifelse(is.na(ihs16.hh$pigs)==TRUE, 0, ihs16.hh$pigs)
ihs16.hh$num.cattle<-ifelse(is.na(ihs16.hh$num.cattle)==TRUE, 0, ihs16.hh$num.cattle)
ihs16.hh$num.pigs<-ifelse(is.na(ihs16.hh$num.pigs)==TRUE, 0, ihs16.hh$num.pigs)

#Number of people per household.
ihs16.indv<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HH_MOD_B.dta", convert.dates = T,convert.factors = F)
ihs16.indv$i<-rep(1, length(nrow(ihs16.indv)))
indv.hh<-ihs16.indv %>% group_by(case_id)%>%
  summarise(hh.members=sum(i, na.rm=T))%>% 
  arrange(case_id)
indv.hh<-as.data.frame(indv.hh)

ihs16.hh<-merge(ihs16.hh, indv.hh, by = "case_id")
ihs16.hh$hhmembers.cattle<-ihs16.hh$hhmembers.pigs<-ihs16.hh$hh.members

#Number of people per EA
ihs16.hh$i<-rep(1, length(nrow(ihs16.hh)))
geo16<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HouseholdGeovariables_stata11/HouseholdGeovariablesIHS4.dta", convert.dates = T,convert.factors = F)

#urban/rural
ihs16.indv2<-read.dta("Data/Exposure_data/Malawi/inputs/IHS_2016/MWI_2016_IHS-IV_v02_M_Stata/HH_MOD_A_FILT.dta", convert.dates = T,convert.factors = F)
ihs16.indv2<-merge(geo16[,c(1,3,4)], ihs16.indv2[,c(1,3,5,6)], by="case_id")

ihs16.hhgeo<-merge(ihs16.hh, ihs16.indv2, by="case_id")
ihs16.hhgeo$urban<-ifelse(ihs16.hhgeo$reside==1,1,0)
ihs16.hhgeo$district<-ifelse(ihs16.hhgeo$district==101, "Chitipa", 
                             ifelse(ihs16.hhgeo$district==102, "Karonga", 
                                    ifelse(ihs16.hhgeo$district==103, "Nkhatabay", 
                                           ifelse(ihs16.hhgeo$district==104, "Rumphi", 
                                                  ifelse(ihs16.hhgeo$district==107, "Mzuzu City", 
                                                         ifelse(ihs16.hhgeo$district==105, "Mzimba", 
                                                                ifelse(ihs16.hhgeo$district==201, "Kasungu", 
                                                                       ifelse(ihs16.hhgeo$district==202, "Nkhota Kota", 
                                                                              ifelse(ihs16.hhgeo$district==203, "Ntchisi", 
                                                                                     ifelse(ihs16.hhgeo$district==204, "Dowa",
                                                                                            ifelse(ihs16.hhgeo$district==205, "Salima", 
                                                                                                   ifelse(ihs16.hhgeo$district==206, "Lilongwe", 
                                                                                                          ifelse(ihs16.hhgeo$district==207, "Mchinji", 
                                                                                                                 ifelse(ihs16.hhgeo$district==208, "Dedza", 
                                                                                                                        ifelse(ihs16.hhgeo$district==209, "Ntcheu", 
                                                                                                                               ifelse(ihs16.hhgeo$district==210, "Lilongwe City", 
                                                                                                                                      ifelse(ihs16.hhgeo$district==301, "Mangochi", 
                                                                                                                                             ifelse(ihs16.hhgeo$district==302, "Machinga", 
                                                                                                                                                    ifelse(ihs16.hhgeo$district==303, "Zomba", 
                                                                                                                                                           ifelse(ihs16.hhgeo$district==304, "Chiradzulu", 
                                                                                                                                                                  ifelse(ihs16.hhgeo$district==305, "Blantyre", 
                                                                                                                                                                         ifelse(ihs16.hhgeo$district==306, "Mwanza", 
                                                                                                                                                                                ifelse(ihs16.hhgeo$district==307, "Thyolo", 
                                                                                                                                                                                       ifelse(ihs16.hhgeo$district==308, "Mulanje", 
                                                                                                                                                                                              ifelse(ihs16.hhgeo$district==309, "Phalombe", 
                                                                                                                                                                                                     ifelse(ihs16.hhgeo$district==310, "Chikwawa", 
                                                                                                                                                                                                            ifelse(ihs16.hhgeo$district==311, "Nsanje", 
                                                                                                                                                                                                                   ifelse(ihs16.hhgeo$district==312, "Balaka", 
                                                                                                                                                                                                                          ifelse(ihs16.hhgeo$district==313, "Neno", 
                                                                                                                                                                                                                                 ifelse(ihs16.hhgeo$district==314, "Zomba City", "Blantyre City"))))))))))))))))))))))))))))))




ihs16.ea<-ihs16.hhgeo %>% group_by(ea_id)%>%
  summarise(hh.size=mean(hh.members, na.rm=T), hhsize.cattle=mean(hhmembers.cattle, na.rm=T),
            hhsize.pigs=mean(hhmembers.pigs, na.rm=T), 
            households=sum(i, na.rm=T), 
            ea.members = sum(hh.members, na.rm=T), eamembers.cattle=sum(hhmembers.cattle, na.rm=T),
            eamembers.pigs=sum(hhmembers.pigs, na.rm=T),
            cattle=max(cattle, na.rm=T), pigs=max(pigs, na.rm=T), num.cattle=sum(num.cattle, na.rm=T), 
            num.pigs = sum (num.pigs, na.rm=T), Latitude=head(lat_modified,1), Longitude=head(lon_modified,1), 
            urban= head(urban,1), district=head(district,1),
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(ea_id)
ihs16.ea<-as.data.frame(ihs16.ea)
ihs16.ea$cattle.dens<-ihs16.ea$num.cattle/ihs16.ea$eamembers.cattle
ihs16.ea$pig.dens<-ihs16.ea$num.pigs/ihs16.ea$eamembers.pigs

ihs16<-ihs16.ea
ihs16$survey <-rep("IHS2016", nrow(ihs16))

############
# DHS 2004
############
dhs.data <- read.csv("Data/Exposure_data/Malawi/inputs/DHS/idhs_00019.csv")

dhs.data$district<-as.factor(dhs.data$GEOALT_MW2010_2016)
levels(dhs.data$district)<-c("Chitipa", "Karonga", 
                             "Nkhatabay and Likoma", "Rumphi", 
                             "Mzimba and Mzuzu City", "Kasungu", 
                             "Nkhotakota", "Ntchisi", "Dowa", 
                             "Salima", "Lilongwe", "Mchinji", 
                             "Dedza", "Ntcheu", "Mangochi", "Machinga", 
                             "Zomba", "Chradzulu", "Blantyre", 
                             "Mwanza", "Thyolo", "Mulange", 
                             "Phalombe", "Chikwawa", "Ndanjge", 
                             "Balaka", "Neno")

DHS04<-subset(dhs.data, YEAR==2004)

DHS04$id.indicator<-rep(1,nrow(DHS04))
DHS04$HHID<-paste(DHS04$DHSID, DHS04$HHNUMALL)
DHS04$URBANHH<-ifelse(DHS04$URBANHH==2,0,1)

DHS04$CATTLENUM<-ifelse(DHS04$CATTLENUM>95,NA, DHS04$CATTLENUM)
DHS04$CBNUM<-ifelse(DHS04$COWBULLNUM>95,NA, DHS04$COWBULLNUM)
DHS04$CATTLENUM<-rowSums(DHS04[,c("CATTLENUM", "CBNUM")], na.rm=TRUE)
DHS04$PIGNUM<-ifelse(DHS04$PIGNUM>95,NA, DHS04$PIGNUM)

#Collapse on household
DHS04.hh<-DHS04 %>% group_by(HHID)%>%
  summarise(num.id=sum(id.indicator, na.rm=T), DHSCLUST= head(DHSID,1), PSUHH=head(PSUHH,1), STRATAHH=head(STRATAHH,1),
            DOMAINHH= head(DOMAINHH,1), HHNUMALL= head(HHNUMALL,1), CLUSTERNOALL=head(CLUSTERNOALL,1), ULTAREAUNITALL=head(ULTAREAUNITALL,1),
            HHWEIGHT= head(HHWEIGHT,1), HHMWEIGHT=head(HHMWEIGHT,1), HHMEMBERS= head(HHMEMBERS,1), urban=head(URBANHH,1), GEO_MW = head(GEO_MW1992_2016,1), 
            GEOALT_MW = head (GEOALT_MW2010_2016,1), cattle.num = mean(CATTLENUM, na.rm=T), district=head(district,1), 
            pig.num = mean(PIGNUM, na.rm=T), 
            livestock = mean (LIVESTOCKYN, na.rm=T))%>% 
  arrange(HHID)
DHS04.hh<-as.data.frame(DHS04.hh)

DHS04.hh$hh.indicator<-rep(1, nrow(DHS04.hh))

#Collapse on cluster
DHS04.cluster<-DHS04.hh %>% group_by(DHSCLUST)%>%
  summarise(cluster.members=sum(HHMEMBERS, na.rm=T), hh.size= mean(HHMEMBERS, na.rm=T), households=sum(hh.indicator, na.rm=T), PSUHH=head(PSUHH,1), STRATAHH=head(STRATAHH,1),
            HHWEIGHT= head(HHWEIGHT,1), HHMWEIGHT=head(HHMWEIGHT,1), urban=head(urban,1), district=head(district,1), 
            num.cattle = sum(cattle.num, na.rm=T), 
            num.pigs = sum(pig.num, na.rm=T), 
            livestock.hhs = sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)

DHS04.cluster<-as.data.frame(DHS04.cluster)
DHS04.cluster$clustermembers.cattle<-DHS04.cluster$clustermembers.pigs<-DHS04.cluster$cluster.members

#Density = # cattle/ # in in house (not #s surveyed!), summed across clusters
DHS04.cluster$cattle.dens<-DHS04.cluster$num.cattle/DHS04.cluster$clustermembers.cattle
DHS04.cluster$pig.dens<-DHS04.cluster$num.pigs/DHS04.cluster$clustermembers.pigs

DHS04.cluster$survey <-rep("DHS2004", nrow(DHS04.cluster))

############
# DHS 2010
############
DHS10<-subset(dhs.data, YEAR==2010)

DHS10$id.indicator<-rep(1,nrow(DHS10))
DHS10$HHID<-paste(DHS10$DHSID, DHS10$HHNUMALL)
DHS10$URBANHH<-ifelse(DHS10$URBANHH==2,0,1)

DHS10$CATTLENUM<-ifelse(DHS10$CATTLENUM>95,NA, DHS10$CATTLENUM)
DHS10$CBNUM<-ifelse(DHS10$COWBULLNUM>95,NA, DHS10$COWBULLNUM)
DHS10$CATTLENUM<-rowSums(DHS10[,c("CATTLENUM", "CBNUM")], na.rm=TRUE)
DHS10$PIGNUM<-ifelse(DHS10$PIGNUM>95,NA, DHS10$PIGNUM)

#Collapse on household
DHS10.hh<-DHS10 %>% group_by(HHID)%>%
  summarise(num.id=sum(id.indicator, na.rm=T), DHSCLUST= head(DHSID,1), PSUHH=head(PSUHH,1), STRATAHH=head(STRATAHH,1),
            DOMAINHH= head(DOMAINHH,1), HHNUMALL= head(HHNUMALL,1), CLUSTERNOALL=head(CLUSTERNOALL,1), ULTAREAUNITALL=head(ULTAREAUNITALL,1),
            HHWEIGHT= head(HHWEIGHT,1), HHMWEIGHT=head(HHMWEIGHT,1), HHMEMBERS= head(HHMEMBERS,1), urban=head(URBANHH,1), GEO_MW = head(GEO_MW1992_2016,1), 
            GEOALT_MW = head (GEOALT_MW2010_2016,1), cattle.num = mean(CATTLENUM, na.rm=T), 
            pig.num = mean(PIGNUM, na.rm=T),livestock = mean (LIVESTOCKYN, na.rm=T), district=head(district,1) )%>% 
  arrange(HHID)
DHS10.hh<-as.data.frame(DHS10.hh)

DHS10.hh$hh.indicator<-rep(1, nrow(DHS10.hh))

#Collapse on cluster
DHS10.cluster<-DHS10.hh %>% group_by(DHSCLUST)%>%
  summarise(cluster.members=sum(HHMEMBERS, na.rm=T), hh.size= mean(HHMEMBERS, na.rm=T), households=sum(hh.indicator, na.rm=T), PSUHH=head(PSUHH,1), STRATAHH=head(STRATAHH,1),
            HHWEIGHT= head(HHWEIGHT,1), HHMWEIGHT=head(HHMWEIGHT,1), urban=head(urban,1), district=head(district,1),
            num.cattle = sum(cattle.num, na.rm=T), 
            num.pigs = sum(pig.num, na.rm=T), 
            livestock.hhs = sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)

DHS10.cluster<-as.data.frame(DHS10.cluster)

DHS10.cluster$clustermembers.cattle<-DHS10.cluster$clustermembers.pigs<-DHS10.cluster$cluster.members

#Density = # cattle/ # in in house (not #s surveyed!), summed across clusters
DHS10.cluster$cattle.dens<-DHS10.cluster$num.cattle/DHS10.cluster$clustermembers.cattle 
DHS10.cluster$pig.dens<-DHS10.cluster$num.pigs/DHS10.cluster$clustermembers.pigs

DHS10.cluster$survey <-rep("DHS2010", nrow(DHS10.cluster))

############
# DHS 2016
############
DHS16<-subset(dhs.data, YEAR==2016)

nlevels(as.factor(DHS16$DHSID))#850; good
DHS16$id.indicator<-rep(1,nrow(DHS16))
DHS16$HHID<-paste(DHS16$DHSID, DHS16$HHNUMALL)
DHS16$URBANHH<-ifelse(DHS16$URBANHH==2,0,1)

DHS16$CATTLENUM<-ifelse(DHS16$CATTLENUM>95,NA, DHS16$CATTLENUM)
DHS16$CBNUM<-ifelse(DHS16$COWBULLNUM>95,NA, DHS16$COWBULLNUM)
DHS16$CATTLENUM<-rowSums(DHS16[,c("CATTLENUM", "CBNUM")], na.rm=TRUE)
DHS16$PIGNUM<-ifelse(DHS16$PIGNUM>95,NA, DHS16$PIGNUM)

#Collapse on household
DHS16.hh<-DHS16 %>% group_by(HHID)%>%
  summarise(num.id=sum(id.indicator, na.rm=T), DHSCLUST= head(DHSID,1), PSUHH=head(PSUHH,1), STRATAHH=head(STRATAHH,1),
            DOMAINHH= head(DOMAINHH,1), HHNUMALL= head(HHNUMALL,1), CLUSTERNOALL=head(CLUSTERNOALL,1), ULTAREAUNITALL=head(ULTAREAUNITALL,1),
            HHWEIGHT= head(HHWEIGHT,1), HHMWEIGHT=head(HHMWEIGHT,1), HHMEMBERS= head(HHMEMBERS,1), urban=head(URBANHH,1), GEO_MW = head(GEO_MW1992_2016,1), 
            GEOALT_MW = head (GEOALT_MW2010_2016,1), cattle.num = mean(CATTLENUM, na.rm=T), district=head(district,1),
            pig.num = mean(PIGNUM, na.rm=T), livestock = mean (LIVESTOCKYN, na.rm=T))%>% 
  arrange(HHID)
DHS16.hh<-as.data.frame(DHS16.hh)

DHS16.hh$hh.indicator<-rep(1, nrow(DHS16.hh))

#Collapse on cluster
DHS16.cluster<-DHS16.hh %>% group_by(DHSCLUST)%>%
  summarise(cluster.members=sum(HHMEMBERS, na.rm=T), hh.size= mean(num.id, na.rm=T), households=sum(HHMEMBERS, na.rm=T), PSUHH=head(PSUHH,1), STRATAHH=head(STRATAHH,1),
            HHWEIGHT= head(HHWEIGHT,1), HHMWEIGHT=head(HHMWEIGHT,1), urban=head(urban,1), 
            num.cattle = sum(cattle.num, na.rm=T), district=head(district,1),
            num.pigs = sum(pig.num, na.rm=T), livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)

DHS16.cluster<-as.data.frame(DHS16.cluster)

DHS16.cluster$clustermembers.cattle<-DHS16.cluster$clustermembers.pigs<-DHS16.cluster$cluster.members

#Density = # cattle/ # household members (not # surveyed), summed across clusters
DHS16.cluster$cattle.dens<-DHS16.cluster$num.cattle/DHS16.cluster$clustermembers.cattle
DHS16.cluster$pig.dens<-DHS16.cluster$num.pigs/DHS16.cluster$clustermembers.pigs

DHS16.cluster$survey <-rep("DHS2016", nrow(DHS16.cluster))

###########
# MIS 2012
###########
household12<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2012/MWHR6HDT/MWHR6HFL.DTA", convert.dates = T,convert.factors = F)
individual12<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2012/MWIR6HDT/MWIR6HFL.DTA", convert.dates = T,convert.factors = F)
geospatial12<-read.csv("Data/Exposure_data/Malawi/inputs/MIS_2012/MWGC6AFL/MWGC6AFL.csv")
geographic12<-readOGR("Data/Exposure_data/Malawi/inputs/MIS_2012/MWGE6AFL", layer="MWGE6AFL") #this just gives EA centroids, N=849

household12$livestock<-ifelse(household12$hv246==1,1,0)
household12$hv246a<-ifelse(household12$hv246a>95, NA, household12$hv246a)
household12$hv246b<-ifelse(household12$hv246b>95, NA, household12$hv246b)
household12$num.cattle<-rowSums(household12[,c("hv246a", "hv246b")], na.rm=TRUE)
household12$num.pigs<-ifelse(household12$sh117b>95, NA, household12$sh117b)
household12$household.members<-household12$hv009
household12$householdmembers.pigs<-ifelse(is.na(household12$num.pigs), NA, household12$household.members)
household12$householdmembers.cattle<-household12$household.members
household12$district<-as.factor(household12$shdistr)
levels(household12$district)<-c("Chitipa", "Karonga", "Nkhatabay", "Rumphi", "Mzimba", "Mzuzu City", 
                                "Kasungu", "Nkhotakota", "Ntchisi", "Dowa", "Salima", "Lilongwe", 
                                "Mchinji", "Dedza", "Ntcheu", "Lilongwe city", "Mangochi", "Manchinga", 
                                "Zomba", "Chiradzulu", "Blantyre", "Thyolo", "Mulanje", "Phalombe", 
                                "Chikwawa", "Nsanje", "Balaka", "Zomba city", "Blantyre city")

household12$DHSCLUST<-household12$hv001
household12<-merge(household12, geographic12@data[,c(1,4,15)], by = "DHSCLUST")
household12$urban<-ifelse(household12$URBAN_RURA=="R",0,1)

#Collapse on household
MIS12.hh<-household12 %>% group_by(hhid)%>%
  summarise(num.id=n(), numid.pigs=mean(householdmembers.pigs, na.rm=T),
            numid.cattle=mean(householdmembers.cattle, na.rm=T), DHSCLUST= head(DHSID,1), urban=head(urban,1), num.cattle = mean(num.cattle, na.rm=T), 
            num.pigs = mean(num.pigs, na.rm=T), 
            livestock=mean(livestock, na.rm=T), district=head(district,1))%>% 
  arrange(hhid)
MIS12.hh<-as.data.frame(MIS12.hh)
MIS12.hh$hh.indicator<-rep(1, nrow(MIS12.hh))

#Collapse on cluster
MIS12.cluster<-MIS12.hh %>% group_by(DHSCLUST)%>%
  summarise(ea.members=sum(num.id, na.rm=T), hh.size= mean(num.id, na.rm=T), 
            eamembers.cattle=sum(numid.cattle, na.rm=T), hhsize.cattle= mean(numid.cattle, na.rm=T), 
            eamembers.pigs=sum(numid.pigs, na.rm=T), hhsize.pigs= mean(numid.pigs, na.rm=T), households=sum(hh.indicator, na.rm=T),
            urban=head(urban,1), 
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), district=head(district,1),
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)

MIS12.cluster<-as.data.frame(MIS12.cluster)

#Density = # cattle/ # in household, summed across clusters
MIS12.cluster$cattle.dens<-MIS12.cluster$num.cattle/MIS12.cluster$eamembers.cattle
MIS12.cluster$pig.dens<-MIS12.cluster$num.pigs/MIS12.cluster$eamembers.pigs

MIS12.cluster$survey <-rep("MIS2012", nrow(MIS12.cluster))

###########
# MIS 2014
###########
household14<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2014/MWHR71DT/MWHR71FL.DTA", convert.dates = T,convert.factors = F)
household.members14<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2014/MWPR71DT/MWPR71FL.DTA", convert.dates = T,convert.factors = F)
individual14<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2014/MWIR71DT/MWIR71FL.DTA", convert.dates = T,convert.factors = F)
geospatial14<-read.csv("Data/Exposure_data/Malawi/inputs/MIS_2014/MWGC71FL/MWGC71FL.csv")
geographic14<-readOGR("Data/Exposure_data/Malawi/inputs/MIS_2014/MWGE71FL", layer="MWGE71FL") 

household14$livestock<-ifelse(household14$hv246==1,1,0)
household14$hv246a<-ifelse(household14$hv246a>95, NA, household14$hv246a)
household14$hv246b<-ifelse(household14$hv246b>95, NA, household14$hv246b)
household14$num.cattle<-rowSums(household14[,c("hv246a", "hv246b")], na.rm=TRUE)
household14$num.pigs<-ifelse(household14$sh117b>95, NA, household12$sh117b)
household14$household.members<-household14$hv009
household14$householdmembers.pigs<-ifelse(is.na(household14$num.pigs), NA, household14$household.members)
household14$householdmembers.cattle<-household14$household.members
household14$district<-as.factor(household14$shdistr)
levels(household14$district)<-c("Chitipa", "Karonga", "Nkhatabay", "Rumphi", "Mzimba", "Mzuzu City", 
                                "Kasungu", "Nkhotakota", "Ntchisi", "Dowa", "Salima", "Lilongwe", 
                                "Mchinji", "Dedza", "Ntcheu", "Lilongwe city", "Mangochi", "Manchinga", 
                                "Zomba", "Chiradzulu", "Blantyre", "Thyolo", "Mulanje", "Phalombe", 
                                "Chikwawa", "Nsanje", "Balaka", "Zomba city", "Blantyre city")

household14$DHSCLUST<-household14$hv001
nlevels(as.factor(household14$DHSCLUST))
nlevels(as.factor(geographic14@data$DHSCLUST))
household14<-merge(household14, geographic14@data[,c(1,4,15)], by = "DHSCLUST")
household14$urban<-ifelse(household14$URBAN_RURA=="R",0,1)

#Collapse on household
MIS14.hh<-household14 %>% group_by(hhid)%>%
  summarise(num.id=n(), numid.cattle=mean(householdmembers.cattle, na.rm=T), numid.pigs=mean(householdmembers.pigs, na.rm=T), 
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
# MIS 2017
###########
household17<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWHR7IDT/MWHR7IFL.DTA", convert.dates = T,convert.factors = F)
household.members17<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWPR7IDT/MWPR7IFL.DTA", convert.dates = T,convert.factors = F)
individual17<-read.dta("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWIR7IDT/MWIR7IFL.DTA", convert.dates = T,convert.factors = F)
geographic17<-readOGR("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWGE7IFL", layer="MWGE7IFL") 
geospatial17<-read.csv("Data/Exposure_data/Malawi/inputs/MIS_2017/MW_2017_MIS_01202020_1748_119657/MWGC7JFL/MWGC7JFL.csv")

household17$livestock<-ifelse(household17$hv246==1,1,0)
household17$hv246a<-ifelse(household17$hv246a>95, NA, household17$hv246a)
household17$hv246b<-ifelse(household17$hv246b>95, NA, household17$hv246b)
household17$num.cattle<-rowSums(household17[,c("hv246a", "hv246b")], na.rm=TRUE)
household17$num.pigs<-ifelse(household17$hv246g>95, NA, household17$hv246g)
household17$household.members<-household17$hv009
household17$householdmembers.pigs<-household17$householdmembers.cattle<-household17$household.members

household17$DHSCLUST<-household17$hv001
nlevels(as.factor(household17$DHSCLUST))
nlevels(as.factor(geographic17@data$DHSCLUST))
household17<-merge(household17, geographic17@data[,c(1,4,15)], by = "DHSCLUST")
household17$urban<-ifelse(household17$URBAN_RURA=="R",0,1)

#Collapse on household
MIS17.hh<-household17 %>% group_by(hhid)%>%
  summarise(num.id=n(), numid.cattle=mean(householdmembers.cattle, na.rm=T), numid.pigs=mean(householdmembers.pigs, na.rm=T), 
            DHSCLUST= head(DHSID,1), urban=head(urban,1), num.cattle = mean(num.cattle, na.rm=T), 
            num.pigs = mean(num.pigs, na.rm=T), 
            livestock=mean(livestock, na.rm=T))%>% 
  arrange(hhid)
MIS17.hh<-as.data.frame(MIS17.hh)
MIS17.hh$hh.indicator<-rep(1, nrow(MIS17.hh))

#Collapse on cluster
MIS17.cluster<-MIS17.hh %>% group_by(DHSCLUST)%>%
  summarise(ea.members=sum(num.id, na.rm=T), hh.size= mean(num.id, na.rm=T), households=sum(hh.indicator, na.rm=T),
            eamembers.cattle=sum(numid.cattle, na.rm=T), hhsize.cattle= mean(numid.cattle, na.rm=T), 
            eamembers.pigs=sum(numid.pigs, na.rm=T), hhsize.pigs= mean(numid.pigs, na.rm=T), 
            urban=head(urban,1), 
            num.cattle = sum(num.cattle, na.rm=T), 
            num.pigs = sum(num.pigs, na.rm=T), 
            livestock.hhs=sum(livestock, na.rm=T))%>% 
  arrange(DHSCLUST)

MIS17.cluster<-as.data.frame(MIS17.cluster)

#Density = # cattle/ # in household, summed across clusters
MIS17.cluster$cattle.dens<-MIS17.cluster$num.cattle/MIS17.cluster$eamembers.cattle
MIS17.cluster$pig.dens<-MIS17.cluster$num.pigs/MIS17.cluster$eamembers.pigs

MIS17.cluster$survey <-rep("MIS2017", nrow(MIS17.cluster))

###########
# MICS 5
###########
mics5<-read.spss("Data/Exposure_data/Malawi/inputs/MICS5/Malawi MICS 2013-14 SPSS Datasets/hh.sav", to.data.frame=TRUE)
mics5<-dplyr::rename(mics5, cluster_num= HH1, hh_num= HH2, livestockbin=HC13, cattle=HC14A, equid=HC14B, goats=HC14C, sheep=HC14D, chicken=HC14E, pigs=HC14F, other=HC14G) 
mics5[mics5 %in% c("<NA>")] <- NA
mics5<-mics5[which(mics5$HH9=="Completed"),]

mics5geo<-read.csv("Data/Exposure_data/Malawi/inputs/MICS5/Malawi MICS 2013-14 SPSS Datasets/2013_14_MICS_EA.csv")

mics5<-merge(mics5, mics5geo, by.x="cluster_num", by.y="CLUSTERNO", all.x=TRUE)

summary(is.na(mics5$Lat))#Only 25 observations missing

summary(is.na(mics5$livestockbin))#no missing

mics5$livestock<-ifelse(mics5$livestockbin=="Yes",1,
                        ifelse(mics5$livestockbin=="No",0, NA))

mics5$urban<-ifelse(mics5$HH6=="Rural", 0,1)
mics5$num.cattle<-ifelse(mics5$cattle=="95+", 95,
                         ifelse(mics5$cattle=="DK"|mics5$cattle=="Missing", NA,
                                as.numeric(as.character(mics5$cattle))))
mics5$num.pigs<-ifelse(mics5$pigs=="95+", 95,
                       ifelse(mics5$pigs=="DK"|mics5$pigs=="Missing", NA,
                              as.numeric(as.character(mics5$pigs))))
mics5$hh.members<-mics5$HH11
mics5$district<-mics5$DISTNAME

#Set hh members to NA for hh with missing cattle, missing pigs
mics5$hhmembers.cattle<-ifelse(!is.na(mics5$num.cattle), mics5$hh.members, NA)
mics5$hhmembers.pigs<-ifelse(!is.na(mics5$num.pigs), mics5$hh.members, NA)
mics5$i <- rep(1, nrow(mics5))

mics5.ea<-mics5 %>% group_by(EACODE)%>%
  summarise(ea.members=sum(hh.members, na.rm=T), eamembers.cattle=sum(hhmembers.cattle, na.rm=T),
            eamembers.pigs=sum(hhmembers.pigs, na.rm=T), 
            num.cattle=sum(num.cattle,na.rm=T), num.pigs = sum(num.pigs, na.rm=T),
            hh.size = mean (hh.members, na.rm=T), 
            hhsize.cattle = mean (hhmembers.cattle, na.rm=T), 
            hhsize.pigs = mean (hhmembers.pigs, na.rm=T), 
            Lat = head(Lat, 1), Long=head(Long,1),
            district=head(district,1),
            urban= head(urban,1), 
            cluster_num=head(cluster_num,1), 
            livestock.hhs = sum(livestock, na.rm=T), 
            households=sum(i, na.rm=T))%>% 
  arrange(EACODE)
mics5.ea<-as.data.frame(mics5.ea)

mics5.ea$cattle.dens<-mics5.ea$num.cattle/mics5.ea$eamembers.cattle
mics5.ea$pig.dens<-mics5.ea$num.pigs/mics5.ea$eamembers.pigs

mics5.ea$survey <-rep("MICS2013", nrow(mics5.ea))
mics5<-mics5.ea
nomiss.loc<-which(!is.na(mics5$Lat))
mics5<-mics5[nomiss.loc,]

############################
# Make column names match up
############################
WHS.cluster$geotype<-rep("cluster", nrow(WHS.cluster))
WHS.cluster$year = 2003
WHS.cluster$ea.members<-WHS.cluster$cluster.members
WHS.cluster$eamembers.cattle<-WHS.cluster$cluster.members
WHS.cluster$eamembers.pigs<-WHS.cluster$cluster.members
WHS.cluster$cluster<-1:nrow(WHS.cluster)
WHS.cluster$survey<-rep("WHS2003", nrow(WHS.cluster))
WHS<-WHS.cluster[,c("cluster", "ea.members", "eamembers.pigs", "eamembers.cattle", "households", "livestock.hhs", "survey", "year")]
WHS03.coords<-data.frame(coords.x1=WHS.cluster$lon, coords.x2=WHS.cluster$lat)

DHS04.cluster$geotype<-rep("cluster", nrow(DHS04.cluster))
DHS04.cluster$year = 2004
DHS04.cluster$ea.members<-DHS04.cluster$cluster.members
DHS04.cluster$eamembers.cattle<-DHS04.cluster$clustermembers.cattle
DHS04.cluster$eamembers.pigs<-DHS04.cluster$clustermembers.pigs
DHS04<-DHS04.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "district", "urban", "num.cattle", "num.pigs", "survey", "year", "households", "livestock.hhs")]
DHS04$baseline_urban<-rep(NA, length(nrow(DHS04)))

DHS10.cluster$geotype<-rep("cluster", nrow(DHS10.cluster))
DHS10.cluster$year = 2010
DHS10.cluster$ea.members<-DHS10.cluster$cluster.members
DHS10.cluster$eamembers.cattle<-DHS10.cluster$clustermembers.cattle
DHS10.cluster$eamembers.pigs<-DHS10.cluster$clustermembers.pigs
DHS10<-DHS10.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs","district", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "households", "livestock.hhs")]
DHS10$baseline_urban<-rep(NA, length(nrow(DHS10)))

DHS16.cluster$geotype<-rep("cluster", nrow(DHS16.cluster))
DHS16.cluster$year= 2016
DHS16.cluster$ea.members<-DHS16.cluster$cluster.members
DHS16.cluster$eamembers.cattle<-DHS16.cluster$clustermembers.cattle
DHS16.cluster$eamembers.pigs<-DHS16.cluster$clustermembers.pigs
DHS16<-DHS16.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs","district", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "households", "livestock.hhs")]
DHS16$baseline_urban<-rep(NA, length(nrow(DHS16)))

MIS12.cluster$geotype<-rep("cluster", nrow(MIS12.cluster))
MIS12.cluster$year = 2012
MIS12<-MIS12.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "district","eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "households", "livestock.hhs")]
MIS12$baseline_urban<-rep(NA, length(nrow(MIS12)))

MIS14.cluster$geotype<-rep("cluster", nrow(MIS14.cluster))
MIS14.cluster$year = 2014
MIS14<-MIS14.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "district", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "households", "livestock.hhs")]
MIS14$baseline_urban<-rep(NA, length(nrow(MIS14)))

MIS17.cluster$geotype<-rep("cluster", nrow(MIS17.cluster))
MIS17.cluster$year = 2017
MIS17<-MIS17.cluster[,c("DHSCLUST", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "households", "livestock.hhs")]
MIS17$baseline_urban<-rep(NA, length(nrow(MIS17)))

TARIS06.coords<-data.frame(coords.x1=TARIS2006v2$coords.x1, coords.x2=TARIS2006v2$coords.x2)
TARIS2006v2$geotype<-rep("village", nrow(TARIS2006v2))
TARIS2006v2$year=2006
TARIS2006v2$survey="TARIS2006"
TARIS06<-TARIS2006v2[,c("village", "ea.members", "eamembers.pigs", "eamembers.cattle", "num.cattle", "num.pigs", "survey", "year", "district", "urban", "households", "livestock.hhs")]

ihps13.coords <- cbind(ihps13$Longitude, ihps13$Latitude)
ihps13$geotype<-rep("ea", nrow(ihps13))
ihps13$year = 2013
ihps13v2<-ihps13[,c("ea_id", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households", "livestock.hhs")]

ihps16.coords <- cbind(ihps16$Longitude, ihps16$Latitude)
ihps16$geotype<-rep("ea", nrow(ihps16))
ihps16$year = 2016
ihps16v2<-ihps16[,c("ea_id", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households", "livestock.hhs")]

ihs04.coords<-cbind(ihs04$Long, ihs04$Lat)
ihs04$geotype<-rep("ea", nrow(ihs04))
ihs04$district<-ihs04$DISTRICT
ihs04$year = 2004
ihs04v2<-ihs04[,c("EACODE", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households", "livestock.hhs")]

ihs10.coords <- cbind(ihs10$Longitude, ihs10$Latitude)
ihs10$geotype<-rep("ea", nrow(ihs10))
ihs10$year = 2010
ihs10v2<-ihs10[,c("ea_id", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "households", "district", "livestock.hhs")]

ihs16.coords <- cbind(ihs16$Longitude, ihs16$Latitude)
ihs16$geotype<-rep("ea", nrow(ihs16))
ihs16$year = 2016
ihs16v2<-ihs16[,c("ea_id", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households", "livestock.hhs")]

mics5.coords<-cbind(mics5$Long, mics5$Lat)
mics5$geotype<-rep("cluster", nrow(mics5))
mics5$year = 2013
mics5<-mics5[,c("EACODE", "ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", "num.pigs", "survey", "year", "district", "households", "livestock.hhs")]

##########################
# Shapefiles
##########################
mshapedhs04<-readOGR("Data/Exposure_data/Malawi/inputs/DHS/shapefiles/MWGE4BFL", layer="MWGE4BFL") #this just gives EA centroids, N=521
mshapedhs10<-readOGR("Data/Exposure_data/Malawi/inputs/DHS/shapefiles/MWGE62FL", layer="MWGE62FL") #this just gives EA centroids, N=849
mshapedhs16<-readOGR("Data/Exposure_data/Malawi/inputs/DHS/shapefiles/MWGE7AFL", layer="MWGE7AFL") #this just gives EA centroids, N=850
ihs04shape<-readOGR("Data/Exposure_data/Malawi/inputs/IHS_2004/2004_IHS_EA_Shapefile", layer="2004_IHS_EA_Shapefile") #this just gives EA centroids, N=1102
mis5shape<-readOGR("Data/Exposure_data/Malawi/inputs/MICS5/Malawi MICS 2013-14 SPSS Datasets/2013_14_MICS_Shapefile", layer="2013_14_MICS_Shapefile") #this just gives EA centroids, N=849

eas<-readOGR("Data/Exposure_data/Malawi/shapefiles/eas_bnd", layer="eas_bnd")#make sure this is the correct EA map. Blank part = Lake Malawi
lake<-readOGR("Data/Exposure_data/Malawi/shapefiles/malawi_lake", layer="malawi_lake")
geographic<-proj4string(eas)
polygon<- polygon_from_extent(eas,proj4string=geographic)
mcrop04 <- crop(mshapedhs04, polygon)
mcrop10 <- crop(mshapedhs10, polygon)
mcrop16 <- crop(mshapedhs16, polygon)

dhs04<-merge(mcrop04, DHS04, by.x="DHSID", by.y="DHSCLUST")#collapses on cluster
dhs10<-merge(mcrop10, DHS10, by.x="DHSID", by.y="DHSCLUST")#collapses on cluster
DHS16$ID<-paste("MW2015", substr(DHS16.cluster$DHSCLUST, 7, 15), sep="")
dhs16<-merge(mcrop16, DHS16, by.x="DHSID", by.y="ID")

mis12<-merge(geographic12, MIS12, by.x="DHSID", by.y="DHSCLUST")
mis14<-merge(geographic14, MIS14, by.x="DHSID", by.y="DHSCLUST")
mis17<-merge(geographic17, MIS17, by.x="DHSID", by.y="DHSCLUST")

proj<-proj4string(dhs10)
whs03.sp<-SpatialPointsDataFrame(coords = WHS03.coords, data = WHS,
                                 proj4string = CRS(proj))

ihps13.sp <- SpatialPointsDataFrame(coords = ihps13.coords, data = ihps13v2,
                                    proj4string = CRS(proj))

ihps16.sp <- SpatialPointsDataFrame(coords = ihps16.coords, data = ihps16v2,
                                    proj4string = CRS(proj))

ihs04.sp <- SpatialPointsDataFrame(coords = ihs04.coords, data = ihs04v2,
                                   proj4string = CRS(proj))

ihs10.sp <- SpatialPointsDataFrame(coords = ihs10.coords, data = ihs10v2,
                                   proj4string = CRS(proj))

ihs16.sp <- SpatialPointsDataFrame(coords = ihs16.coords, data = ihs16v2,
                                   proj4string = CRS(proj))

mics5.sp <- SpatialPointsDataFrame(coords = mics5.coords, data = mics5,
                                   proj4string = CRS(proj))


TARIS06.sp <- SpatialPointsDataFrame(coords = TARIS06.coords, data = TARIS06,
                                      proj4string = CRS(proj))

library(raster)

#Combine surveys
merged<-bind(dhs04, dhs10 , dhs16, ihps13.sp, ihps16.sp, ihs04.sp, ihs10.sp, ihs16.sp, mis12, mis14, mis17, mics5.sp, TARIS06.sp, whs03.sp)
merged@data$cluster = rownames(merged@data)

merged@data$district<-as.factor(merged@data$district)
levels(merged@data$district)<-c("Balaka","Blantyre", "Blantyre City", "Blantyre City", 
                                "Blantyre City", "Blantyre","Chikwawa", "Chiradzulu",
                                "Chitipa", "Chiradzulu", "Dedza","Dowa",
                                "Karonga","Kasungu","Lilongwe",  "Lilongwe City",   
                                 "Lilongwe City",  "Lilongwe", "Machinga", 
                                "Machinga", "Mangochi","Mchinji","Mulanje", 
                                "Mulanje", "Mwanza", "Mzimba", "Mzimba and Mzuzu City", 
                                "Mzuzu City", "Nsanje","Neno","Nkhata Bay", 
                                "Nkhata Bay", "Nkhata Bay","Nkhotakota", "Nkhotakota", 
                                "Nsanje","Ntcheu","Ntchisi","Phalombe","Rumphi","Salima","Thyolo",       
                                 "Zomba", "Zomba City", "Zomba City")

##########################
# Add predictors
##########################
proj<-proj4string(merged)

#Districts
merge<-over(merged, districts)
merge$districts<-merge$ADM2_EN
merged@data$district<-ifelse(is.na(merged@data$district), as.character(merge$ADM2_EN), as.character(merged@data$district))#Still missing 267: 15 from DHS 2004, rest from WHS 2003
idxNA = which (is.na(merged@data$district))
missing_district = merged[idxNA,]
missing_district = missing_district[which(missing_district@data$survey!="WHS2003"),]
plot(districts)
plot(missing_district, add=TRUE, col="red")

dists <- dist2Line(missing_district,districts)
missing_district@data$district<-as.factor(districts@data$ADM2_EN[dists[,4]])
merged<-merged[which(merged@data$survey!="WHS2003"),]
merged@data$district<-ifelse(is.na(merged@data$district), as.character(missing_district@data$district), as.character(merged@data$district))

#Cities
#the four major cities are Blantyre, Lilongwe, Mzuzu, and Zomba. These are separate layers (e.g. Lilongwe vs. Lilongwe City)
table(merged@data$district, merged@data$urban)

#Protected areas
protected<-readOGR("Data/Predictor_data/Protected_areas/WDPA_Jan2019-shapefile", layer="WDPA_Jan2019-shapefile-polygons", p4s=proj)
protected.crop<-crop(protected, polygon)
merge<-over(merged, protected.crop)
merge$protected<-ifelse(is.na(merge$WDPAID)==FALSE, 1,0)
merged@data$protected<-merge$protected

#Elevation
gmted3 <- raster("Data/Predictor_data/Elevation/30s030e_20101117_gmted_med075.tif") 
gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
gmted3.crop<-crop(gmted3, polygon)
gmted6.crop<-crop(gmted6, polygon)
merged.trans <- spTransform(merged, crs(gmted3.crop))
elevation.1<-raster::extract(gmted3.crop, merged.trans)
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
merge<-over(merged, bind.lakes)
merge$water.body<-ifelse(is.na(merge$TYPE)==FALSE, 1,0)
merged@data$water.body<-merge$water.body

merged@data$baseline_urban<-as.factor(merged@data$baseline_urban)
merged@data<-merged@data[,c("ea.members", "eamembers.pigs", "eamembers.cattle", "urban", "num.cattle", 
                            "num.pigs", "survey", "year", "baseline_urban", "district", "cluster", "water.body",
                            "protected", "elevation")]
writeOGR(obj=merged, dsn="Data/Exposure_data/Malawi/Created_datasets/merged_data", layer="merged_map", driver="ESRI Shapefile", overwrite_layer = TRUE)

#----------------------#
# Outcome data
#----------------------#
HAT<-read.csv("Data/Outcome_data/Malawi/HAT_Malawi_2000_2014.csv")

HAT<-HAT[which(!is.na(HAT$Latitude)),]#Remove missing location
HAT_coords<-cbind(HAT$Longitude, HAT$Latitude)

HAT$District<-as.factor(HAT$District)

levels(HAT$District)<-c("", "Kasungu", "Mzimba", "Mzimba", "Nkhotakota", "Rumphi")


HAT.sp <- SpatialPointsDataFrame(coords = HAT_coords, data = HAT,
                                 proj4string = crs(geographic))

#(1) Protected areas
merge.protected<-over(HAT.sp, protected)
protectedv2<-ifelse(is.na(merge.protected$WDPAID)==FALSE, 1,0)

#(2) Elevation
merge.trans <- spTransform(HAT.sp, crs(gmted3.crop))
gmted3.crop<-crop(gmted3, districts)
gmted6.crop<-crop(gmted6, districts)
elevation.1<-raster::extract(gmted3, merge.trans)
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


writeOGR(obj=HAT.sp, dsn="Data/Outcome_data/Malawi/HAT_Data", layer="Malawi/HAT_data", driver="ESRI Shapefile", overwrite_layer=TRUE)

#----------------------#
#Random field prediction
#----------------------#
mapExt=extent(districts)
nx=200
ny=450
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
elevation.1<-raster::extract(gmted3, merge.trans)
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

writeOGR(obj=gridP.sp, dsn="Data/Exposure_data/Malawi/Created_datasets/prediction_data", layer="prediction_data", driver="ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(obj=protected.crop, dsn="Data/Exposure_data/Malawi/Created_datasets/prediction_data_protected", layer="prediction_data_protected", driver="ESRI Shapefile", overwrite_layer=TRUE)
writeRaster(gmted3.crop, "Data/Exposure_data/Malawi/Created_datasets/prediction_data_elevation1", format="GTiff", overwrite=TRUE)
writeRaster(gmted6.crop, "Data/Exposure_data/Malawi/Created_datasets/prediction_data_elevation2", format="GTiff", overwrite=TRUE)
writeOGR(obj=lakes.crop, dsn="Data/Exposure_data/Malawi/Created_datasets/prediction_data_lakes", layer="prediction_data_lakes", driver="ESRI Shapefile", overwrite_layer=TRUE)

#---------------#
#  LandScan
#---------------#
#library(raster)
#for (year in 2000:2018){
#  dpath<-paste0("Data/denominator data/LandScan Global ",year, "/", "lspop", year)
#  x <- new("GDALReadOnlyDataset", dpath)
#  getDriver(x)
#  getDriverLongName(getDriver(x))
#  xx<-asSGDF_GROD(x)
#  r <- raster(xx)
#  landScancrop2<-crop(r, districts)
#  filename<-paste0("Data/denominator data/LandScan/Malawi/LandScan_", year)
#  writeRaster(landScancrop2, filename, format="GTiff", overwrite=TRUE)
#}