library(readstata13)
library(dplyr)
library(uwIntroStats)
library(psych)
library(survey)
library(foreign)
library(rgdal)
library(INLA)
library(RColorBrewer)
library(classInt)

source('codes/my_functions.R')

binaryYN<-function(x){
  if(nlevels(as.factor(x))>2){
    ifelse(x=="No",0,
           ifelse(x=="Yes",1,NA))
  }else{
    ifelse(x=="Yes",1,0)
  }
}

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

geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
states<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/gadm36_SSD_shp", layer="gadm36_SSD_1")
counties<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/gadm36_SSD_shp", layer="gadm36_SSD_2")
subcounties<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/gadm36_SSD_shp", layer="gadm36_SSD_3")
proj<-proj4string(states)

all<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/ssd_admbndl_admall_ssnbs_itos_20160114", layer="ssd_admbndl_admALL_ssnbs_itos_20160114")

humdata<-readOGR("Data/Exposure_data/South_Sudan/shapefiles/ssd_admbnda_adm2_imwg_nbs_20180817", layer="ssd_admbnda_adm2_imwg_nbs_20180817")

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
                      bicycle=SS2008A_BICYCLE, canoe_boat=SS2008A_BOAT, animal_cart=SS2008A_ANIMTRAN, tractor=SS2008A_TRACTOR, 
                      tv=SS2008A_TV, radio=SS2008A_RADIO, mobile=SS2008A_MOBILE, landline=SS2008A_PHONE, 
                      computer=SS2008A_COMPUTER, fridge=SS2008A_REFRIDGE, satellite=SS2008A_SATELL, fan=SS2008A_FAN, 
                      aircon=SS2008A_AIRCON, land=SS2008A_LANDOWN, landarea=SS2008A_AREACULT,
                      )

#---------------------#
# Data prep
#---------------------#
census$rooms<-as.numeric(census$rooms)
census$memrooms<-census$rooms/census$members

census$landarea<-as.numeric(census$landarea)

census$toilet<-as.factor(census$toilet)
levels(census$toilet)<-c("Pit_latrine_private", "Pit_latrine_shared", "Flush_private", 
                         "Flush_shared", "Bucket", "None")
levels(census$toilet)<-c(3,2,5,4,1,0)
census$toilet<-as.numeric(as.character(census$toilet))

census$water<-as.factor(census$water)
levels(census$water)<-c("Shared_standpipe", "Shared_standpipe", "Borehole", "Borehole", 
                        "Hand_pump", "Shared_standpipe", "Dug_well", "Open", "Open_w.filter", 
                        "Open", "Open", "Tanker_borehole", "Tanker_surface")
levels(census$water)<-c(3,3,2,2,2,3,1,0,1,0,0,2,1)
census$water<-as.numeric(as.character(census$water))

census$fuel<-as.factor(census$fuel)
levels(census$fuel)<-c("Firewood", "Charcoal", "Gas", "Electricity", "Dung", "Grass", "No_cooking")
levels(census$fuel)<-c(1,2,3,4,0,0,NA)
census$fuel<-as.numeric(as.character(census$fuel))

census$adm1<-as.factor(census$adm1)
levels(census$adm1)<-c("Upper Nile", "Jonglei", "Unity", 
                       "Warab", "Northern Bahr El Ghazal", 
                       "Western Barh El Ghazal", "Lakes", "Western Equatoria", 
                       "Central Equatoria", "Eastern Equatoria")

census$adm2<-as.factor(census$adm2)
levels(census$adm2)<-c("Renk", "Manyo, Melut", "Fashoda", 
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

census$land<-ifelse(census$land==1,1,
                    ifelse(census$land==2|census$land==3|census$land==4,0,NA))

binvars<-as.data.frame(apply(census[,c("own", "car", "motorcycle", "bicycle", "canoe_boat", "animal_cart", 
                                      "tractor", "tv", "radio", "mobile", "landline","computer", "fridge", 
                                      "satellite", "fan", "aircon", "urban")], MARGIN=2, binary))

binvars$SS2008A_DWNUM<-census$SS2008A_DWNUM

binvars$vehicle<-ifelse(binvars$car==1,3, 
                        ifelse(binvars$motorcycle_scooter==1,2,
                               ifelse(binvars$bicycle==1|binvars$animal_cart==1,1,0)))

census<-census[,c("memrooms", "land", "landarea", "toilet", "water", "fuel", "stratum", "adm1", "adm2", "SS2008A_DWNUM", "HHWT")]

census<-merge(census, binvars[,c("SS2008A_DWNUM","own", "canoe_boat", "tractor", "tv", "radio", "mobile", "landline", "computer",
                                 "fridge", "satellite", "fan", "aircon", "urban")], by="SS2008A_DWNUM")

census$landarea<-ifelse(is.na(census$land)==TRUE, NA, census$landarea)

#Make sure all variables for factor analysis are numeric or integer
is.fact<-sapply(census, is.factor)
summary(is.fact)

is.char<-sapply(census, is.character)
summary(is.char) 

is.log<-sapply(census, is.logical) 
summary(is.log) 

is.num<-sapply(census, is.numeric)
summary(is.num)

#----------------#
# Factor analysis
#----------------#
#(1) Prepare data
#(1a) Remove variables with high misssingness
remove<-c("SS2008A_DWNUM", "stratum", "HHWT", "adm1", "adm2", "urban")
idx<-which(names(census)%ni%remove)
data=census[,idx]
miss<-propmiss(data)
lowmiss<-miss[which(miss$propmiss<=0.05),] 
names.old<-names(data)
names.new<-lowmiss$variable
index<-which(names.old %in% names.new)
data<-data[,index]

#(1b) Remove variables with 0 variance
is.nonzeroSD<-sapply(data, nonzeroSD)
scratch<-as.data.frame(is.nonzeroSD)
scratch$nonzeroSD<-ifelse(scratch$is.nonzeroSD==TRUE,1,0)
scratch2<-scratch[which(scratch$nonzeroSD==1),]
names.old<-names(data)
names.new<-rownames(scratch2) 
index<-which(names.old %in% names.new)
new_dat<-data[,index]
new_dat$survey<-data$survey

#(2) Create correlation matrix
data=new_dat
matrix<-mixedCor(data= data, ncat=2)

#(3) Do factor analysis (from documentation: principal components extraction using correlation method with one factor
#extracted, substitution of mean for missing values, estimation of the factor scores using the regression
#method)
fit <- fa(matrix$rho, 1, n.obs=nrow(data), fm="pa", scores="regression", missing="TRUE", weight="NULL") 

score<-as.data.frame(factor.scores(data,fit, impute="mean")$scores)$PA1

census$comb.score<-score

#--------------------#
# Check wealth scores
#--------------------#
keep<-which(names(census) %ni% remove)
census_check<-census[,keep]

for(n in 1:ncol(census_check)){
  variable<-names(census_check)[n]
  model<-lm(census$comb.score~census_check[,n])
  res<-model$coefficients[2]
  if(!is.na(res)){
    if(res<0){
      miss<-(as.numeric((summary(is.na(census_check[,n]))[3]))/nrow(census_check))*100
      message(paste0(variable, ":(miss ",round(miss),"%)"))
      print(res)
    }
  }
} 

#---------------------#
# Get direct estimates
#---------------------#
des<-svydesign(ids=~1, weights=~HHWT, data=census)
wealth.census<-svyby(~comb.score, ~adm2, des, svymean)
wealth.adm2<-census %>% group_by(adm2)%>%
  summarise(adm1=head(adm1,1))%>% 
  arrange(adm2)
wealth.adm2<-as.data.frame(wealth.adm2)

wealth.census<-merge(wealth.census, wealth.adm2[,c("adm1", "adm2")], by="adm2")

wealth.census<-dplyr::rename(wealth.census, wealth.se=se, wealth=comb.score, county=adm2, state=adm1)
wealth.census$year<-rep(2008, nrow(wealth.census))
wealth.census$survey<-rep("Census", nrow(wealth.census))

#----------------------------------#
# Merge with census data and export
#----------------------------------#
census<-merge(census, wealth.census, by.x="new_adm2", by.y="county")
writeOGR(obj=census, dsn="Data/Confounder_data/South_Sudan/census_wealth", layer="census_wealth", driver="ESRI Shapefile", overwrite_layer = TRUE)



