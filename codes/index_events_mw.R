library(stringr)
library(rgdal)
library(cobalt)
library(rgeos)
library(sp)
UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")

conflict<-read.csv("Data/Confounder_data/ged191.csv", header=T)
conflict_malawi<-subset(conflict, country=="Malawi")

natural<-read.csv("Data/Confounder_data/Julianne_Meisner_2019-11-15.csv", header=F)
names(natural) <- lapply(natural[2, ], as.character)
natural <- natural[-1,] 
natural <- natural[-1,] 
natural_malawi<-subset(natural, Country=="Malawi")

natural_malawi$start<-str_sub(natural_malawi$`Start date`, start= -4)
natural_malawi$end<-str_sub(natural_malawi$`End date`, start= -4)
natural_malawi<-subset(natural_malawi, end>=2000)

natural_malawi$affected<-natural_malawi$`Total affected`
natural_malawi$affected<-ifelse(natural_malawi$affected=="Total affected", NA, natural_malawi$affected)
natural_malawi$affected<-as.numeric(as.character(natural_malawi$affected))
#Subset to natural disasters affecting 10 or more peoples
natural_malawi<-natural_malawi[which(natural_malawi$affected>10),]

natural_malawi$deaths<-natural_malawi$`Total deaths`
natural_malawi$type<-natural_malawi$`Disaster type`
natural_malawi$location<-droplevels(natural_malawi$Location)
locations<-natural_malawi$location
natural_malawi$location<-as.character(natural_malawi$location)
natural_malawi$location[6]<-("Balaka, Blantyre, Chikhwawa, Chiradzulu, Machinga, Mangochi, Mulanje, Mwanza, Neno, Nsanje, Phalombe, Thyolo, Zomba")
natural_malawi$location[17]<-("Dedza, Dowa, Kasungu, Lilongwe, Mchinji, Nkhotakota, Ntcheu, Ntchisi, Salima, Balaka, Blantyre, Chikhwawa, Chiradzulu, Machinga, Mangochi, Mulanje, Mwanza, Neno, Nsanje, Phalombe, Thyolo, Zomba")
natural_malawi$location[43]<-("Chitipa, Karonga, Likomam, Mzimba, Nkhotakota, Rumphi, Dedza, Dowa, Kasungu, Lilongwe, Mchinji, Nkhotakota, Ntcheu, Ntchisi, Salima, Balaka, Blantyre, Chikhwawa, Chiradzulu, Machinga, Mangochi, Mulanje, Mwanza, Neno, Nsanje, Phalombe, Thyolo, Zomba")
natural_malawi$location[42]<-("Chitipa, Karonga, Likomam, Mzimba, Nkhotakota, Rumphi, Dedza, Dowa, Kasungu, Lilongwe, Mchinji, Nkhotakota, Ntcheu, Ntchisi, Salima, Balaka, Blantyre, Chikhwawa, Chiradzulu, Machinga, Mangochi, Mulanje, Mwanza, Neno, Nsanje, Phalombe, Thyolo, Zomba")

#Read in outcome data to figure out which locations I need
HAT<-readOGR("Data/Outcome_data/Malawi/HAT_regress", layer="HAT_regress", p4s=proj)
colnames(HAT@data)<-c("District","year","focus","cases","density_med_c", "density_sd_c",  
                      "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                      "density_q025_p", "density_q975_p", "counts_med_c", "counts_sd_c",
                      "counts_q025_c", "counts_q975_c", "counts_med_p", "counts_sd_p", 
                      "counts_q025_p", "counts_q975_p", "source", "cluster", "intercept",
                      "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                      "wealth_975","prox1km", "travel_hrs","protected5km")

#--------------------------------------------#
#Subset to districts with cases
#--------------------------------------------#

pos<-HAT[which(HAT@data$cases>1),]
district<-unique(pos@data$District)
HAT<-HAT[which(HAT@data$District%in%district),]
HAT@data$District<-droplevels(HAT@data$District)

keep=c()
for(i in 1:nrow(natural_malawi)){
  loc=natural_malawi$location[i]
  d1<-grepl(district[1], loc, fixed=TRUE)
  d2<-grepl(district[2], loc, fixed=TRUE)
  d3<-grepl(district[3], loc, fixed=TRUE)
  d4<-grepl(district[4], loc, fixed=TRUE)
  if(d1==TRUE|d2==TRUE|d3==TRUE|d4==TRUE){
    keep=c(keep,i)
  }
}

natural_malawi<-natural_malawi[keep,]

natural_malawi$Kasungu<-rep(NA, nrow(natural_malawi))
natural_malawi$Mzimba<-rep(NA, nrow(natural_malawi))
natural_malawi$Nkhotakota<-rep(NA, nrow(natural_malawi))
natural_malawi$Rumphi<-rep(NA, nrow(natural_malawi))

for(i in 1:nrow(natural_malawi)){
  loc=natural_malawi$location[i]
  if(grepl(district[1], loc, fixed=TRUE)==TRUE){
    natural_malawi$Kasungu[i]=1
  }else{
    natural_malawi$Kasungu[i]=0
  }
  if(grepl(district[2], loc, fixed=TRUE)==TRUE){
    natural_malawi$Mzimba[i]=1
  }else{
    natural_malawi$Mzimba[i]=0
  }
  if(grepl(district[3], loc, fixed=TRUE)==TRUE){
    natural_malawi$Nkhotakota[i]=1 
  }else{
    natural_malawi$Nkhotakota[i]=0
  }
  if(grepl(district[4], loc, fixed=TRUE)==TRUE){
    natural_malawi$Rumphi[i]=1
  }else{
    natural_malawi$Rumphi[i]=0
  }
}

#Scan location, make sure there's nothing finer than districts included
natural_malawi<-natural_malawi[,c("start", "end", "affected", "deaths", "type", "Kasungu", "Mzimba", "Nkhotakota", "Rumphi")]

HAT@data$disaster<-rep(NA, nrow(HAT@data))
HAT@data$dis_num<-rep(NA, nrow(HAT@data))
HAT@data$dis_deaths<-rep(NA, nrow(HAT@data))

for(i in 1:nrow(HAT@data)){
  district<-HAT@data$District[i]
  year<-HAT@data$year[i]
  cols<-which(names(natural_malawi)==district)
  rows<-which(natural_malawi[,cols]==1&natural_malawi$start<=year&natural_malawi$end>=year)
  if(length(rows)!=0){
    affected<-natural_malawi$affected[rows]
    deaths<-natural_malawi$deaths[rows]
    type<-natural_malawi$type[rows]
  HAT@data$disaster[i]=as.character(type)
  HAT@data$dis_num[i]=as.numeric(as.character(affected))
  HAT@data$dis_deaths[i]=as.numeric(as.character(deaths))
  }else{
  HAT@data$disaster[i]="None"
  HAT@data$dis_num[i]="0"
  HAT@data$dis_deaths[i]="0"
  }
}

writeOGR(obj=HAT, dsn="Outcome_data/Malawi/HATData_dis", layer="HATData_dis", driver="ESRI Shapefile", overwrite_layer=TRUE)
