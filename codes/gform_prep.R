rm(list = ls())
library(INLA)
library(rgdal)
library(sp)
library(splines)
library(raster)
library(rgeos)
library(dismo)
library(dplyr)
library(geosphere)
library(geodist)
set.seed(1154)
inla.setOption("num.threads", 4)
'%ni%'<-Negate('%in%')

country<-"DRC"

geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

if(country=="Malawi"){
  UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
  districts<-readOGR("Data/Exposure_data/Malawi/shapefiles/mwi_admbnda_adm2_nso_20181016", layer="mwi_admbnda_adm2_nso_20181016")
  proj<-proj4string(districts)
  espg<-("+init=epsg:20936")
  
  HAT<-readOGR(paste0("Data/Outcome_data/", country, "/HAT_final"), layer="HAT_final")
  colnames(HAT@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                          "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                          "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                          "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                          "wealth_975","district", "travel_hrs","protected5km", "parish", "conflict", "disaster", "dis_num", 
                          "dis_deaths", "NDVI", "LST", "elevation")
  years=c(2000:2018)
}

if(country=="Uganda"){
  UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
  districts<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_4")
  proj<-proj4string(districts)
  espg<-("+init=epsg:21095")
  
  HAT.g<-readOGR(paste0("Data/Outcome_data/", country, "/HAT_final"), layer="HAT_final_gamb")
  colnames(HAT.g@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                        "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                        "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                        "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                        "wealth_975","district", "travel_hrs","protected5km", "parish", "conflict", "disaster", "dis_num", 
                        "dis_deaths", "NDVI", "LST", "elevation")
  
  HAT.r<-readOGR(paste0("Data/Outcome_data/", country, "/HAT_final"), layer="HAT_final_rhod")
  colnames(HAT.r@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                          "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                          "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                          "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                          "wealth_975","district", "travel_hrs","protected5km", "parish", "conflict", "disaster", "dis_num", 
                          "dis_deaths", "NDVI", "LST", "elevation")
  years=c(2000:2018)
}

if(country=="DRC"){
  UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
  districts<-readOGR("Data/Exposure_data/DRC/shapefiles/gadm36_COD_shp", layer="gadm36_COD_2") 
  proj<-proj4string(districts)
  espg<-("+init=epsg:4049")

  HAT<-readOGR(paste0("Data/Outcome_data/", country, "/HAT_final"), layer="HAT_final_gamb")
  colnames(HAT@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                        "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                        "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                        "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                        "wealth_975", "travel_hours", "protected5km", "adm2", "conflict", "district", "disaster", "dis_num", 
                        "dis_deaths", "interaction", "NDVI", "LST", "elevation")
  years=c(2008:2015)
}

#------------------------------------------------------------#
# Interference set definition: Euclidean distance (total #)
#------------------------------------------------------------#
if(country=="DRC"){
  HAT@data$prox1km<-rep(NA, nrow(HAT@data))
}

if(country!="Uganda"){
  HAT@data$st.cluster<-1:nrow(HAT@data)
  HAT@data$cluster<-rep(NA, nrow(HAT@data))
  
  years<-unique(HAT@data$year)
  
  for (year in years){
    sub<-HAT[which(HAT@data$year==year),]
    idx<-sub@data$st.cluster
    cluster<-1:nrow(sub@data)
    HAT@data$cluster[which(HAT@data$st.cluster%in%idx)]<-cluster
    print(year)
  }
  
  #PO clusters
  if(country!="DRC"){
    HAT.small<-HAT[which(HAT@data$year==min(years)),]
    mdist <- distm(HAT.small)
    hc <- hclust(as.dist(mdist), method="complete")
    HAT.small@data$dist_clust <- cutree(hc, h=5000)
  }else{
    HAT.small<-HAT[which(HAT@data$year==min(years)),]
    HAT.small@data$dist_clust<-rep(NA, nrow(HAT.small@data))
    for(i in 1:length(unique(HAT.small@data$adm2))){
      idx<-which(HAT.small@data$adm2==unique(HAT.small@data$adm2)[i])
      sub<-HAT.small[idx,]
      if(nrow(sub)>1){
        mdist <- distm(sub)
        hc <- hclust(as.dist(mdist), method="complete")
        dist_clust <- cutree(hc, h=5000)
        dist_clust<-paste(unique(HAT.small@data$adm2)[i], dist_clust, sep=".")
        HAT.small@data$dist_clust[idx]<-dist_clust
        print(i)
      }
    }
  }
  
  HAT<-merge(HAT, HAT.small@data[,c("cluster", "dist_clust")], by="cluster")
  
  HAT@data$interaction<-interaction(HAT@data$dist_clust, HAT@data$year)
  
  HAT.collapse<-HAT@data%>%group_by(interaction)%>%summarise(cattle_POmean=mean(density_med_c, na.rm=T), cattle_POsum=sum(density_med_c, na.rm=T),
                                                              pigs_POmean=mean(density_med_p, na.rm=T), pigs_POsum=sum(density_med_p, na.rm=T))%>% as.data.frame()

  HAT<-merge(HAT, HAT.collapse[,c("interaction", "cattle_POmean", "cattle_POsum", 
                                  "pigs_POmean", "pigs_POsum")], by="interaction")
  
  #Proximity to HAT Atlas data
  if(country=="DRC"){
    HAT.small<-HAT[which(HAT@data$year==min(years)),]
    HATsmall.pred<-HAT.small[which(HAT.small@data$source=="pred"),]
    HATsmall.pred<-spTransform(HATsmall.pred, UTM)
    HATsmall.pred@data$prox1km<-rep(1, nrow(HATsmall.pred))
    
    HATsmall.atlas<-HAT.small[which(HAT.small@data$source=="WHO"),]
    HATsmall.atlas<-spTransform(HATsmall.atlas, UTM)
    HATsmall.atlas@data$prox1km<-rep(0, nrow(HATsmall.atlas))
    
    for(i in 1:length(unique(HAT.small@data$adm2))){
      idx.pred<-which(HATsmall.pred@data$adm2==unique(HATsmall.pred@data$adm2)[i])
      idx.atlas<-which(HATsmall.atlas@data$adm2==unique(HATsmall.atlas@data$adm2)[i])
      sub.pred<-HATsmall.pred[idx.pred,]
      sub.atlas<-HATsmall.atlas[idx.atlas,]
      
      if((nrow(sub.pred)>0 & nrow(sub.atlas)>0)){
        points_matrix<-gWithinDistance(sub.atlas, sub.pred, dist=1, byid=TRUE)
        points_matrix[lower.tri(points_matrix, diag=TRUE)]<-NA
        v<-rowSums(points_matrix, na.rm=T)==0
        sub.pred@data$prox1km[v]<-0
        HATsmall.pred@data$prox1km[idx.pred]<-sub.pred@data$prox1km
      }
    }
    
    HATsmall<-rbind(HATsmall.pred, HATsmall.atlas)

    HAT<-merge(HAT, HATsmall@data[,c("cluster", "prox1km")], by="cluster")
    
    HAT@data$prox1km<-HAT@data$prox1km.y
    
    cIdx<-c(which(names(HAT@data)=="prox1km.x"), which(names(HAT@data)=="prox1km.y"))
    HAT@data<-HAT@data[,-cIdx]
    
    HAT<-HAT[which(HAT@data$prox1km==0),]
  }
  
  writeOGR(obj=HAT, dsn=paste0("Data/Outcome_data/", country, "/HAT_final"), layer="HAT_finalPO", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
}

if(country=="Uganda"){
  
  HAT.g@data$st.cluster<-1:nrow(HAT.g@data)
  HAT.g@data$cluster<-rep(NA, nrow(HAT.g@data))
  
  years<-unique(HAT.g@data$year)
  
  for (year in years){
    sub<-HAT.g[which(HAT.g@data$year==year),]
    idx<-sub@data$st.cluster
    cluster<-1:nrow(sub@data)
    HAT.g@data$cluster[which(HAT.g@data$st.cluster%in%idx)]<-cluster
  }
  
  HAT.small<-HAT.g[which(HAT.g@data$year==min(years)),]
  mdist.g <- distm(HAT.small)
  hc.g <- hclust(as.dist(mdist.g), method="complete")
  HAT.small@data$dist_clust <- cutree(hc.g, h=5000)

  HAT.g<-merge(HAT.g, HAT.small@data[,c("cluster", "dist_clust")], by="cluster")
    
  HAT.g@data$cattle_POmean<-HAT.g@data$cattle_POsum<-HAT.g@data$pigs_POmean<-HAT.g@data$pigs_POsum<-rep(NA, nrow(HAT.g@data))
  
  for(i in 1:nrow(HAT.g@data)){
    sub<-subset(HAT.g@data, dist_clust==HAT.g@data$dist_clust[i])
    sub<-subset(sub, cluster!=HAT.g@data$cluster[i])
    HAT.g@data$cattle_POmean[i]<-mean(sub$density_med_c, na.rm=T)
    HAT.g@data$cattle_POsum[i]<-sum(sub$density_med_c, na.rm=T)
    HAT.g@data$pigs_POmean[i]<-mean(sub$density_med_p, na.rm=T)
    HAT.g@data$pigs_POsum[i]<-sum(sub$density_med_p, na.rm=T)
  }
  
  writeOGR(obj=HAT.g, dsn=paste0("Data/Outcome_data/", country, "/HAT_final"), layer="HAT_finalPO_gamb", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  HAT.r@data$st.cluster<-1:nrow(HAT.r@data)
  HAT.r@data$cluster<-rep(NA, nrow(HAT.r@data))
  
  years<-unique(HAT.r@data$year)
  
  for (year in years){
    sub<-HAT.r[which(HAT.r@data$year==year),]
    idx<-sub@data$st.cluster
    cluster<-1:nrow(sub@data)
    HAT.r@data$cluster[which(HAT.r@data$st.cluster%in%idx)]<-cluster
  }
  
  HAT.small<-HAT.r[which(HAT.r@data$year==2000),]
  mdist.r <- distm(HAT.small)
  hc.r <- hclust(as.dist(mdist.r), method="complete")
  HAT.small@data$dist_clust <- cutree(hc.r, h=5000)
  
  HAT.r<-merge(HAT.r, HAT.small@data[,c("cluster", "dist_clust")], by.x="cluster", by.y="cluster")
  
  HAT.r@data$cattle_POmean<-HAT.r@data$cattle_POsum<-HAT.r@data$pigs_POmean<-HAT.r@data$pigs_POsum<-rep(NA, nrow(HAT.r@data))
  
  for(i in 1:nrow(HAT.r@data)){
    sub<-subset(HAT.r@data, dist_clust==HAT.r@data$dist_clust[i])
    sub<-subset(sub, cluster!=HAT.r@data$cluster[i])
    HAT.r@data$cattle_POmean[i]<-mean(sub$density_med_c, na.rm=T)
    HAT.r@data$cattle_POsum[i]<-sum(sub$density_med_c, na.rm=T)
    HAT.r@data$pigs_POmean[i]<-mean(sub$density_med_p, na.rm=T)
    HAT.r@data$pigs_POsum[i]<-sum(sub$density_med_p, na.rm=T)
  }
  
  writeOGR(obj=HAT.r, dsn=paste0("Data/Outcome_data/", country, "/HAT_final"), layer="HAT_finalPO_rhod", driver="ESRI Shapefile", overwrite_layer=TRUE)
}



