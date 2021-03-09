rm(list = ls())
library(rgdal)
library(INLA)
library(dplyr)
library(ggmap)
library(ggplot2)
library(gridExtra)
library(gdalUtils)
library(classInt)
library(RColorBrewer)
library(xtable)
library(stringr)

simple<-FALSE

source("codes/functions_prediction.R")

country<-"DRC"

folderData<-paste0("plots/", country, "/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

if(country=="Malawi"){
  plot.years<-all.years<-c(2000:2020)
  m=3
  height<-800
  width<-420
  nx=200
  ny = 450
  modisIdx<-c("h21v09", "h21v10")
  country.short="mw"
  districts<-readOGR("Data/Exposure_data/Malawi/shapefiles/mwi_admbnda_adm2_nso_20181016", layer="mwi_admbnda_adm2_nso_20181016")
  humdata<-readOGR("Data/Exposure_data/Malawi/shapefiles/mwi_popa_adm3_tradauth_geonode_nso2008_ocha", layer="mwi_popa_adm3_tradauth_geonode_nso2008_ocha")
  HAT<-readOGR("Data/Outcome_data/Malawi/HAT_final", layer="HAT_finalPO", p4s=proj4string(humdata))
  colnames(HAT@data)<-c("cluster", "year","focus","cases","density_med_c", "density_sd_c",  
                        "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                        "density_q025_p", "density_q975_p", "source", "prox1km", "intercept",
                        "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                        "wealth_q975", "travel_hrs","protected","district", "disaster", "dis_num", "dis_deaths", "NDVI", "LST", "elevation", 
                        "st_cluster", "dist_cluster", "cattle_PO_mean", "cattle_PO_sum", "pigs_PO_mean", "pigs_PO_sum")
  HAT.districts<-over(HAT, humdata)
  HAT_dists<-unique(HAT.districts$TRAD_AUTH)
  HAT_districts<-unique(HAT.districts$DISTRICT)
  dists<-humdata[which(humdata@data$TRAD_AUTH%in%HAT_dists),]
  
  HAT_data<-readOGR(paste0("Data/Outcome_data/HAT_Data/", country), layer="HAT_data")
  
  merged_data<-readOGR(paste("Data/Exposure_data/Malawi/Created_datasets/merged_data", sep=""), layer="merged_map") 
  colnames(merged_data@data)<-c("members", "mem.pigs", "mem.cattle",
                                "urban", "n.cattle", "n.pigs", "survey", "year", "base_urb",
                                "district", "cluster", "water.body", "protected", "elevation")
}

if(country=="DRC"){
  plot.years<-all.years<-c(2008:2015)
  modisIdx<-c("h20v08", "h20v09", "h21v08", "h21v09", "h21v10", "h19v09", "h19v08", "h19v10", "h20v10")
  m=2
  height<-800
  width<-800
  nx= 1170
  ny= 1090
  simple<-TRUE
  country.short="drc"
  districts<-readOGR("Data/Exposure_data/DRC/shapefiles/gadm36_COD_shp", layer="gadm36_COD_2") 
  humdata<-districts
  
  merged_data<-readOGR(paste("Data/Exposure_data/DRC/Created_datasets/merged_data", sep=""), layer="merged_map") 
  colnames(merged_data@data)<-c("members", "mem.pigs", "mem.cattle",
                                "urban", "n.cattle", "n.pigs", "survey", "year", "base_urb",
                                "district", "cluster", "water.body", "protected", "elevation")
  
  HAT<-readOGR("Data/Outcome_data/DRC/HAT_final", layer="HAT_finalPO", p4s=proj4string(humdata))
  colnames(HAT@data)<-c("cluster", "interaction", "year","focus","cases","density_med_c", "density_sd_c",  
                        "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                        "density_q025_p", "density_q975_p", "source", "prox_1km", "intercept",
                        "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                        "wealth_q975", "travel_hrs","protected","adm2", "conflict", "district", "disaster", "dis_num", "dis_deaths", "NDVI", "LST", "elevation", 
                        "st_cluster", "dist_cluster", "cattle_PO_mean", "cattle_PO_sum", "pigs_PO_mean", "pigs_PO_sum", "prox1km")
  
  HAT.districts<-over(HAT, humdata)

  HAT_dists<-unique(HAT.districts$NAME_2)

  dists<-humdata[which(humdata@data$NAME_2%in%HAT_dists),]
  
  HAT_data<-readOGR(paste0("Data/Outcome_data/HAT_Data/", country), layer="HAT_data")
}

if(country=="Uganda"){
  plot.years<-all.years<-c(2000:2020)
  m=6
  country.short="ug"
  height<-420
  width<-420
  nx=330
  ny=330
  modisIdx<-c("h20v08", "h20v09", "h21v08", "h21v09")
  districts<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_4")
  humdata<-readOGR("Data/Exposure_data/Uganda/shapefiles/gadm36_UGA_shp", layer="gadm36_UGA_2")
  HAT.g<-readOGR("Data/Outcome_data/Uganda/HAT_final", layer="HAT_final_gamb", p4s=proj4string(humdata))
  colnames(HAT.g@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                          "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                          "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                          "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                          "wealth_q975","district", "travel_hrs","protected","parish", "conflict","disaster", "dis_num", "dis_deaths", "NDVI", "LST", "elevation")
  HAT.r<-readOGR("Data/Outcome_data/Uganda/HAT_final", layer="HAT_final_rhod", p4s=proj4string(humdata))
  colnames(HAT.r@data)<-c("year","focus","cases","density_med_c", "density_sd_c",  
                          "density_q025_c","density_q975_c","density_med_p","density_sd_p",  
                          "density_q025_p", "density_q975_p", "source", "prox1km", "cluster", "intercept",
                          "WorldPop","LandScan","urban","nightlights","wealth_med","wealth_sd","wealth_q025",   
                          "wealth_q975","district", "travel_hrs","protected","parish", "conflict","disaster", "dis_num", "dis_deaths", "NDVI", "LST", "elevation")
  HAT.districts.r<-over(HAT.r, humdata)
  HAT.districts.g<-over(HAT.g, humdata)
  
  HAT_dists.r<-unique(HAT.districts.r$NAME_2)
  HAT_dists.g<-unique(HAT.districts.g$NAME_2)
  
  dists.r<-humdata[which(humdata@data$NAME_2%in%HAT_dists.r),]
  dists.g<-humdata[which(humdata@data$NAME_2%in%HAT_dists.g),]
  
  HAT_data<-readOGR(paste0("Data/Outcome_data/HAT_Data/", country), layer="HAT_data")
  
  merged_data<-readOGR(paste("Data/Exposure_data/Uganda/Created_datasets/merged_data", sep=""), layer="merged_map") 
  colnames(merged_data@data)<-c("members", "mem.pigs", "mem.cattle",
                                "urban", "n.cattle", "n.pigs", "survey", "year", "base_urb",
                                "district", "cluster", "water.body", "protected", "elevation")
}

if(country=="South_Sudan"){
  country.short="SS"
  width<-420
  height<-300
  HAT_data<-readOGR("Data/Outcome_data/South_Sudan/SS_HAT", layer="SS_HAT")
  census<-readOGR("Data/Exposure_data/South_Sudan/Created_datasets/census_data", layer="census_data")
  humdata<-census
  active<-readOGR("Data/Outcome_data/South_Sudan/HAT_final_SS", layer="HAT_final_SS_active")
  HAT<-readOGR("Data/Outcome_data/South_Sudan/HAT_final_SS", layer="HAT_final_SS")
  dists<-HAT
  
  maps<-readOGR("Data/Exposure_data/South_Sudan/Created_datasets/census_data", layer="SS_HAT_distribute")
  colnames(maps@data)<-c("new_adm2", "shape_area", "ADM1_EN", "ADM1_PC", "cattle_median", "cattle_sd", "pigs_median", "pigs_sd")
}

if(country!="South_Sudan"){
  wealth_data<-readOGR(paste0("Data/Confounder_data/", country, "/wealth"), layer="wealth")
  cattle_data<-readOGR(paste("Data/Exposure_data/Created_datasets/cattle_data", sep=""), layer="cattle_data")
  pig_data<-readOGR(paste("Data/Exposure_dataCreated_datasets/pig_data", sep=""), layer="pig_data")
  proj<-proj4string(cattle_data)
  
  load(paste0("results/livestock_mapping/", country, "/prediction_cattle_", country.short, ".Rdata"))
  load(paste0("results/livestock_mapping/", country, "/prediction_pigs_", country.short, ".Rdata"))
  load(paste0("results/livesotck_wealth/", country, "/models/wealth4.Rdata"))
  wealth.model<-res
} 

#----------------#
#Plot observed data (SPDE countries only)
#----------------#
if(country!="South_Sudan"){
  png(paste0("plots/", country, "/observed_", country.short, ".png"), height=height, width=width)
  if(country=="Uganda"){
    plot(humdata)
  }else{
    plot(districts)
  }
  if(country!="DRC"){
    plot(merged_data, add=T, col="red", pch=".", cex=2)
  }else{
    plot(merged_data, add=T, col="red", pch=".", cex=5)
  }
  dev.off()
}

#----------------#
#Plot study area
#----------------#
if(country!="Uganda"){
  if(country=="South_Sudan"){
    png(paste0("plots/", country, "/Study_", country.short, ".png"), height=8, width=10, res=300, units="in")
  }else{
    png(paste0("plots/", country, "/Study_", country.short, ".png"), height=height, width=width)
  }
  plot(humdata)
  plot(dists, add=TRUE, col="red")
  dev.off()
}
if(country=="Uganda"){
  png(paste0("plots/", country, "/Study_", country.short, "_r.png"), height=height, width=width)
  plot(humdata)
  plot(dists.r, add=TRUE, col="red")
  dev.off()
  
  png(paste0("plots/", country, "/Study_", country.short, "_g.png"), height=height, width=width)
  plot(humdata)
  plot(dists.g, add=TRUE, col="red")
  dev.off()
}

#--------------------------------------------------#
#Plot active surveillance numerator and denominator
#--------------------------------------------------#
if(country=="South_Sudan"){
  active<-active[which(active@data$n2008>0),]
  vars<-paste0("cs", 2000:2014)
  nvars<-paste0("n", 2000:2014)
  active@data$cases<-(rowSums(active@data[,c(vars)]))
  active@data$n<-(rowSums(active@data[,c(nvars)]))
  
  med.palette <- brewer.pal(n = 7, name = "Reds")
  med.int <- classIntervals(active@data$cases,
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  
  png(paste0("plots/", country, "/counties_active", country.short, ".png"), height=8, width=15, res=300, units="in")
  plot(active, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
  
  #Active surveillance denominator
  med.palette <- brewer.pal(n = 7, name = "Blues")
  med.int <- classIntervals(active@data$n,
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  
  png(paste0("plots/", country, "/counties_active_sampled", country.short, ".png"), height=8, width=15, res=300, units="in")
  plot(active, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
  
  #WorldPop denominator
  med.palette <- brewer.pal(n = 7, name = "Blues")
  med.int <- classIntervals(HAT@data$WP2008,
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/counties_WP_denom", country.short, ".png"), height=8, width=15, res=300, units="in")
  plot(HAT, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
  
  #LandScan denominator
  med.palette <- brewer.pal(n = 7, name = "Blues")
  med.int <- classIntervals(HAT@data$LS2008,
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/counties_LS_denom", country.short, ".png"), height=8, width=15, res=300, units="in")
  plot(HAT, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
}

#---------------------------------#
#Plot grid and mesh
#---------------------------------#
xmin<-districts@bbox[1,1]
xmax<-districts@bbox[1,2]
ymin<-districts@bbox[2,1]
ymax<-districts@bbox[2,2]
xrange<-xmax-xmin
xmed<-xmin+(xrange/2)
yrange<-ymax-ymin
ymed<-ymin+(yrange/2)

png(paste0("plots/", country, "/grid_", country.short, ".png"), height=height, width=width)
plot(1000, xlim = c(xmed-xrange/5, xmed), ylim = c(ymed, ymed+yrange/5), cex.axis = 1.5, cex.lab = 1.5, xlab = "Longitude", ylab = "Latitude", asp = 1)
grid(nx = nx/5, ny = ny/5, lty = 1, col = "black", lwd = 0.2)
plot(humdata, add = TRUE, lwd = 1.5)
dev.off()

formedMap<-districts

if(country=="Uganda"){
  formedMap<-spTransform(districts, CRS("+proj=utm +zone=35N +ellps=WGS84"))
  formedMap2<-rgeos::gBuffer(formedMap, byid=TRUE, width=0)
  formedMap<-spTransform(formedMap2, CRS("+proj=longlat +datum=WGS84"))
}

spUnion = maptools::unionSpatialPolygons(formedMap, rep(1, dim(formedMap)[1]))
loc.mesh = matrix(0, nrow = 0, ncol = 2)
for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
  loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
}

mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                    max.edge = c(0.3,0.6), 
                    offset = c(1, 1))

png(paste0("plots/", country, "/mesh_", country.short, ".png"), height=height, width=width)
if(ymin<0){
if(country=="DRC"){
plot(1000, xlim = c(xmin, xmax), ylim = c(ymin*1.15,ymax*1.3), cex.axis = 1, cex.lab = 1, xlab = "Longitude", ylab = "Latitude", asp = 1)
}else{
plot(1000, xlim = c(xmin*0.92, xmax*1.06), ylim = c(ymin*1.1,ymax*1), cex.axis = 1, cex.lab = 1, xlab = "Longitude", ylab = "Latitude", asp = 1)
}
}else{
plot(1000, xlim = c(xmin*0.90, xmax*1.1), ylim = c(ymin*0.90,ymax*1.1), cex.axis = 1, cex.lab = 1, xlab = "Longitude", ylab = "Latitude", asp = 1)
}
plot(mesh, add=TRUE)
if(country!="DRC"){
  plot(humdata, add=TRUE, lwd=1.5)
}else{
plot(humdata, add=TRUE, lwd=2.5)
}
dev.off()

#-----------------------------#
# Descriptive statistics table
#-----------------------------#
if(country=="Malawi"){
  vars<-c("elevation", "cases", "WorldPop", "LandScan", "LST", "NDVI", "wealth", "cattle", "pigs", "drought", "flood", "storm")
  descrip_table<-as.data.frame(matrix(nrow=length(vars), ncol=2))
  colnames(descrip_table)<-c("Mean", "(sd)")
  rownames(descrip_table)<-vars
  descrip_table[1,1]<-mean(HAT@data$elevation, na.rm=T)
  descrip_table[1,2]<-sd(HAT@data$elevation, na.rm=T)
  descrip_table[2,1]<-mean(HAT@data$cases, na.rm=T)
  descrip_table[2,2]<-sd(HAT@data$cases, na.rm=T)
  descrip_table[3,1]<-mean(HAT@data$WorldPop, na.rm=T)
  descrip_table[3,2]<-sd(HAT@data$WorldPop, na.rm=T)
  descrip_table[4,1]<-mean(HAT@data$LandScan, na.rm=T)
  descrip_table[4,2]<-sd(HAT@data$LandScan, na.rm=T)
  descrip_table[5,1]<-mean(HAT@data$LST, na.rm=T)
  descrip_table[5,2]<-sd(HAT@data$LST, na.rm=T)
  descrip_table[6,1]<-mean(HAT@data$NDVI, na.rm=T)
  descrip_table[6,2]<-sd(HAT@data$NDVI, na.rm=T)
  descrip_table[7,1]<-mean(HAT@data$wealth_med, na.rm=T)
  descrip_table[7,2]<-mean(HAT@data$wealth_q975-HAT@data$wealth_q025, na.rm=T)
  descrip_table[8,1]<-mean(HAT@data$density_med_c, na.rm=T)
  descrip_table[8,2]<-mean(HAT@data$density_q975_c-HAT@data$density_q025_c, na.rm=T)
  descrip_table[9,1]<-mean(HAT@data$density_med_p, na.rm=T)
  descrip_table[9,2]<-mean(HAT@data$density_q975_p-HAT@data$density_q025_p, na.rm=T)
  descrip_table[10,1]<-table(HAT@data$disaster)[1]#Drought
  descrip_table[10,2]<-round(((table(HAT@data$disaster)[1])/nrow(HAT@data)*100))
  descrip_table[11,1]<-table(HAT@data$disaster)[2]#Flood
  descrip_table[11,2]<-round(((table(HAT@data$disaster)[2])/nrow(HAT@data)*100))
  descrip_table[12,1]<-table(HAT@data$disaster)[3]#Storm
  descrip_table[12,2]<-round(((table(HAT@data$disaster)[3])/nrow(HAT@data)*100))
  xtab<-xtable(descrip_table)
  save(file="plots/Malawi/descrip_table_rhod.Rdata",xtab)
}

if(country=="Uganda"){
  vars<-c("elevation", "cases", "WorldPop", "LandScan", "LST", "NDVI", "wealth", "cattle", "pigs", "conflict", "epidemic","flood")
  descrip_table<-as.data.frame(matrix(nrow=length(vars), ncol=2))
  colnames(descrip_table)<-c("Mean", "(sd)")
  rownames(descrip_table)<-vars
  descrip_table[1,1]<-mean(HAT.g@data$elevation, na.rm=T)
  descrip_table[1,2]<-sd(HAT.g@data$elevation, na.rm=T)
  descrip_table[2,1]<-mean(HAT.g@data$cases, na.rm=T)
  descrip_table[2,2]<-sd(HAT.g@data$cases, na.rm=T)
  descrip_table[3,1]<-mean(HAT.g@data$WorldPop, na.rm=T)
  descrip_table[3,2]<-sd(HAT.g@data$WorldPop, na.rm=T)
  descrip_table[4,1]<-mean(HAT.g@data$LandScan, na.rm=T)
  descrip_table[4,2]<-sd(HAT.g@data$LandScan, na.rm=T)
  descrip_table[5,1]<-mean(HAT.g@data$LST, na.rm=T)
  descrip_table[5,2]<-sd(HAT.g@data$LST, na.rm=T)
  descrip_table[6,1]<-mean(HAT.g@data$NDVI, na.rm=T)
  descrip_table[6,2]<-sd(HAT.g@data$NDVI, na.rm=T)
  descrip_table[7,1]<-mean(HAT.g@data$wealth_med, na.rm=T)
  descrip_table[7,2]<-sd(HAT.g@data$wealth_sd, na.rm=T)
  descrip_table[8,1]<-mean(HAT.g@data$density_med_c, na.rm=T)
  descrip_table[8,2]<-mean(HAT.g@data$density_sd_c, na.rm=T)
  descrip_table[9,1]<-mean(HAT.g@data$density_med_p, na.rm=T)
  descrip_table[9,2]<-mean(HAT.g@data$density_sd_p, na.rm=T)
  descrip_table[10,1]<-sum(HAT.g@data$conflict, na.rm=T)
  descrip_table[10,2]<-round((sum(HAT.g@data$conflict, na.rm=T)/nrow(HAT.g@data)*100),2)
  descrip_table[11,1]<-as.numeric(table(HAT.g@data$disaster)[1])#Epidemic
  descrip_table[11,2]<-round((as.numeric(table(HAT.g@data$disaster)[1])/nrow(HAT.g@data)*100),2)
  descrip_table[12,1]<-as.numeric(table(HAT.g@data$disaster)[2])#Flood
  descrip_table[12,2]<-round((as.numeric(table(HAT.g@data$disaster)[2])/nrow(HAT.g@data)*100),2)
  descrip_table.g<-descrip_table
  xtab.g<-xtable(descrip_table.g)
  save(file="plots/Uganda/descrip_table_gamb.Rdata",xtab.g)
  
  vars<-c("elevation", "cases", "WorldPop", "LandScan", "LST", "NDVI", "wealth", "cattle", "pigs", "conflict", "drought", "epidemic","flood", "landslide")
  descrip_table<-as.data.frame(matrix(nrow=length(vars), ncol=2))
  colnames(descrip_table)<-c("Mean", "(sd)")
  rownames(descrip_table)<-vars
  descrip_table[1,1]<-mean(HAT.r@data$elevation, na.rm=T)
  descrip_table[1,2]<-sd(HAT.r@data$elevation, na.rm=T)
  descrip_table[2,1]<-mean(HAT.r@data$cases, na.rm=T)
  descrip_table[2,2]<-sd(HAT.r@data$cases, na.rm=T)
  descrip_table[3,1]<-mean(HAT.r@data$WorldPop, na.rm=T)
  descrip_table[3,2]<-sd(HAT.r@data$WorldPop, na.rm=T)
  descrip_table[4,1]<-mean(HAT.r@data$LandScan, na.rm=T)
  descrip_table[4,2]<-sd(HAT.r@data$LandScan, na.rm=T)
  descrip_table[5,1]<-mean(HAT.r@data$LST, na.rm=T)
  descrip_table[5,2]<-sd(HAT.r@data$LST, na.rm=T)
  descrip_table[6,1]<-mean(HAT.r@data$NDVI, na.rm=T)
  descrip_table[6,2]<-sd(HAT.r@data$NDVI, na.rm=T)
  descrip_table[7,1]<-mean(HAT.r@data$wealth_med, na.rm=T)
  descrip_table[7,2]<-sd(HAT.r@data$wealth_sd, na.rm=T)
  descrip_table[8,1]<-mean(HAT.r@data$density_med_c, na.rm=T)
  descrip_table[8,2]<-mean(HAT.r@data$density_sd_c, na.rm=T)
  descrip_table[9,1]<-mean(HAT.r@data$density_med_p, na.rm=T)
  descrip_table[9,2]<-mean(HAT.r@data$density_sd_p, na.rm=T)
  descrip_table[10,1]<-sum(HAT.r@data$conflict, na.rm=T)
  descrip_table[10,2]<-round((sum(HAT.r@data$conflict, na.rm=T)/nrow(HAT.r@data)*100),2)
  descrip_table[11,1]<-as.numeric(table(HAT.r@data$disaster)[1])#Drought
  descrip_table[11,2]<-round((as.numeric(table(HAT.r@data$disaster)[1])/nrow(HAT.r@data)*100),2)
  descrip_table[12,1]<-as.numeric(table(HAT.r@data$disaster)[2])#Epidemic
  descrip_table[12,2]<-round((as.numeric(table(HAT.r@data$disaster)[2])/nrow(HAT.r@data)*100),2)
  descrip_table[13,1]<-as.numeric(table(HAT.r@data$disaster)[3])#Flood
  descrip_table[13,2]<-round((as.numeric(table(HAT.r@data$disaster)[3])/nrow(HAT.r@data)*100),2)
  descrip_table[14,1]<-as.numeric(table(HAT.r@data$disaster)[4])#Landslide
  descrip_table[14,2]<-round((as.numeric(table(HAT.r@data$disaster)[4])/nrow(HAT.r@data)*100),2)
  descrip_table.r<-descrip_table
  xtab.r<-xtable(descrip_table.r)
  save(file="plots/Uganda/descrip_table_rhod.Rdata",xtab.r)
}

if(country=="South_Sudan"){
  vars<-c("elevation", "cases", "n_screened", "WorldPop", "LandScan", "LST", "NDVI", "wealth", "cattle", "pigs", "war")
  descrip_table<-as.data.frame(matrix(nrow=11, ncol=2))
  colnames(descrip_table)<-c("Mean", "(sd)")
  rownames(descrip_table)<-vars
  descrip_table[1,1]<-mean(HAT@data$elevation, na.rm=T)
  descrip_table[1,2]<-sd(HAT@data$elevation, na.rm=T)
  descrip_table[3,1]<-mean(active@data$n2008, na.rm=T)
  descrip_table[3,2]<-sd(active@data$n2008, na.rm=T)
  descrip_table[2,1]<-mean(HAT@data$cs2008, na.rm=T)
  descrip_table[2,2]<-sd(HAT@data$cs2008, na.rm=T)
  descrip_table[4,1]<-mean(HAT@data$WP2008, na.rm=T)
  descrip_table[4,2]<-sd(HAT@data$WP2008, na.rm=T)
  descrip_table[5,1]<-mean(HAT@data$LS2008, na.rm=T)
  descrip_table[5,2]<-sd(HAT@data$LS2008, na.rm=T)
  descrip_table[6,1]<-mean(HAT@data$LST_2007, na.rm=T)
  descrip_table[6,2]<-sd(HAT@data$LST_2007, na.rm=T)
  descrip_table[7,1]<-mean(HAT@data$NDVI2007, na.rm=T)
  descrip_table[7,2]<-sd(HAT@data$NDVI2007, na.rm=T)
  descrip_table[8,1]<-mean(HAT@data$wealth, na.rm=T)
  descrip_table[8,2]<-mean(HAT@data$wealth_se, na.rm=T)
  descrip_table[9,1]<-mean(HAT@data$pred_c, na.rm=T)
  descrip_table[9,2]<-mean(HAT@data$pred_c_sd, na.rm=T)
  descrip_table[10,1]<-mean(HAT@data$pred_p, na.rm=T)
  descrip_table[10,2]<-mean(HAT@data$pred_p_sd, na.rm=T)
  descrip_table[11,1]<-table(HAT@data$war2008)[2]
  descrip_table[11,2]<-round(((table(HAT@data$war2008)[2])/nrow(HAT@data)*100))
  
  descriptab<-xtable(descrip_table)
  save(file="plots/South_Sudan/descrip_table.Rdata",descriptab)
}

if(country=="DRC"){
  vars<-c("elevation", "cases", "WorldPop", "LandScan", "LST", "NDVI", "wealth", "cattle", "pigs", "conflict", "earthquake", "epidemic","flood", "wildfire")
  descrip_table<-as.data.frame(matrix(nrow=length(vars), ncol=2))
  colnames(descrip_table)<-c("Mean", "(sd)")
  rownames(descrip_table)<-vars
  descrip_table[1,1]<-mean(HAT@data$elevation, na.rm=T)
  descrip_table[1,2]<-sd(HAT@data$elevation, na.rm=T)
  descrip_table[2,1]<-mean(HAT@data$cases, na.rm=T)
  descrip_table[2,2]<-sd(HAT@data$cases, na.rm=T)
  descrip_table[3,1]<-mean(HAT@data$WorldPop, na.rm=T)
  descrip_table[3,2]<-sd(HAT@data$WorldPop, na.rm=T)
  descrip_table[4,1]<-mean(HAT@data$LandScan, na.rm=T)
  descrip_table[4,2]<-sd(HAT@data$LandScan, na.rm=T)
  descrip_table[5,1]<-mean(HAT@data$LST, na.rm=T)
  descrip_table[5,2]<-sd(HAT@data$LST, na.rm=T)
  descrip_table[6,1]<-mean(HAT@data$NDVI, na.rm=T)
  descrip_table[6,2]<-sd(HAT@data$NDVI, na.rm=T)
  descrip_table[7,1]<-mean(HAT@data$wealth_med, na.rm=T)
  descrip_table[7,2]<-sd(HAT@data$wealth_med, na.rm=T)
  descrip_table[8,1]<-mean(HAT@data$density_med_c, na.rm=T)
  descrip_table[8,2]<-sd(HAT@data$density_med_c, na.rm=T)
  descrip_table[9,1]<-mean(HAT@data$density_med_p, na.rm=T)
  descrip_table[9,2]<-sd(HAT@data$density_med_p, na.rm=T)
  descrip_table[10,1]<-sum(HAT@data$conflict, na.rm=T)
  descrip_table[10,2]<-round((sum(HAT@data$conflict, na.rm=T)/nrow(HAT@data)*100),2)
  descrip_table[11,1]<-as.numeric(table(HAT@data$disaster)[1])#Earthquake
  descrip_table[11,2]<-round((as.numeric(table(HAT@data$disaster)[1])/nrow(HAT@data)*100),2)
  descrip_table[12,1]<-as.numeric(table(HAT@data$disaster)[2])#Epidemic
  descrip_table[12,2]<-round((as.numeric(table(HAT@data$disaster)[2])/nrow(HAT@data)*100),2)
  descrip_table[13,1]<-as.numeric(table(HAT@data$disaster)[3])#Flood
  descrip_table[13,2]<-round((as.numeric(table(HAT@data$disaster)[3])/nrow(HAT@data)*100),2)
  descrip_table[14,1]<-as.numeric(table(HAT@data$disaster)[5])#Wildfire
  descrip_table[14,2]<-round((as.numeric(table(HAT@data$disaster)[5])/nrow(HAT@data)*100),2)
  xtab<-xtable(descrip_table)
  save(file="plots/DRC/descrip_table.Rdata",xtab)
}

#---------#
#png map
#---------#
if(country=="South_Sudan"){
  #Cattle median
  med.palette <- brewer.pal(n = 7, name = "Purples")
  med.int <- classIntervals(round(maps@data$cattle_median,2),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/Cattle_density_median_", country.short, ".png"), width=12, height=8, units="in", res=300)
  plot(maps, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
  
  #Cattle SE
  med.palette <- brewer.pal(n = 7, name = "Purples")
  med.int <- classIntervals(round(maps@data$cattle_sd,2),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/Cattle_density_sd_", country.short, ".png"), width=12, height=8, units="in", res=300)
  plot(maps, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
  
  #Pig median
  med.palette <- brewer.pal(n = 7, name = "Purples")
  med.int <- classIntervals(round(maps@data$pigs_median,4),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/Pig_density_median_", country.short, ".png"), width=12, height=8, units="in", res=300)
  plot(maps, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
  
  #Pig SE
  med.palette <- brewer.pal(n = 7, name = "Purples")
  med.int <- classIntervals(round(maps@data$pigs_sd,2),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/Pig_density_sd_", country.short, ".png"), width=12, height=8, units="in", res=300)
  plot(maps, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
}else{
  if(country=="Uganda"){
    lakes<-readOGR("Data/Exposure_data/Uganda/Created_datasets/prediction_data_lakes", layer="prediction_data_lakes")
    lakes1<-readOGR("Data/Lakes_wetlands/GLWD-level1", layer="glwd_1", p4s=proj)
    lakes2<-readOGR("Data/Lakes_wetlands/GLWD-level2", layer="glwd_2", p4s=proj)
    lakes1.crop<-crop(lakes1, districts)
    lakes2.crop<-crop(lakes2, districts)
    bind.lakes<-bind(lakes1.crop, lakes2.crop)
    victoria<-bind.lakes[which(grepl("Victoria", bind.lakes@data$LAKE_NAME)==TRUE),]
  }
  
  spde_RE<-inla.spde2.pcmatern(mesh, prior.range=c(0.3, 0.05), prior.sigma=c(1,0.05))
  
  gridP<-readOGR(paste0("Data/Exposure_data/", country, "/Created_datasets/prediction_data", sep=""), layer = "prediction_data")
  gridP@data$water.body<-gridP@data$water_body
  
  mapExt = extent(formedMap)
  nx = nx
  ny = ny
  nPred = nx*ny
  x = seq(mapExt@xmin-0.1, mapExt@xmax+0.1, length.out = nx)
  xs = rep(x, each = ny)
  y = seq(mapExt@ymin-0.1, mapExt@ymax+0.1, length.out = ny) 
  ys = rep(y, nx)
  loc.pred = cbind(xs, ys)
  
  xMat = matrix(xs, ncol = nx)
  yMat = matrix(ys, ncol = nx)
  
  nSamples=1000
  inla.seed = as.integer(runif(1)*.Machine$integer.max)
  prevSample = NULL
  
  urbanNum = getUrbanicity(year=2010, loc=loc.pred, regMap=districts, country=country)
  urbanFac = factor(x = urbanNum,
                    levels = c(0, 1))
  data.pred<-gridP
  data.pred$urban = urbanFac
  data.pred$year = rep(2010, nrow(data.pred))
  data.pred@data$intercept<-rep(1, nrow(data.pred))
  
  load(paste0("results/livestock_mapping/", country, "/cattle_results/cattle", m, ".Rdata"))
  
  predictors<-paste(res$names.fixed, collapse=" + ")
  form.cov <- as.formula(paste("~ 0 +", predictors))
  
  for(year in plot.years){
    kIdx<-which(all.years==year)
    liveSamples = sampling(res.zp = res,
                           samp.dat= cattle_data,
                           data.pred = data.pred,
                           form.cov =form.cov,
                           nSamp = 1000,
                           spde= spde_RE,
                           simple=simple,
                           kIdx=kIdx,
                           year = year,
                           loc.pred=loc.pred,
                           prevSample=NULL)
    
    LiveRes = list(live.median = matrix(0, nrow = ny, ncol = nx),
                   live.l     = matrix(0, nrow = ny, ncol = nx),
                   live.u     = matrix(0, nrow = ny, ncol = nx))
    
    live.median = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = median, na.rm=T)
    live.median = matrix(live.median, ncol = dim(xMat)[2])
    
    live.upper = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs=0.975, na.rm=T)
    live.upper = matrix(live.upper, ncol = dim(xMat)[2])
    
    live.lower = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs=0.025, na.rm=T)
    live.lower = matrix(live.lower, ncol = dim(xMat)[2])
    
    gridM = data.frame(Longitude = gridP@coords[,1],
                       Latitude = gridP@coords[,2])
    coordinates(gridM) = ~ Longitude + Latitude
    proj4string(gridM) = proj4string(formedMap)
    mask = as.numeric(over(gridM, formedMap)[,1])
    mask[!is.na(mask)] = 1
    mask = matrix(mask, ncol = nx)
    
    live.median<-live.median*mask
    ras.median<-raster(live.median, crs=proj, xmn=min(xMat), xmx=max(xMat),
                       ymn=min(yMat), ymx=max(yMat))
    ras.median<-flip(ras.median, direction="y")
    ras.median<-crop(ras.median, districts)
    
    live.ci<-live.upper-live.lower
    live.ci<-live.ci*mask
    ras.ci<-raster(live.ci, crs=proj, xmn=min(xMat), xmx=max(xMat),
                   ymn=min(yMat), ymx=max(yMat))
    ras.ci<-flip(ras.ci, direction="y")
    ras.ci<-crop(ras.ci, districts)
    
    if(country=="Uganda"){
      ras.median<-mask(ras.median, victoria, inverse=TRUE)
      ras.ci<-mask(ras.ci, victoria, inverse=TRUE)
    }
    
    png(paste0("plots/", country, "/Cattle_density_mean_", year, "_", country.short, ".png"), height=height, width=width)
    plot(ras.median, main="")
    plot(humdata, add=T, lwd=0.5)
    dev.off()
    
    png(paste0("plots/", country, "/Cattle_density_95ci_", year, "_", country.short, ".png"), height=height, width=width)
    plot(ras.ci, main="")
    plot(humdata, add=T, lwd=0.5)
    dev.off()
    
    writeRaster(ras.median, paste0("livestock_rasters/", country, "/cattle_median_", year), format="GTiff", overwrite=TRUE)
    writeRaster(ras.ci, paste0("livestock_rasters/", country, "/cattle_95ci_", year), format="GTiff", overwrite=TRUE)
    print(year)
  }
  
  load(paste0("results/livestock_mapping/", country, "/pig_results/pigs", m, ".Rdata"))
  
  predictors<-paste(res$names.fixed, collapse=" + ")
  form.cov <- as.formula(paste("~ 0 +", predictors))
  
  for(year in plot.years){
    kIdx<-which(all.years==year)
    liveSamples = sampling(res.zp = res,
                           samp.dat= pig_data,
                           data.pred = data.pred,
                           form.cov =form.cov,
                           nSamp = 1000,
                           spde= spde_RE,
                           simple=simple,
                           kIdx=kIdx,
                           year = year,
                           loc.pred=loc.pred,
                           prevSample=NULL)
    
    LiveRes = list(live.median = matrix(0, nrow = ny, ncol = nx),
                   live.l     = matrix(0, nrow = ny, ncol = nx),
                   live.u     = matrix(0, nrow = ny, ncol = nx))
    
    live.median = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = median)
    live.median = matrix(live.median, ncol = dim(xMat)[2])
    
    live.upper = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs=0.975, na.rm=T)
    live.upper = matrix(live.upper, ncol = dim(xMat)[2])
    
    live.lower = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs=0.025, na.rm=T)
    live.lower = matrix(live.lower, ncol = dim(xMat)[2])
    
    gridM = data.frame(Longitude = gridP@coords[,1],
                       Latitude = gridP@coords[,2])
    coordinates(gridM) = ~ Longitude + Latitude
    proj4string(gridM) = proj4string(formedMap)
    mask = as.numeric(over(gridM, formedMap)[,1])
    mask[!is.na(mask)] = 1
    mask = matrix(mask, ncol = nx)
    
    live.median<-live.median*mask
    ras.median<-raster(live.median, crs=proj, xmn=min(xMat), xmx=max(xMat),
                       ymn=min(yMat), ymx=max(yMat))
    ras.median<-flip(ras.median, direction="y")
    ras.median<-crop(ras.median, districts)
    
    live.ci<-live.upper-live.lower
    live.ci<-live.ci*mask
    ras.ci<-raster(live.ci, crs=proj, xmn=min(xMat), xmx=max(xMat),
                   ymn=min(yMat), ymx=max(yMat))
    ras.ci<-flip(ras.ci, direction="y")
    ras.ci<-crop(ras.ci, districts)
    
    if(country=="Uganda"){
      ras.median<-mask(ras.median, victoria, inverse=TRUE)
      ras.ci<-mask(ras.ci, victoria, inverse=TRUE)
    }
    
    png(paste0("plots/", country, "/Pig_density_mean_", year, "_", country.short, ".png"), height=height, width=width)
    plot(ras.median, main="")
    plot(humdata, add=T, lwd=0.5)
    dev.off()
    
    png(paste0("plots/", country, "/Pig_density_95ci_", year, "_", country.short, ".png"), height=height, width=width)
    plot(ras.ci, main="")
    plot(humdata, add=T, lwd=0.5)
    dev.off()
    
    writeRaster(ras.median, paste0("livestock_rasters/", country, "/pig_median_", year), format="GTiff", overwrite=TRUE)
    writeRaster(ras.ci, paste0("livestock_rasters/", country, "/pig_95ci_", year), format="GTiff", overwrite=TRUE)
    print(year)
  }
}
#--------------------------------#
# Cases
#--------------------------------#
if(country=="Malawi"){
  png(paste0("plots/", country, "/cases_map_", country.short, ".png"), height=height, width=width)
  plot(humdata)
  plot(HAT_data, opacity=0.005, pch=19, cex=(HAT_data@data$Nw_HAT)/2, col="red", add=TRUE )
  dev.off()
}

if(country=="DRC"){
  png(paste0("plots/", country, "/cases_map_", country.short, ".png"), height=height, width=width)
  plot(humdata)
  plot(HAT_data, opacity=0.005, pch=19, cex=(HAT_data@data$Nw_HAT)/10, col="red", add=TRUE )
  dev.off()
}

if(country=="Uganda"){
  type<-c("g", "r")
  for(i in 1:2){
    types<-unique(HAT_data@data$Parasit)
    data<-HAT_data[which(HAT_data@data$Parasit==types[i]),]
    png(paste0("plots/", country, "/cases_map_", country.short, "_", type[i], ".png"), height=height, width=width)
    plot(humdata)
    plot(data, opacity=0.005, pch=19, cex=(data@data$Nw_HAT_)/10, col="red", add=TRUE )
    dev.off()
  }
}

if(country=="South_Sudan"){
  vars<-paste0("cs", 2000:2014)
  HAT@data$cases<-(rowSums(HAT@data[,c(vars)]))
  
  med.palette <- brewer.pal(n = 7, name = "Reds")
  med.int <- classIntervals(HAT@data$cases,
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/cases_map_", country.short, ".png"), height=8, width=11, units="in", res=300)
  plot(census)
  plot(HAT, border = T, col = med.col,
       axes = F, main="", add=TRUE)
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
}

#--------------------------------#
# Raster plot for NDVI (2010)
#--------------------------------#
NDVI <- raster("Data/Mediator/LTDR5_2010.tif")
NDVI.crop<-crop(NDVI, districts)
NDVI.sub<-mask(NDVI.crop, districts)

if(country=="Uganda"){
  NDVI.sub<-mask(NDVI.sub, victoria, inverse=TRUE)
}

png(paste0("plots/", country, "/NDVI_", country.short, ".png"), height=height, width=width)
plot(NDVI.sub)
plot(humdata, add=TRUE)
dev.off()

if(country=="South_Sudan"){
  med.palette <- brewer.pal(n = 7, name = "BrBG")
  med.int <- classIntervals(round(HAT@data$NDVI2008,2),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  
  png(paste0("plots/", country, "/NDVI_", country.short, ".png"), height=8, width=20, units="in", res=300)
  plot(HAT, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
}

#--------------------------------#
# Raster plot for LST (2010)
#--------------------------------#
LSTs<-list()
for(i in 1:length(modisIdx)){
  filenames<-list.files('Data/Mediator/MODIS/original/')
  filenames.short<-str_sub(filenames, start=10, end= 13)#Just do time in years
  search<-which(filenames.short==2010)
  filenames.year<-filenames[search]
  use<-filenames.year[which(grepl(modisIdx[i], filenames.year)==TRUE)][1]
  sds <- get_subdatasets(paste0('Data/Mediator/MODIS/original/', use))
  filename <- paste0('Data/Mediator/MOD21_',year,'_', i, '.tif')
  gdal_translate(sds[1], dst_dataset = filename)
  res <- raster(filename)
  LSTs<-c(LSTs, res)
}

LST<-Reduce(merge, LSTs)
LST<-projectRaster(LST, crs=crs(humdata))
LST.crop<-crop(LST, humdata)
LST.sub<-mask(LST.crop, humdata)
if(country=="Uganda"){
  LST.sub<-mask(LST.sub, victoria, inverse=TRUE)
}

png(paste0("plots/", country, "/LST_", country.short, ".png"), height=height, width=height)
plot(LST.sub)
plot(humdata, add=TRUE)
dev.off()

if(country=="South_Sudan"){
  med.palette <- brewer.pal(n = 7, name = "RdBu")
  med.palette<-rev(med.palette)
  med.int <- classIntervals(round(HAT@data$LST_2008,2),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  
  png(paste0("plots/", country, "/LST_", country.short, ".png"), height=8, width=20, units="in", res=300)
  plot(HAT, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
}

#--------------------------------#
# Raster plot for elevation
#--------------------------------#
if(country=="Malawi"){
  gmted3 <- raster("Data/Predictor_data/Elevation/30s030e_20101117_gmted_med075.tif") 
  gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
  gmted<-merge(gmted3, gmted6)
  gmted.crop<-crop(gmted, extent(humdata))
  gmted.sub<-mask(gmted.crop, humdata)
}

if(country=="Uganda"){
  gmted2 <- raster("Data/Predictor_data/Elevation/10s000e_20101117_gmted_med075.tif") 
  gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
  gmted<-merge(gmted2, gmted6)
  gmted.crop<-crop(gmted, extent(humdata))
  gmted.sub<-mask(gmted.crop, humdata)
  gmted.sub<-mask(gmted.sub, victoria, inverse=TRUE)
}

if(country=="DRC"){
  gmted1 <- raster("Data/Predictor_data/10n000e_20101117_gmted_med075.tif")
  gmted2 <- raster("Data/Predictor_data/10s000e_20101117_gmted_med075.tif") 
  gmted3 <- raster("Data/Predictor_data/30s030e_20101117_gmted_med075.tif") 
  gmted4 <- raster("Data/Predictor_data/30s000e_20101117_gmted_med075.tif") 
  gmted5 <- raster("Data/Predictor_data/Elevation/10n030e_20101117_gmted_med075.tif") 
  gmted6 <- raster("Data/Predictor_data/Elevation/10s030e_20101117_gmted_med075.tif") 
  gmted<-merge(gmted1, gmted2, gmted3, gmted4, gmted5, gmted6)
  gmted.crop<-crop(gmted, extent(humdata))
  gmted.sub<-mask(gmted.crop, humdata)
}

png(paste0("plots/", country, "/elevation_", country.short, ".png"), height=height, width=width)
plot(gmted.sub)
plot(humdata, add=TRUE)
dev.off()

if(country=="South_Sudan"){
  med.palette <- brewer.pal(n = 7, name = "Greys")
  med.int <- classIntervals(round(HAT@data$LST_2008,2),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/elevation_", country.short, ".png"), width=20, height=8, units="in", res=300)
  plot(HAT, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
}
#--------------------------------#
# Raster plot for wealth 
#--------------------------------#
if(country=="South_Sudan"){
  med.palette <- brewer.pal(n = 7, name = "Greens")
  med.int <- classIntervals(round(HAT@data$wealth,2),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/wealth_", country.short, ".png"), width=20, height=8, units="in", res=300)
  plot(HAT, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
  
  #Wealth SE
  med.palette <- brewer.pal(n = 7, name = "Greens")
  med.int <- classIntervals(round(HAT@data$wealth_se,4),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  png(paste0("plots/", country, "/wealth_se_", country.short, ".png"), width=20, height=8, units="in", res=300)
  plot(HAT, border = T, col = med.col,
       axes = F, main="")
  legend(x = "right",inset = 0,
         legend = names(attr(med.col, 'table')),
         fill = med.palette, cex= 1.25, horiz = FALSE, bty = 'n')
  dev.off()
}else{
  m=4
  
  model<-wealth.model
  nightRaster = raster(paste('Data/Predictor_data/Nighttime_lights/nightlights_2010.tif', sep=""))
  nightlights.crop<-crop(nightRaster, data.pred)
  merged.trans <- spTransform(data.pred, crs(nightlights.crop))
  nightlights.extract<-raster::extract(nightlights.crop, merged.trans)
  data.pred@data$nightlights<-nightlights.extract
  data.pred@data$nightlights<-ifelse(is.na(data.pred@data$nightlights)==TRUE, round(mean(data.pred@data$nightlights, na.rm=TRUE)), data.pred@data$nightlights)
  
  wealth.predictors<-paste(model$names.fixed, collapse=" + ")
  wealth.form.cov <- as.formula(paste("~ 0 +", wealth.predictors))
  kIdx<-which(all.years==2010)
  wealthSamples = samplingGauss(res.zp = model,
                                samp.dat= wealth_data,
                                data.pred = data.pred,
                                form.cov = wealth.form.cov,
                                nSamp = 1000,
                                spde= spde_RE,
                                year = 2010,
                                kIdx= kIdx,
                                loc.pred=loc.pred,
                                prevSample=NULL)
  
  WealthRes = list(wealth.median = matrix(0, nrow = ny, ncol = nx),
                   wealth.l     = matrix(0, nrow = ny, ncol = nx),
                   wealth.u     = matrix(0, nrow = ny, ncol = nx))
  
  wealth.median = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = median)
  wealth.median = matrix(wealth.median, ncol = dim(xMat)[2])
  WealthRes$wealth.median = wealth.median
  
  wealth.upper = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = quantile, probs=0.975, na.rm=T)
  wealth.upper = matrix(wealth.upper, ncol = dim(xMat)[2])
  
  wealth.lower = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = quantile, probs=0.025, na.rm=T)
  wealth.lower = matrix(wealth.lower, ncol = dim(xMat)[2])
  
  WealthRes$wealth.lower = wealth.lower
  WealthRes$wealth.upper = wealth.upper
  
  WealthRes$CI<-WealthRes$wealth.upper - WealthRes$wealth.lower
  
  gridM = data.frame(Longitude = gridP@coords[,1],
                     Latitude = gridP@coords[,2])
  coordinates(gridM) = ~ Longitude + Latitude
  proj4string(gridM) = proj4string(formedMap)
  mask = as.numeric(over(gridM, formedMap)[,1])
  mask[!is.na(mask)] = 1
  mask = matrix(mask, ncol = nx)
  
  wealth.median<-wealth.median*mask
  
  ras.median<-raster(wealth.median, crs=proj, xmn=min(xMat), xmx=max(xMat),
                     ymn=min(yMat), ymx=max(yMat))
  ras.median<-flip(ras.median, direction="y")
  ras.median<-crop(ras.median, districts)
  
  if(country=="Uganda"){
    ras.median<-mask(ras.median, victoria, inverse=TRUE)
  }
  
  png(paste0("plots/", country, "/wealth_", country.short, "10.png"), height=height, width=width)
  plot(ras.median, main="")
  plot(humdata, add=T, lwd=0.5)
  dev.off()
  
  wealth.lower<-wealth.lower*mask
  wealth.upper<-wealth.upper*mask
  wealth.ci<-wealth.upper-wealth.lower
  
  ras.ci<-raster(wealth.ci, crs=proj, xmn=min(xMat), xmx=max(xMat),
                 ymn=min(yMat), ymx=max(yMat))
  ras.ci<-flip(ras.ci, direction="y")
  ras.ci<-crop(ras.ci, districts)
  
  if(country=="Uganda"){
    ras.ci<-mask(ras.ci, victoria, inverse=TRUE)
  }
  
  png(paste0("plots/", country, "/wealth_95ci_", country.short, "10.png"), height=height, width=width)
  plot(ras.ci, main="")
  plot(humdata, add=T, lwd=0.5)
  dev.off()  
}

#-----------------------------#
# Time plots
#-----------------------------#
quant<-function(x,bound){
  sd<-sd(x, na.rm=T)
  mean=mean(x, na.rm=T)
  if(bound=="upper"){
    res<-mean+1.96*sd
  }else{
    res<-mean-1.96*sd
  }
  return(res)
}

if(country=="South_Sudan"){
  HAT.time<-as.data.frame(matrix(nrow=length(2000:2014), ncol=13))
  colnames(HAT.time)<-c("year", "cases", "cases_low", "cases_high", "n", "n_low", "n_high",
                        "NDVI", "NDVI_low", "NDVI_high", "LST", "LST_low", "LST_high")
  HAT.time$year<-c(2000:2014)
  for(i in 1:length(2000:2014)){
    year=HAT.time$year[i]
    colIdx<-which(grepl(year, colnames(HAT@data)))
    
    varIdx<-which(grepl("cs", colnames(HAT@data)))
    caseIdx<-varIdx[which(varIdx %in% colIdx)]
    HAT.time$cases[i]<-median(HAT@data[,caseIdx], na.rm=T)
    HAT.time$cases_low[i]<-quantile(HAT@data[,caseIdx], probs=0.025)
    HAT.time$cases_high[i]<-quantile(HAT@data[,caseIdx], probs=0.975)
    
    varIdx<-which(grepl("NDVI", colnames(HAT@data)))
    caseIdx<-varIdx[which(varIdx %in% colIdx)]
    HAT.time$NDVI[i]<-median(HAT@data[,caseIdx], na.rm=T)
    HAT.time$NDVI_low[i]<-quantile(HAT@data[,caseIdx], probs=0.025)
    HAT.time$NDVI_high[i]<-quantile(HAT@data[,caseIdx], probs=0.975)
    
    varIdx<-which(grepl("LST", colnames(HAT@data)))
    caseIdx<-varIdx[which(varIdx %in% colIdx)]
    HAT.time$LST[i]<-median(HAT@data[,caseIdx],na.rm=T)
    HAT.time$LST_low[i]<-quantile(HAT@data[,caseIdx], probs=0.025)
    HAT.time$LST_high[i]<-quantile(HAT@data[,caseIdx], probs=0.975)
    
    colIdx<-which(grepl(year, colnames(active@data)))
    varIdx<-which(grepl("n2", colnames(active@data)))
    caseIdx<-varIdx[which(varIdx %in% colIdx)]
    HAT.time$n[i]<-median(active@data[,caseIdx], na.rm=T)
    HAT.time$n_low[i]<-quantile(active@data[,caseIdx], probs=0.025)
    HAT.time$n_high[i]<-quantile(active@data[,caseIdx], probs=0.975)
  }
  
  png(paste0("plots/", country, "/cases_time_", country.short, ".png"), width=20, height=8, units="in", res=300)
  plot(NA,
       xlab="Year", ylab="All cases", 
       ylim=c(0, max(HAT.time$cases_high)*1.1),
       xlim=c(min(HAT.time$year), max(HAT.time$year)),
       type='l', lwd=2)
  lines(HAT.time$year, HAT.time$cases, lwd=2)
  lines(HAT.time$year, HAT.time$cases_low, lwd=1, lty=2 )
  lines(HAT.time$year, HAT.time$cases_high, lwd=1, lty=2 )
  dev.off()
  
  png(paste0("plots/", country, "/sampled_time_", country.short, ".png"), width=20, height=8, units="in", res=300)
  plot(NA,
       xlab="Year", ylab="Number sampled by active surveillance", 
       ylim=c(0, max(HAT.time$n_high)*1.1),
       xlim=c(min(HAT.time$year), max(HAT.time$year)),
       type='l', lwd=2)
  lines(HAT.time$year, HAT.time$n, lwd=2)
  lines(HAT.time$year, HAT.time$n_low, lwd=1, lty=2 )
  lines(HAT.time$year, HAT.time$n_high, lwd=1, lty=2 )
  dev.off()
  
  png(paste0("plots/", country, "/NDVI_time_", country.short, ".png"), width=20, height=8, units="in", res=300)
  plot(NA,
       xlab="Year", ylab="NDVI", 
       ylim=c(min(HAT.time$NDVI_low, na.rm=T)*0.9, max(HAT.time$NDVI_high)*1.1),
       xlim=c(min(HAT.time$year), max(HAT.time$year)),
       type='l', lwd=2)
  lines(HAT.time$year, HAT.time$NDVI, lwd=2)
  lines(HAT.time$year, HAT.time$NDVI_low, lwd=1, lty=2 )
  lines(HAT.time$year, HAT.time$NDVI_high, lwd=1, lty=2 )
  dev.off()
  
  png(paste0("plots/", country, "/LST_time_", country.short, ".png"), width=20, height=8, units="in", res=300)
  plot(NA,
       xlab="Year", ylab="LST", 
       ylim=c(min(HAT.time$LST_low)*0.9, max(HAT.time$LST_high)*1.1),
       xlim=c(min(HAT.time$year), max(HAT.time$year)),
       type='l', lwd=2)
  lines(HAT.time$year, HAT.time$LST, lwd=2)
  lines(HAT.time$year, HAT.time$LST_low, lwd=1, lty=2 )
  lines(HAT.time$year, HAT.time$LST_high, lwd=1, lty=2 )
  dev.off()
}else{
  if(country=="Malawi"|country=="DRC"){
    #1: livestock density
    data_time = HAT@data %>%
      dplyr::group_by(year) %>% 
      dplyr::summarize(density= mean(density_med_p, na.rm=T), 
                       low_q=mean(density_q025_p, na.rm=T),
                       high_q=mean(density_q975_p, na.rm=T))%>%
      as.data.frame() 
    
    par(mfrow = c(1,1))
    png(paste0("plots/", country, "/pigs_time_", country.short, ".png"), height=height, width=width)
    plot(NA,
         xlab="Year", ylab="Density", 
         ylim=c(0, max(data_time$high_q)*1.1),
         xlim=c(min(data_time$year), max(data_time$year)),
         type='l', lwd=2)
    lines(data_time$year, data_time$density, lwd=2)
    lines(data_time$year, data_time$low_q, lwd=1, lty=2 )
    lines(data_time$year, data_time$high_q, lwd=1, lty=2 )
    dev.off()
    
    data_time = HAT@data %>%
      dplyr::group_by(year) %>% 
      dplyr::summarise(density= mean(density_med_c, na.rm=T), 
                       low_q=mean(density_q025_c, na.rm=T),
                       high_q=mean(density_q975_c, na.rm=T))%>%
      as.data.frame() 
    
    par(mfrow = c(1,1))
    png(paste0("plots/", country, "/cattle_time_", country.short, ".png"), height=height, width=width)
    plot(NA,
         xlab="Year", ylab="Density", 
         ylim=c(0, max(data_time$high_q)*1.1),
         xlim=c(min(data_time$year), max(data_time$year)),
         type='l', lwd=2)
    lines(data_time$year, data_time$density, lwd=2)
    lines(data_time$year, data_time$low_q, lwd=1, lty=2 )
    lines(data_time$year, data_time$high_q, lwd=1, lty=2 )
    dev.off()
    
    #2: Cases
    data_time = HAT@data %>%
      dplyr::group_by(year) %>% 
      dplyr::summarise(cases_mean= mean(cases, na.rm=T), 
                       cases_low=quant(cases, bound="lower"),
                       cases_high=quant(cases, bound="upper"),
                       cases_sum=sum(cases, na.rm=T))%>%
      as.data.frame() 
    
    par(mfrow = c(1,1))
    png(paste0("plots/", country, "/cases_time_mean_", country.short, ".png"), height=height, width=width)
    plot(NA,
         xlab="Year", ylab="Average cases per cluster", 
         ylim=c(min(data_time$cases_low)*0.9, max(data_time$cases_high)*1.1),
         xlim=c(min(data_time$year), max(data_time$year)),
         type='l', lwd=2)
    lines(data_time$year, data_time$cases_mean, lwd=2)
    lines(data_time$year, data_time$cases_low, lwd=1, lty=2 )
    lines(data_time$year, data_time$cases_high, lwd=1, lty=2 )
    dev.off()
    
    par(mfrow = c(1,1))
    png(paste0("plots/", country, "/cases_time_sum_", country.short, ".png"), height=height, width=width)
    plot(NA,
         xlab="Year", ylab="Total cases", 
         ylim=c(min(data_time$cases_sum)*0.9, max(data_time$cases_sum)*1.1),
         xlim=c(min(data_time$year), max(data_time$year)),
         type='l', lwd=2)
    lines(data_time$year, data_time$cases_sum, lwd=2)
    dev.off()
    
    #3: NDVI and LST
    data_time = HAT@data %>%
      dplyr::group_by(year) %>% 
      dplyr::summarise(NDVI_mean= mean(NDVI, na.rm=T), 
                       NDVI_low=quant(NDVI, bound="lower"),
                       NDVI_high=quant(NDVI, bound="upper"),
                       LST_mean= mean(LST, na.rm=T), 
                       LST_low=quant(LST, bound="lower"),
                       LST_high=quant(LST, bound="upper"))%>%
      as.data.frame() 
    
    par(mfrow = c(1,1))
    png(paste0("plots/", country, "/NDVI_time_", country.short, ".png"), height=height, width=width)
    plot(NA,
         xlab="Year", ylab="NDVI", 
         ylim=c(min(data_time$NDVI_low)*0.9, max(data_time$NDVI_high)*1.1),
         xlim=c(min(data_time$year), max(data_time$year)),
         type='l', lwd=2)
    lines(data_time$year, data_time$NDVI_mean, lwd=2)
    lines(data_time$year, data_time$NDVI_low, lwd=1, lty=2 )
    lines(data_time$year, data_time$NDVI_high, lwd=1, lty=2 )
    dev.off()
    
    par(mfrow = c(1,1))
    png(paste0("plots/", country, "/LST_time_", country.short, ".png"), height=height, width=width)
    plot(NA,
         xlab="Year", ylab="LST", 
         ylim=c(min(data_time$LST_low)*0.9, max(data_time$LST_high)*1.1),
         xlim=c(min(data_time$year), max(data_time$year)),
         type='l', lwd=2)
    lines(data_time$year, data_time$LST_mean, lwd=2)
    lines(data_time$year, data_time$LST_low, lwd=1, lty=2 )
    lines(data_time$year, data_time$LST_high, lwd=1, lty=2 )
    dev.off()
  }
  
  if(country=="Uganda"){
    for(i in 1:2){
      #1: Livestock density
      dat<-list(HAT.g, HAT.r)
      if(i==1){
        t="g"
      }else{
        t="r"
      }
      data_time = dat[[i]]@data %>%
        dplyr::group_by(year) %>% 
        dplyr::summarise(density= mean(density_med_p, na.rm=T), 
                         low_q=mean(density_q025_p, na.rm=T),
                         high_q=mean(density_q975_p, na.rm=T))%>%
        as.data.frame() 
      
      par(mfrow = c(1,1))
      png(paste0("plots/", country, "/pigs_time_", country.short, "_", t, ".png"), height=height, width=width)
      plot(NA,
           xlab="Year", ylab="Density", 
           ylim=c(0, max(data_time$high_q)*1.1),
           xlim=c(min(data_time$year), max(data_time$year)),
           type='l', lwd=2)
      lines(data_time$year, data_time$density, lwd=2)
      lines(data_time$year, data_time$low_q, lwd=1, lty=2 )
      lines(data_time$year, data_time$high_q, lwd=1, lty=2 )
      dev.off()
      
      data_time = dat[[i]]@data %>%
        dplyr::group_by(year) %>% 
        dplyr::summarise(density= mean(density_med_c, na.rm=T), 
                         low_q=mean(density_q025_c, na.rm=T),
                         high_q=mean(density_q975_c, na.rm=T))%>%
        as.data.frame() 
      
      par(mfrow = c(1,1))
      png(paste0("plots/", country, "/cattle_time_", country.short, "_", t, ".png"), height=height, width=width)
      plot(NA,
           xlab="Year", ylab="Density", 
           ylim=c(0, max(data_time$high_q)*1.1),
           xlim=c(min(data_time$year), max(data_time$year)),
           type='l', lwd=2)
      lines(data_time$year, data_time$density, lwd=2)
      lines(data_time$year, data_time$low_q, lwd=1, lty=2 )
      lines(data_time$year, data_time$high_q, lwd=1, lty=2 )
      dev.off()
      
      #2: Cases
      data_time = dat[[i]]@data %>%
        dplyr::group_by(year) %>% 
        dplyr::summarise(cases_mean= mean(cases, na.rm=T), 
                         cases_low=quant(cases, bound="lower"),
                         cases_high=quant(cases, bound="upper"),
                         cases_sum=sum(cases, na.rm=T))%>%
        as.data.frame() 
      
      par(mfrow = c(1,1))
      png(paste0("plots/", country, "/cases_time_mean_", country.short, "_", t, ".png"), height=height, width=width)
      plot(NA,
           xlab="Year", ylab="Average cases per cluster", 
           ylim=c(min(data_time$cases_low)*0.9, max(data_time$cases_high)*1.1),
           xlim=c(min(data_time$year), max(data_time$year)),
           type='l', lwd=2)
      lines(data_time$year, data_time$cases_mean, lwd=2)
      lines(data_time$year, data_time$cases_low, lwd=1, lty=2 )
      lines(data_time$year, data_time$cases_high, lwd=1, lty=2 )
      dev.off()
      
      par(mfrow = c(1,1))
      png(paste0("plots/", country, "/cases_time_sum_", country.short, "_", t, ".png"), height=height, width=width)
      plot(NA,
           xlab="Year", ylab="Total cases", 
           ylim=c(min(data_time$cases_sum)*0.9, max(data_time$cases_sum)*1.1),
           xlim=c(min(data_time$year), max(data_time$year)),
           type='l', lwd=2)
      lines(data_time$year, data_time$cases_sum, lwd=2)
      dev.off()
      
      #3: NDVI and LST
      data_time = dat[[i]]@data %>%
        dplyr::group_by(year) %>% 
        dplyr::summarise(NDVI_mean= mean(NDVI, na.rm=T), 
                         NDVI_low=quant(NDVI, bound="lower"),
                         NDVI_high=quant(NDVI, bound="upper"),
                         LST_mean= mean(LST, na.rm=T), 
                         LST_low=quant(LST, bound="lower"),
                         LST_high=quant(LST, bound="upper"))%>%
        as.data.frame() 
      
      par(mfrow = c(1,1))
      png(paste0("plots/", country, "/NDVI_time_", country.short, "_", t, ".png"), height=height, width=width)
      plot(NA,
           xlab="Year", ylab="NDVI", 
           ylim=c(min(data_time$NDVI_low)*0.9, max(data_time$NDVI_high)*1.1),
           xlim=c(min(data_time$year), max(data_time$year)),
           type='l', lwd=2)
      lines(data_time$year, data_time$NDVI_mean, lwd=2)
      lines(data_time$year, data_time$NDVI_low, lwd=1, lty=2 )
      lines(data_time$year, data_time$NDVI_high, lwd=1, lty=2 )
      dev.off()
      
      par(mfrow = c(1,1))
      png(paste0("plots/", country, "/LST_time_", country.short, "_", t, ".png"), height=height, width=width)
      plot(NA,
           xlab="Year", ylab="LST", 
           ylim=c(min(data_time$LST_low)*0.9, max(data_time$LST_high)*1.1),
           xlim=c(min(data_time$year), max(data_time$year)),
           type='l', lwd=2)
      lines(data_time$year, data_time$LST_mean, lwd=2)
      lines(data_time$year, data_time$LST_low, lwd=1, lty=2 )
      lines(data_time$year, data_time$LST_high, lwd=1, lty=2 )
      dev.off()
    }
  }
}

#------------------#
# Data availability
#------------------#
availability_cattle<-as.data.frame(matrix(NA, nrow=length(2000:2020), ncol=4))
colnames(availability_cattle)<-c("Year", "Malawi", "Uganda","DRC")
availability_cattle$Year<-2000:2020
availability_cattle$Malawi<-c(rep(0, 4),767, 0, 143, 0,0,0,1575,0, 140, 1461, 140, 0, 1954, 150, 0,0,0)
availability_cattle$Uganda<-c(rep(0,6), 336, 0,0, 394, 216, 618, 0,0, 210, 0,685, 0, 810, 0,0)
availability_cattle$DRC<-c(rep(0,10), 350, 0,0, 492, rep(0,7))

max_cattle<-max(c(max(availability_cattle$Uganda), max(availability_cattle$DRC),max(availability_cattle$Malawi)))

png(paste0("plots/Malawi/availability_cattle_mw.png"), height=400, width=400)
p<-ggplot(data=availability_cattle, aes(x=Year, y=Malawi))+geom_bar(stat="identity")+theme_bw()+ylab("Clusters")+ coord_cartesian(ylim=c(0, max_cattle))
p
dev.off()

png(paste0("plots/DRC/availability_cattle_drc.png"), height=400, width=400)
p<-ggplot(data=availability_cattle, aes(x=Year, y=DRC))+geom_bar(stat="identity")+theme_bw()+ylab("Clusters")+ coord_cartesian(ylim=c(0, max_cattle))
p
dev.off()

png(paste0("plots/Uganda/availability_cattle_ug.png"), height=400, width=400)
p<-ggplot(data=availability_cattle, aes(x=Year, y=Uganda))+geom_bar(stat="identity")+theme_bw()+ylab("Clusters")+ coord_cartesian(ylim=c(0, max_cattle))
p
dev.off()

availability_pigs<-as.data.frame(matrix(NA, nrow=length(2000:2020), ncol=4))
colnames(availability_pigs)<-c("Year", "Malawi", "Uganda","DRC")
availability_pigs$Year<-2000:2020
availability_pigs$Malawi<-c(rep(0,4), 793, 0, 142, 0,0,0,1575, 0, 137, 1461, 130, 0,1954, 150,0,0,0)
availability_pigs$Uganda<-c(rep(0,6), 336, 0,0, 394, 216, 659, 0,0, 210, 0,685, 0, 810, 0,0)
availability_pigs$DRC<-c(rep(0,10), 350, 0,0, 492, rep(0,7))

max_pigs<-max(c(max(availability_pigs$Uganda), max(availability_pigs$DRC),max(availability_pigs$Malawi)))

png(paste0("plots/Malawi/availability_pigs_mw.png"), height=400, width=400)
p<-ggplot(data=availability_pigs, aes(x=Year, y=Malawi))+geom_bar(stat="identity")+theme_bw()+ylab("Clusters")+ coord_cartesian(ylim=c(0, max_pigs))
p
dev.off()

png(paste0("plots/DRC/availability_pigs_drc.png"), height=400, width=400)
p<-ggplot(data=availability_pigs, aes(x=Year, y=DRC))+geom_bar(stat="identity")+theme_bw()+ylab("Clusters")+ coord_cartesian(ylim=c(0, max_pigs))
p
dev.off()

png(paste0("plots/Uganda/availability_pigs_ug.png"), height=400, width=400)
p<-ggplot(data=availability_pigs, aes(x=Year, y=Uganda))+geom_bar(stat="identity")+theme_bw()+ylab("Clusters")+ coord_cartesian(ylim=c(0, max_pigs))
p
dev.off()




