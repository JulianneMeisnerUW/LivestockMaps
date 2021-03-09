#------------Code adapted from Yunhan Wu----------------#
rm(list = ls())
library(rgdal)
library(spdep)
library(geosphere)
library(raster)

country<-"Malawi"

folderData<-paste0("Data/Exposure_data/", country, "/Created_datasets/Urbanicity/")

if(country=="Malawi"){
  s.frames<-c(1998, 2008)
  country.short<-"mwi"
  country.gadm<-"MWI"
  UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
  humdata.code<-"nso_20181016"
}

if(country=="Uganda"){
  s.frames<-c(1991, 2002, 2014)
  country.short<-"uga"
  country.gadm<-"UGA"
  UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
  humdata.code<-"ubos_20200824"
}

if(country=="DRC"){
  s.frames<-c(2003, 2010)
  country.short<-"cod"
  country.gadm<-"COD"
  UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
  humdata.code<-"rgc_20190911"
}

## load data
dat<-readRDS(paste0(folderData, "/", country.short, '_crc.rds'))

if(country=="Malawi"){
  dat<-dat[dat$adm1!='Likoma',] # exclude Likoma
}

#load validation data
valid<-readRDS(paste0(folderData, "/validation.rds"))
l_formula<-valid$l_formula
c_list<-valid$c_list

#load shapefile
shapefile<-readOGR(paste0("Data/Predictor_data/Urbanicity/", country, "/shapeFiles_humdata"), layer=paste0(country.short, "_admbnda_adm2_", humdata.code)) 

################################################################
#########   fit model and predict on national grid
################################################################
natl_pred<-function(obs_dat,natl_grid,l_formula,c_list){
  
  # keep only the used covariates
  pop_vec<-natl_grid$pop_den
  samp_dat<-obs_dat[,c('urban',c_list)]
  natl_dat<-natl_grid[,c(c_list,'complete')]
  
  # fit national model
  natl_fitted <-glm(l_formula,
                    data = samp_dat, family = "binomial")
  
  # predict indicator surface
  pred_prob<-predict(natl_fitted, natl_dat, na.action = na.pass,
                     type='response')
  
  pred_label<-ifelse(pred_prob > 0.5, "urban", "rural")
  
  natl_dat$pred_prob<-pred_prob
  natl_dat$pred_indicator<-as.numeric(pred_label=='urban')
  
  natl_dat$pred_prob<-ifelse(natl_dat$complete=="FALSE", NA, natl_dat$pred_prob)
  natl_dat$pred_indicator<-ifelse(natl_dat$complete=="FALSE", NA, natl_dat$pred_indicator)
  
  ind_urb_frac<-sum(pop_vec*natl_dat$pred_indicator,na.rm=TRUE)/
    sum(pop_vec,na.rm=TRUE)
  
  prob_urb_frac<-sum(pop_vec*natl_dat$pred_prob,na.rm=TRUE)/
    sum(pop_vec,na.rm=TRUE)
  
  res_list<-list()
  
  {
    res_list$ind_urb_frac<-ind_urb_frac
    res_list$prob_urb_frac<-prob_urb_frac
    res_list$fit_dat<-natl_dat
    res_list$fit_glm<-natl_fitted
  }
  
  return(res_list)
  
}

## load national grid
for(i in 1:length(s.frames)){
  if(i<length(s.frames)){
    subset<-dat[which(dat$year>=s.frames[i]&dat$year<s.frames[i+1]),]
  }else{
    subset<-dat[which(dat$year>=s.frames[i]),]
  }
  
  urb_dat<-readRDS(paste0(folderData, '/nat_grid_', s.frames[i], '.rds'))
  if(country=="Malawi"){
    urb_dat<-urb_dat[urb_dat$adm1!='Likoma',] # exclude Likoma
  }
  res_main<-natl_pred(obs_dat=subset,
                      natl_grid=urb_dat,
                      l_formula=l_formula,
                      c_list=c_list)
  ################################################################
  #########   produce urban/rural indicator surface
  ################################################################
  if(s.frames[i]<2000){
    wp.year<-2000
  }else{
    wp.year<-s.frames[i]
  }
  
  # indicator surface
  results<-res_main$fit_dat[!is.na(res_main$fit_dat$x),]
  rast.sp <- SpatialPointsDataFrame(results[c('x','y')], data=results[c('pred_indicator')], proj4string = crs(shapefile))
  ind_surf<-rasterFromXYZ(rast.sp, crs=crs(rast.sp))
  
  writeRaster(ind_surf, paste0(folderData, "/urbanicity_", s.frames[i]), format="GTiff", overwrite=TRUE)
}