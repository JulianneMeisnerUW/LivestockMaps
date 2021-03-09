#------------Code adapted from Yunhan Wu----------------#
rm(list = ls())
library(rgdal)
library(spdep)
library(geosphere)
library(raster)
library(dplyr)

'%ni%' <- Negate('%in%')

country<-"Uganda"

folderData<-paste0("Data/Exposure_data/", country, "/Created_datasets/Urbanicity/")

if(country=="Malawi"){
  s.frames<-c(1998, 2008)
  country.short<-"mwi"
  country.gadm<-"MWI"
  UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
}

if(country=="Uganda"){
  s.frames<-c(1991, 2002, 2014)
  country.short<-"uga"
  country.gadm<-"UGA"
  UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
}

if(country=="DRC"){
  s.frames<-c(2003, 2010)
  country.short<-"cod"
  country.gadm<-"COD"
  UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
}

## load data
dat<-readRDS(paste0(folderData, "/", country.short, '_crc.rds'))

if(country=="Malawi"){
  dat<-dat[dat$adm1!='Likoma',] # exclude Likoma
}

################################################################
#########   fit model and predict on national grid
################################################################
## Splits the data into a test/training set to check urbanicity model
validate<-function(obs_dat, l_formula,c_list){
  
  #subset obs_dat
  test_rows<-sample(1:nrow(obs_dat), nrow(obs_dat)/3)
  train_rows<-which(1:nrow(obs_dat)%ni%test_rows)  
  test<-obs_dat[test_rows,]
  train<-obs_dat[train_rows,]
  
  # fit model
  model <-glm(l_formula,
                    data = train, family = "binomial")
    
  # predict indicator surface
  pred_prob<-predict(model, test, na.action = na.pass,
                     type='response')
  
  pred_label<-ifelse(pred_prob > 0.5, "urban", "rural")
  pred_bin<-ifelse(pred_prob > 0.5, 1, 0)
  
  test$predicted<-pred_bin
  
  #Sensitivity (Pr(pred urban | true urban))
  test_urban<-test[which(test$urban==1),]
  n.right.urban<-as.numeric(table(test_urban$predicted)[2])
  sens<-(n.right.urban/sum(table(test_urban$predicted)))*100
  print(paste0("Sensitivity: ", sens))
  
  #Specificity (Pr(pred rural | true rural))
  test_rural<-test[which(test$urban==0),]
  n.right.rural<-as.numeric(table(test_rural$predicted)[1])
  spec<-(n.right.rural/sum(table(test_rural$predicted)))*100
  print(paste0("Specificity: ", spec))
  
  #Proportion correct
  test$dif<-abs(as.numeric(as.character(test$urban))-test$predicted)
  n.right<-as.numeric(table(test$dif)[1])
  prop.right<-(n.right/sum(table(test$dif)))*100
  print(paste0("Validity: ", prop.right))
  
  res_list=c(list(specificity=spec, sensitivity=sens, proportion_right=prop.right, l_formula=l_formula, c_list=c_list))
    
  return(res_list)
}

cov_list<-c('adm0.5', 'adm1', 'pop_den', 'noaa', 'x', 'y')

f1<-as.formula(paste('urban~pop_den+noaa+adm0.5+adm0.5*pop_den'))#Rank-deficient error when I include interaction between adm0.5 and pop dens and adm0.5 and noaa, DRC only

res_main<-validate(obs_dat=dat,
                      l_formula=f1,
                      c_list=cov_list)
  
saveRDS(res_main, file = paste0("Data/Exposure_data", country, "Created_datasets/Urbanicity/validation.rds"))


