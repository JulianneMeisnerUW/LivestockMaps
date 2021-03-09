library(INLA)
library(rgdal)
library(sp)
library(splines)
library(raster)
library(rgeos)
library(dplyr)
country<-"Uganda"

mappingFolderData<-paste0("results/livestock_mapping/", country, "/")
wealthFolderData<-paste0("results/wealth/", country, "/")
folderData<-paste0("results/regression/", country, "/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

if(country=="Malawi"){
    UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
    shapefile.name<-"mwi_admbnda_adm2_nso_20181016"
    layer.name<-"mwi_admbnda_adm2_nso_20181016"
    country.short<-"mw"
    wp.short<-"mwi"
}
if(country=="Uganda"){
    UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
    shapefile.name<-"gadm36_UGA_shp"
    layer.name<-"gadm36_UGA_4"
    country.short<-"ug"
    wp.short<-"uga"
}

source('codes/functions_prediction.R')
inla.setOption("num.threads", 4)
'%ni%'<-Negate('%in%')

#----------------------#
# Load models
#----------------------#
m=6
load(paste0(mappingFolderData,"cattle_results/models/cattle", m, ".Rdata"))
model<-res

w=4
load(paste0(wealthFolderData,"wealth", w, ".Rdata"))
wealth.model<-res

load(paste0(mappingFolderData,"prediction_cattle_", country.short, ".Rdata"))
load(paste0(mappingFolderData,"prediction_pigs_", country.short, ".Rdata"))

#----------------------#
# Population at risk
#----------------------#
cattle_data<-readOGR(paste("Data/Exposure_data/", country, "/Created_datasets/cattle_data", sep=""), layer="cattle_data")
colnames(cattle_data@data)<-c("n.cattle", "mem.cattle", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")
pig_data<-readOGR(paste("Data/Exposure_data/", country, "/Created_datasets/pig_data", sep=""), layer="pig_data")
colnames(pig_data@data)<-c( "n.pigs", "mem.pigs", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")
wealth_data<-readOGR(paste("Data/Confounder_data/", country, "/wealth", sep=""), layer="wealth")
colnames(wealth_data@data)<-c("cluster", "survey", "coords.x1", "coords.x2", "comb.score", "year", "urban", "district", "nightlights", "prop_u2", "time.unstruct", "time")

districts<-readOGR(paste0("Data/Exposure_data/", country, "/shapefiles/", shapefile.name), layer=layer.name) 

HAT<-readOGR(paste0("Data/Outcome_data/HAT_Data/", country), layer="HAT_data")
if(country=="Malawi"){
    colnames(HAT@data)<-c("district", "TA", "location_name", "lon", "lat", "focus", "year", "n.screened", "surv.type", "cases", "S1", "S2", "Sna", "type", "exported_transboundary", "protected", "elevation", "water.body")
}
if(country=="Uganda"){
    colnames(HAT@data)<-c("district", "Parish", "location_name", "lon", "lat", "focus", "year", "n.screened", "surv.type", "cases", "S1", "S2", "Sna", "type", "exported_transboundary", "protected", "elevation", "water.body")
}

HAT@data<-HAT@data[,c("district", "location_name", "focus", "year", "cases", "protected","elevation", "water.body")]

HATyears<-sort(unique(HAT@data$year))
HAT@data$id<-1:nrow(HAT@data)
HAT@data$urban<-rep(NA, nrow(HAT@data))

#Add urban variable
for (year in HATyears){
  subset<-HAT[which(HAT@data$year==year),]
  urbanNum = getUrbanicity(year=year, loc=subset@coords, regMap=districts, country=country)
  urbanFac = factor(x = urbanNum,
                    levels = c(0, 1))
  subset$urban<-urbanFac
  subsetID<-subset$id
  HAT@data$urban[subsetID]<-subset$urban
}

all.years<-2000:2020

spUnion = maptools::unionSpatialPolygons(districts, rep(1, dim(districts)[1]))
loc.mesh = matrix(0, nrow = 0, ncol = 2)
for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
  loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
}

mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                    max.edge = c(0.3,0.6), 
                    offset = c(1, 1))

spde_RE<-inla.spde2.pcmatern(mesh, prior.range=c(0.3, 0.05), prior.sigma=c(1,0.05))

pc.u=1
pc.alpha=0.01

hyperpc1<-list(prec=list(prior="pc.prec", param=c(pc.u, pc.alpha))) 

regressDat=list()

prevSample.cattle=prevSample.pigs=prevSample.wealth=NULL

for(year in HATyears){
  #Get gridded density
  model_name<-paste0("Model_",m)
  period=c(rep(1:4,each=5),5) 
  idx=which(all.years%in%year)
  pIdx=period[idx]
  ###############Cattle predictions##########
  ###########################################
  medIdx = (idx-1)*6 + 1
  sdIdx = (idx-1)*6 + 2
  lowerIdx = (idx-1)*6 + 3
  upperIdx = (idx-1)*6 + 4
  xMatIdx = (idx-1)*6 + 8
  yMatIdx = (idx-1)*6 + 9
  
  median<-Predict_nST.results.cattle[[m]]$prediction[[medIdx]]
  sd<-Predict_nST.results.cattle[[m]]$prediction[[sdIdx]]
  lower<-Predict_nST.results.cattle[[m]]$prediction[[lowerIdx]]
  upper<-Predict_nST.results.cattle[[m]]$prediction[[upperIdx]]
  xMat<-Predict_nST.results.cattle[[m]]$prediction[[xMatIdx]]
  yMat<-Predict_nST.results.cattle[[m]]$prediction[[yMatIdx]]
  
  ras.median<-raster(median, crs=geographic, xmn=min(xMat), xmx=max(xMat),
                     ymn=min(yMat), ymx=max(yMat))
  ras.median<-flip(ras.median, direction="y")
  ras.median<-crop(ras.median, districts)
  
  density.dat<-rasterToPoints(ras.median, fun=NULL, spatial=TRUE)
  
  #Check everything looks okay
  par(mfrow = c(1,1))
  pdf(file=paste0("results/regression/", country, "/Plot_Check_Cattle_",year, ".pdf"))
  plot(ras.median)
  HAT_subset<-HAT[which(HAT@data$year==year),]
  plot(HAT, add=TRUE)
  dev.off()
  
  colnames(density.dat@data)<-c("density_med_c")
  
  ras.sd<-raster(sd, crs=geographic, xmn=min(xMat), xmx=max(xMat),
                 ymn=min(yMat), ymx=max(yMat))
  ras.sd<-flip(ras.sd, direction="y")
  ras.sd<-crop(ras.sd, districts)
  density.sd<-rasterToPoints(ras.sd, fun=NULL, spatial=TRUE)
  colnames(density.sd@data)<-c("density_sd_c")
  
  ras.q025<-raster(lower, crs=geographic, xmn=min(xMat), xmx=max(xMat),
                   ymn=min(yMat), ymx=max(yMat))
  ras.q025<-flip(ras.q025, direction="y")
  ras.q025<-crop(ras.q025, districts)
  density.q025<-rasterToPoints(ras.q025, fun=NULL, spatial=TRUE)
  colnames(density.q025@data)<-c("density_q025_c")
  
  ras.q975<-raster(upper, crs=geographic, xmn=min(xMat), xmx=max(xMat),
                   ymn=min(yMat), ymx=max(yMat))
  ras.q975<-flip(ras.q975, direction="y")
  ras.q975<-crop(ras.q975, districts)
  density.q975<-rasterToPoints(ras.q975, fun=NULL, spatial=TRUE)
  colnames(density.q975@data)<-c("density_q975_c")
  
  density.dat@data$density_sd_c<-density.sd@data$density_sd
  density.dat@data$density_q025_c<-density.q025@data$density_q025
  density.dat@data$density_q975_c<-density.q975@data$density_q975
  
  merge.dist<-over(density.dat, districts)
  density.dat@data$district<-merge.dist$ADM2_EN
  density.dat@data$cases<-rep(0, nrow(density.dat@data))
  
  ###############Pig predictions############
  ###########################################
  model_name<-paste0("Model_",m)
  
  median<-Predict_nST.results.pigs[[m]]$prediction[[medIdx]]
  sd<-Predict_nST.results.pigs[[m]]$prediction[[sdIdx]]
  lower<-Predict_nST.results.pigs[[m]]$prediction[[lowerIdx]]
  upper<-Predict_nST.results.pigs[[m]]$prediction[[upperIdx]]
  xMat<-Predict_nST.results.pigs[[m]]$prediction[[xMatIdx]]
  yMat<-Predict_nST.results.pigs[[m]]$prediction[[yMatIdx]]
  
  ras.median<-raster(median, crs=geographic, xmn=min(xMat), xmx=max(xMat),
                     ymn=min(yMat), ymx=max(yMat))
  ras.median<-flip(ras.median, direction="y")
  ras.median<-crop(ras.median, districts)
  density.medp<-rasterToPoints(ras.median, fun=NULL, spatial=TRUE)
  colnames(density.medp@data)<-c("density_med_p")
  
  #Check everything looks okay
  par(mfrow = c(1,1))
  pdf(file=paste0("results/regression/", country, "/Plot_Check_Pigs_",year, ".pdf"))
  plot(ras.median)
  HAT_subset<-HAT[which(HAT@data$year==year),]
  plot(HAT, add=TRUE)
  dev.off()
  
  ras.sd<-raster(sd, crs=geographic, xmn=min(xMat), xmx=max(xMat),
                 ymn=min(yMat), ymx=max(yMat))
  ras.sd<-flip(ras.sd, direction="y")
  ras.sd<-crop(ras.sd, districts)
  density.sd<-rasterToPoints(ras.sd, fun=NULL, spatial=TRUE)
  colnames(density.sd@data)<-c("density_sd_p")
  
  ras.q025<-raster(lower, crs=geographic, xmn=min(xMat), xmx=max(xMat),ymn=min(yMat), ymx=max(yMat))
  ras.q025<-flip(ras.q025, direction="y")
  ras.q025<-crop(ras.q025, districts)
  density.q025<-rasterToPoints(ras.q025, fun=NULL, spatial=TRUE)
  colnames(density.q025@data)<-c("density_q025_p")
  
  ras.q975<-raster(upper, crs=geographic, xmn=min(xMat), xmx=max(xMat),
                   ymn=min(yMat), ymx=max(yMat))
  ras.q975<-flip(ras.q975, direction="y")
  ras.q975<-crop(ras.q975, districts)
  density.q975<-rasterToPoints(ras.q975, fun=NULL, spatial=TRUE)
  colnames(density.q975@data)<-c("density_q975_p")
  
  density.dat@data$density_med_p<-density.medp@data$density_med_p
  density.dat@data$density_sd_p<-density.sd@data$density_sd
  density.dat@data$density_q025_p<-density.q025@data$density_q025
  density.dat@data$density_q975_p<-density.q975@data$density_q975
  
  density.dat@data$focus<-rep(NA, nrow(density.dat@data))
  density.dat@data$year<-rep(year, nrow(density.dat@data))
  density.dat@data$source<-rep("pred", nrow(density.dat@data))
  
  #Predict to locations from final model
  HAT$intercept<-rep(1, nrow(HAT))
  inla.seed = as.integer(runif(1)*.Machine$integer.max)
  nSamples = 1000 
  prevSample = NULL
  predictors<-paste(model$names.fixed, collapse=" + ")
  form.cov <- as.formula(paste("~ 0 +", predictors))
  
  #Cattle first
  liveSamples = sampling(res.zp = model,
                         samp.dat= cattle_data,
                         data.pred = HAT,
                         form.cov = form.cov,
                         kIdx=idx,
                         nSamp = nSamples,
                         spde= spde_RE,
                         year = year,
                         loc.pred=HAT@coords,
			 prevSample=prevSample.cattle)
  prevSample.cattle=liveSamples$sample
  
  # Post-process to only store median and standard deviations
  LiveRes = list(live.median = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000),
                 live.sd     = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000),
                 live.q025   = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000),
                 live.q975 = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000))
  
  live.median = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = median)
  live.sd = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = sd)
  live.q025 = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs = c(0.025))
  live.q975 = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs = c(0.975))
  
  HAT.dat<-cbind(live.median, live.sd, live.q025, live.q975)
  casesIdx<-which(HAT@data$year==year & HAT@data$cases>0)
  cases<-rep(0, nrow(HAT.dat))
  cases[casesIdx]<-HAT@data$cases[casesIdx]
  HAT.dat<-cbind(HAT.dat, cases)
  HAT.dat<-as.data.frame(HAT.dat)
  colnames(HAT.dat)<-c("density_med_c", "density_sd_c", "density_q025_c", "density_q975_c", "cases")
  HAT.dat$district<-HAT@data$district
  HAT.dat$focus<-HAT@data$focus
  HAT.dat$year<-rep(year, nrow(HAT.dat))
  HAT.dat$source<-rep("WHO", nrow(HAT.dat))
  
  #Add pig predictions
  liveSamples = sampling(res.zp = model,
                         samp.dat= pig_data,
                         data.pred = HAT,
                         form.cov = form.cov,
                         nSamp = nSamples,
                         spde= spde_RE,
                         kIdx=idx,
                         year = year,
                         debug = FALSE,
                         inla.seed = inla.seed,
                         loc.pred=HAT@coords,
                         prevSample=prevSample.pigs)
  prevSample.pigs=liveSamples$sample

  # Post-process to only store median and standard deviations
  LiveRes = list(live.median = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000),
                 live.sd     = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000),
                 live.q025   = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000),
                 live.q975 = matrix(0, nrow = dim(HAT@coords)[1], ncol = 1000))
  
  live.median = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = median)
  live.sd = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = sd)
  live.q025 = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs = c(0.025))
  live.q975 = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs = c(0.975))
  
  HAT.dat$density_med_p<-live.median
  HAT.dat$density_sd_p<-live.sd
  HAT.dat$density_q025_p<-live.q025
  HAT.dat$density_q975_p<-live.q975
  
  #Label predicted points which are within 1km of observed points
  #Project in UTM first
  density.test<-SpatialPointsDataFrame(coords= density.dat@coords, data=density.dat@data, proj4string=UTM)
  HAT.test<-SpatialPointsDataFrame(coords= HAT@coords, data=HAT@data, proj4string=UTM)
  points_matrix <- gWithinDistance(HAT.test, density.test, dist = 1, byid = TRUE)
  points_matrix[lower.tri(points_matrix, diag=TRUE)] <- NA
  v <- rowSums(points_matrix, na.rm=TRUE) == 0
  
  density.dat@data$prox1km<-rep(1, nrow(density.dat))
  density.dat@data$prox1km[v]<-0
  
  HAT.dat$prox1km<-rep(0, nrow(HAT.dat))
  
  #Re-arrange columns
  HAT.dat<-HAT.dat[,c("district", "year", "focus", "cases", "density_med_c", 
                      "density_sd_c", "density_q025_c", "density_q975_c","density_med_p", 
                      "density_sd_p", "density_q025_p", "density_q975_p", "source", "prox1km")]
  density.dat@data<-density.dat@data[,c("district", "year", "focus", "cases", "density_med_c", 
                                        "density_sd_c", "density_q025_c", "density_q975_c","density_med_p", 
                                        "density_sd_p", "density_q025_p", "density_q975_p", "source", "prox1km")]
  
  Data<-rbind(density.dat@data, HAT.dat)
  
  Data$cluster<-1:nrow(Data)
  
  Data.coords<-rbind(density.dat@coords, HAT@coords)
  
  Data.sp <- SpatialPointsDataFrame(coords = Data.coords, data = Data,
                                    proj4string = geographic)

  Data.sp@data$intercept<-rep(1, nrow(Data.sp@data))
  
  #Add WorldPop
  popRasterW = raster(paste('Data/Denominator_data/WorldPop/,' wp.short, '_ppp_', year, '.tif', sep = ""))
  PW.crop<-crop(popRasterW, districts)
  trans.WP <- spTransform(Data.sp, crs(PW.crop))
  merge.WP<-raster::extract(PW.crop, trans.WP)
  Data.sp@data$WorldPop<-round(merge.WP,0) #Rounded to the nearest person
  
  #Add LandScan
  popRasterL = raster(paste('LandScan/Denominator_data/LandScan_', year, '.tif', sep = "")) 
  PL.crop<-crop(popRasterL, districts)
  crs(PL.crop)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  trans.LP <- spTransform(Data.sp, crs(PL.crop))
  merge.LP<-raster::extract(PL.crop, trans.LP)
  Data.sp@data$LandScan<-round(merge.LP,0)
  
  urbanNum<-getUrbanicity(year=year, loc=Data.sp@coords, regMap=districts, country=country)
  urbanFac = factor(x = urbanNum,
                    levels = c(0, 1))
  Data.sp@data$urban<-urbanFac

  #Add nightlights
  if(year<=2013){
  nightYear=year
  }else{
  nightYear=2013
  }
  nightRaster = raster(paste('Data/Nighttime_lights/nightlights_',nightYear,'.tif', sep=""))
  nightlights.crop<-crop(nightRaster, Data.sp)
  merged.trans <- spTransform(Data.sp, crs(nightlights.crop))
  nightlights.extract<-raster::extract(nightlights.crop, merged.trans)
  Data.sp@data$nightlights<-nightlights.extract
  Data.sp@data$nightlights<-ifelse(is.na(Data.sp@data$nightlights)==TRUE, round(mean(Data.sp@data$nightlights, na.rm=TRUE)), Data.sp@data$nightlights)
  
  #Wealth predictions
  wealth.predictors<-paste(wealth.model$names.fixed, collapse=" + ")
  wealth.form.cov <- as.formula(paste("~ 0 +", wealth.predictors))
  wealthSamples = samplingGauss(res.zp = wealth.model,
                                samp.dat= wealth_data,
                                data.pred = Data.sp,
                                form.cov = wealth.form.cov,
                                nSamp = nSamples,
                                kIdx=idx,
                                spde= spde_RE,
                                year = year,
                                loc.pred=Data.sp@coords,
				prevSample=prevSample.wealth)
  prevSample.wealth=wealthSamples$sample

  wealth.median = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = median)
  wealth.sd = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = sd)
  wealth.q025 = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = quantile, probs = c(0.025), na.rm=TRUE)
  wealth.q975 = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = quantile, probs = c(0.975), na.rm=TRUE)
  
  Data.sp@data$wealth_med<-wealth.median
  Data.sp@data$wealth_sd<-wealth.sd
  Data.sp@data$wealth_q025_p<-wealth.q025
  Data.sp@data$wealth_q975_p<-wealth.q975
  
  regressDat=c(regressDat, Data.sp)
}

save(file = paste0("results/regression/", country, "/RegressionDataMod", m, ".RData"), regressDat)

Data<-regressDat[[1]]@data
Coords<-regressDat[[1]]@coords
proj<-proj4string(regressDat[[1]])

for(i in 2:length(regressDat)){
  new_data<-regressDat[[i]]@data
  new_coords<-regressDat[[i]]@coords
  Data<-rbind(Data,new_data)
  Coords<-rbind(Coords, new_coords)
}

finalData.sp <- SpatialPointsDataFrame(coords = Coords, data = Data,
                                       proj4string = CRS(proj))

writeOGR(obj=finalData.sp, dsn=paste0("Data/Outcome_data/", country, "/HATData"), layer=paste0(country.short, "_regress"), driver="ESRI Shapefile", overwrite_layer=TRUE)
