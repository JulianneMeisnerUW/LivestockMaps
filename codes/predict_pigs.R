#-----------------Code adapted from Geir-Arne Fuglstad----------------#

PredictFn <- function(m){
  #--------------------#
  # Load packages
  #--------------------#
  library(rgdal)
  library(sp)
  library(raster)
  library(foreign)
  library(SpatialEpi)
  library(survey)
  library(RColorBrewer)
  library(INLA)
  inla.setOption(mkl=TRUE)
  inla.setOption("num.threads", 8)
  #INLA:::inla.dynload.workaround()
  library(uwIntroStats)
  library(gpclib)
  library(ggplot2)
  library(splancs)
  library(fields)
  library(maptools)
  library(parallel)
  
  geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  source("codes/functions_prediction.R")
  #--------------------#
  # Load files
  #--------------------#
  set.seed(1154)
  country<-"Malawi"
  
  if(country=="Malawi"){
	country.short<-"mw"
	shapefile.name<-"mwi_admbnda_adm2_nso_20181016"
	layer.name<-"mwi_admbnda_adm2_nso_20181016"
  }
  if(country=="Uganda"){
	country.short<-"ug"
  	shapefile.name<-"gadm36_UGA_shp"
  	layer.name<-"gadm36_UGA_4"
  }

 if(country=="DRC"){
  country.short<-"drc"
  shapefile.name<-"gadm36_COD_shp"
  layer.name<-"gadm36_COD_2"
 }
  
  directory<-paste0("results/livestock_mapping/", country, "/raster_output/")
  if(file.exists(directory)==FALSE){
    dir.create(file.path(directory))
  }
  
  folderData<-paste0("results/", country, "/pig_results/models/")
  
  pig_data<-readOGR(paste("Data/Exposure_data/", country, "/Created_datasets/pig_data", sep=""), layer="pig_data")
  message("pig data loaded")
  colnames(pig_data@data)<-c("n.pigs", "mem.pigs", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")
  pig_data@data<-pig_data@data[,c("n.pigs", "mem.pigs", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")]
  proj=proj4string(pig_data)
  coords_pig<-pig_data@coords
  load(paste0(folderData, "pig_results_", country.short,".Rdata"))
  maptools::gpclibPermit()

  districts<-readOGR(paste0("Data/Exposure_data/", country, "/shapefiles/", shapefile.name), layer=layer.name) 
  
  formedMap<- districts
  if(country=="Uganda"){
	formedMap<-spTransform(formedMap, CRS("+proj=utm +zone=35N +ellps=WGS84"))
	formedMap2<-rgeos::gBuffer(formedMap, byid=TRUE, width=0)
	formedMap<-spTransform(formedMap2, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84-0,0,0"))
	}
 
  spUnion = maptools::unionSpatialPolygons(formedMap, rep(1, dim(formedMap)[1]))
  loc.mesh = matrix(0, nrow = 0, ncol = 2)
  for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
    loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
  }
 
  mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                      max.edge = c(0.3,0.6), 
                      offset = c(1, 1))

  spde_RE<-inla.spde2.pcmatern(mesh, prior.range=c(0.3, 0.05), prior.sigma=c(1,0.05))
  
  gridP<-readOGR(paste("Data/Exposure_data/", country, "Created_datasets/prediction_data", sep=""), layer = "prediction_data")
  
  mapExt = extent(formedMap)
  if(country=="Malawi"){
  nx = 200
  ny = 450
  }
  if(country=="Uganda"){
  nx = 330
  ny = 330
  }
  if(country=="DRC"){
  nx = 1170
  ny = 1090
  }

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
  
  name=paste("Model ", m, sep="")
  model<-ST.results.pigs[[m]]
  if(class(model)!="try-error"){
  
  all.years=2000:2020
  
  prediction=list()
  district_res=list()
  region_res=list()
  
  for (idxP in 1:4){ #Separate model fit for 4 x 5 year periods
    years = (all.years[1] + (idxP-1)*5 + (1:5))-1
    prevSample = NULL
    idx<-which(years<=2020)
    years<-years[idx]
    for(year in years){
      kIdx = which(all.years%in%year)
      message(paste0("Year index ", kIdx))
      print(paste('Year:', year, '/ Model:', name))
      
      maxyear=2015 
      if (year > maxyear){
        urbyear = maxyear
      }else{
        urbyear=year
      }
      colnames(gridP@data)<-c("protected", "water.body", "elevation")
      
      # Make prediction data
      urbanNum = getUrbanicity(year=urbyear, loc=loc.pred, regMap=districts, country=country)
      urbanFac = factor(x = urbanNum,
                        levels = c(0, 1))
      data.pred<-gridP
      data.pred$urban = urbanFac
      data.pred$year = rep(year, nrow(data.pred))
      data.pred@data$intercept<-rep(1, nrow(data.pred))
      data.pred@data$lon = coordinates(gridP)[,1]
      data.pred@data$lat = coordinates(gridP)[,2]
      
      predictors<-paste(model$names.fixed, collapse=" + ")
      
      form.cov <- as.formula(paste("~ 0 +", predictors))
      
      message(paste0("Prediction for model ", m, ", year ", year))
      
      liveSamples = sampling(res.zp = model,
                             samp.dat= pig_data,
                             data.pred = data.pred,
                             form.cov = form.cov,
                             nSamp = nSamples,
                             spde= spde_RE,
                             kIdx= kIdx,
                             year = unique(data.pred@data$year),
                             debug = FALSE,
                             inla.seed = inla.seed,
                             loc.pred=loc.pred,
                             prevSample = prevSample)
      
      message("prediction finished")
      # Store sample for later use
      prevSample = liveSamples$sample
      
      # Post-process to only store median and standard deviations
      LiveRes = list(live.median = matrix(0, nrow = ny, ncol = nx),
                     live.sd     = matrix(0, nrow = ny, ncol = nx))
      
      live.median = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = median, na.rm=TRUE)
      live.median = matrix(live.median, ncol = dim(xMat)[2])
      LiveRes$live.median = live.median
      
      live.sd = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = sd, na.rm=TRUE)
      live.sd = matrix(live.sd, ncol = dim(xMat)[2])
      LiveRes$live.sd = live.sd
      
      live.q025 =apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs = c(0.025), na.rm=TRUE)
      live.q025 =matrix(live.q025, ncol = dim(xMat)[2])
      LiveRes$live.q025<-live.q025

	live.q975 =apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = quantile, probs = c(0.975), na.rm=TRUE)
      live.q975 =matrix(live.q975, ncol = dim(xMat)[2])
      LiveRes$live.q975<-live.q975

      # Create mask to remove stuff outside national borders
      gridM = data.frame(Longitude = gridP@coords[,1],
                         Latitude = gridP@coords[,2])
      coordinates(gridM) = ~ Longitude + Latitude
      proj4string(gridM) = proj4string(formedMap)
      mask = as.numeric(over(gridM, formedMap)[,1])
      mask[!is.na(mask)] = 1
      mask = matrix(mask, ncol = nx)
      
      # Remove stuff
      LiveRes$live.median = LiveRes$live.median*mask
      LiveRes$live.sd = LiveRes$live.sd*mask
      LiveRes$live.q025 = LiveRes$live.q025*mask
      LiveRes$live.q975 = LiveRes$live.q975*mask
      
      message("Locations outside national borders removed")
      # Add grid locations
      LiveRes$xMat = xMat
      LiveRes$yMat = yMat
      
      # Store results
      prediction=c(prediction, LiveRes)
      
      #Write raster
      ras.median<-raster(LiveRes$live.median, crs=geographic, xmn=min(LiveRes$xMat), xmx=max(LiveRes$xMat),
                         ymn=min(LiveRes$yMat), ymx=max(LiveRes$yMat))
      ras.median<-flip(ras.median, direction="y")
      
      writeRaster(ras.median, paste0(directory, "PigPrediction_", year, "_Model_", m), format="GTiff", overwrite=TRUE)
      message("Raster done")
      
     }
  }
  res=list(prediction=prediction)
  }else{
res = "error"
}
 return(res)
}



