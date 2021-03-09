#-----------------Code adapted from Geir-Arne Fuglstad----------------#

validFn <- function(surv){
  '%ni%'<-Negate('%in%')
  #--------------------#
  # Load packages
  #--------------------#
  library(tidyr)
  library(rgdal)
  library(sp)
  library(raster)
  library(foreign)
  library(SpatialEpi)
  library(survey)
  library(RColorBrewer)
  library(INLA)
  inla.setOption(mkl=TRUE)
  inla.setOption("num.threads",8)
  #inla.setOption(pardiso.license= "~/pardiso.lic")
  library(uwIntroStats)
  library(gpclib)
  library(ggplot2)
  library(splancs)
  library(fields)
  library(maptools)
  library(parallel)
  maptools::gpclibPermit()

  message("packages loaded")
  
  country<-"Malawi"
  if(country=="Malawi"){
	shapefile.name<-"mwi_admbnda_adm2_nso_20181016"
	layer.name<-"mwi_admbnda_adm2_nso_20181016"
  }
  if(country=="Uganda"){
  	shapefile.name<-"gadm36_UGA_shp"
  	layer.name<-"gadm36_UGA_4"
  }
 if(country=="DRC"){
  shapefile.name<-"gadm36_COD_shp"
  layer.name<-"gadm36_COD_2"
  }

  source("codes/functions_prediction.R")
  
  geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84-0,0,0")
  
  #--------------------#
  # Load files
  #--------------------#
  cattle_data<-readOGR(paste("Data/Exposure_data/", country, "/Created_datasets/cattle_data", sep=""), layer="cattle_data")
  colnames(cattle_data@data)<-c("n.cattle", "mem.cattle", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")
  cattle_data@data<-cattle_data@data[,c("n.cattle", "mem.cattle", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")]
  load(paste0(folderData,"District_vector.Rdata"))
  cattle_data@data$district<-districts_cattle
  cattle_data@data$urban<-as.numeric(as.character(cattle_data@data$urban))
  cattle_data@data$density<-cattle_data@data$n.cattle/cattle_data@data$mem.cattle
  cattle_data@data$survnum<-as.numeric(cattle_data@data$survey)
  proj=proj4string(cattle_data)
  
  districts<-readOGR(paste0("Data/Exposure_data/", country, "/shapefiles/", shapefile.name), layer=layer.name) 
  formedMap<- districts
  formedMap<-spTransform(formedMap, UTM)
  formedMap2<-rgeos::gBuffer(formedMap, byid=TRUE, width=0)
  formedMap<-spTransform(formedMap2, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84-0,0,0"))
 

  gridP<-readOGR(paste("Data/Exposure_data/", country, "Created_datasets/prediction_data", sep=""), layer = "prediction_data")
  protected.crop<-readOGR(paste("Data/Exposure_data/", country, "Created_datasets/prediction_data_protected", sep=""), layer = "prediction_data_protected")
  lakes.crop<-readOGR(paste("Data/Exposure_data/", country, "Created_datasets/prediction_data_lakes", sep=""),layer = "prediction_data_lakes")
  elevation.1 <- raster(paste("Data/Exposure_data/", country, "Created_datasets/prediction_data_elevation1.tif", sep=""))
  elevation.2 <- raster(paste("Data/Exposure_data/", country, "Created_datasets/prediction_data_elevation2.tif", sep=""))
  message(paste0("Covariate data loaded for survey ", surv))
  
  #-----------------------#
  # Set-up for INLA models
  #-----------------------#
  spUnion = maptools::unionSpatialPolygons(formedMap, rep(1, dim(formedMap)[1]))
  loc.mesh = matrix(0, nrow = 0, ncol = 2)
  for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
    loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
  }
  
  mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                      max.edge = c(0.3,0.6), 
                      offset = c(1, 1))
  
  spde_RE1<-inla.spde2.pcmatern(mesh, prior.range=c(0.3, 0.05), prior.sigma=c(1,0.05))

  #-----------------------#
  # Set-up for prediction
  #-----------------------#
  nSamples = 1000 
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
  inla.seed = as.integer(runif(1)*.Machine$integer.max)
  prevSample = NULL
  
  #-----------------------#
  # Subset data
  #-----------------------#
  surveys<-unique(cattle_data@data$survey)
  surveys<-surveys[!is.na(surveys)]
  survNum=1:length(surveys)
  
  Validation=list()
  
  message(paste0("Preparing data for survey ", surv) )
  
  if(country=="Uganda"){
  region<-over(cattle_data, districts)
  cattle_data@data$region<-region$NAME_1
  }
  
  cattle_s<-cattle_data[which(cattle_data@data$survnum==surv),]
  
  if(country=="Uganda"){
  cattle_s@data$stratum = interaction(cattle_s@data$region, cattle_s@data$urban)
  }else{
  cattle_s@data$stratum = interaction(cattle_s@data$district, cattle_s@data$urban)
  }
  
  numInStr = vector('numeric', nlevels(as.factor(cattle_s@data$stratum)))
  for(i in 1:nlevels(as.factor(cattle_s@data$stratum))){
    numInStr[i] = sum(as.numeric(cattle_s@data$stratum) == i, na.rm = TRUE)
  }
  dataSet = list(a = c(), b = c(), c = c(), d = c())
  for(i in 1:length(numInStr)){
    # Choose which cluster to put in training and testing data
    idx.all = which(as.numeric(cattle_s@data$stratum) == i)
    idx.tmp = sample.int(n = numInStr[i], replace = FALSE)
    # Assign equal number to all if possible
    nEqual = floor(numInStr[i]/4)
    dataSet[[1]] = c(dataSet[[1]], idx.all[idx.tmp[1:nEqual]])
    idx.tmp = idx.tmp[-(1:nEqual)]
    dataSet[[2]] = c(dataSet[[2]], idx.all[idx.tmp[1:nEqual]])
    idx.tmp = idx.tmp[-(1:nEqual)]
    dataSet[[3]] = c(dataSet[[3]], idx.all[idx.tmp[1:nEqual]])
    idx.tmp = idx.tmp[-(1:nEqual)]
    dataSet[[4]] = c(dataSet[[4]], idx.all[idx.tmp[1:nEqual]])
    idx.tmp = idx.tmp[-(1:nEqual)]
    
    # Assign remaining at random
    if(length(idx.tmp) == 0)
      next
    idx.left = sample.int(n = 4, size = length(idx.tmp))
    for(j in 1:length(idx.left)){
      dataSet[[idx.left[j]]] = c(dataSet[[idx.left[j]]], idx.all[idx.tmp[j]])
    }
  }
  for(i in 1:4){
    dataSet[[i]] = cattle_s@data$cluster[dataSet[[i]]]
    dataSet[[i]] = sort(dataSet[[i]])
  }
  # Select dataset to use
  idxKeep = sample(1:4,1)
  
  idxRM = which(cattle_s@data$cluster %in% dataSet[[idxKeep]])
  cattle_test = cattle_s[-idxRM, ]
  message(paste0("Nrows, cattle_test: ", nrow(cattle_test@data)))
  
  cattle_train<-cattle_data[which(cattle_data@data$cluster%ni%cattle_test@data$cluster),]
  
  coords_train<-cattle_train@coords
  
  iset <- inla.spde.make.index('i', n.spde=spde_RE1$n.spde)
  
  A<-inla.spde.make.A(mesh, loc=coords_train)
  
  #Set up stack
  N= nrow(cattle_train@data)
  
  stack<-inla.stack(
      data=list(resp=cattle_train@data$n.cattle), 
      A=list(A,1, diag(N), diag(N), diag(N), diag(N)), 
      effects=list(iset, list(intercept=rep(1, N), urban=cattle_train@data$urban, protected=cattle_train@data$protected, water.body=cattle_train@data$water.body, elevation = cattle_train@data$elevation), time=cattle_train@data$time, cluster=cattle_train@data$cluster, survey=cattle_train@data$survey, time.unstruct=cattle_train@data$time), 
      tag='est')
  
  hyper<-list(prec=list(prior="pc.prec", param=c(0.5, 0.01)))
  
  message(paste0("Training data, survey ", surv))
  n.models<-length(list.files(paste0(folderData, "results/cattle_results/models")))
  for (m in 1:n.models){
  message(paste0("Running model for survey ", surv, " and model ", m) )
  
  file.name<-list.files(paste0("results/livestock_mapping/", country, "/cattle_results/models"))[m]
  load(paste0("results/livestock_mapping/", country, "/cattle_results/models", file.name))
  model<-res
  if(class(model)=="try-error"){
  resVal=list(survey=surveys[surv], model=paste0("Model_", m), result="model_error")
  Validation=c(Validation, resVal)
  }else{
  formula<-as.formula(model$.args$formula)
  
  res <- try(inla (formula, data=inla.stack.data(stack), family="zeroinflatedpoisson1", offset=log(cattle_train@data$mem.cattle), control.predictor=list(A=inla.stack.A(stack), link=1), control.inla=list(strategy='adaptive', int.strategy='eb', cmin=0),
                 control.compute=list(config = TRUE), num.threads=inla.getOption("num.threads"), control.fixed=list(prec=1, prec.intercept=1)),TRUE)#, openmp.strategy='pardiso.parallel'
  
  if(class(res)=="try-error"){
  message(paste0("Model ",m," for survey no. ", surv, " exited with error"))
  resVal=list(survey=surveys[surv], model=paste0("Model_", m), result="validation_error")
  Validation=c(Validation, resVal)
  }else{
  message(paste0("Completed model ", m, " for survey no. ", surv) )
    
  year=cattle_test@data$year[1]

  #######################
  #Point-level prediction
  #######################
  pred.locs=cattle_test@coords
  predictors<-paste(res$names.fixed, collapse=" + ")
  form.cov = as.formula(paste("~ 0 +", predictors))
  all.years<-c(2000:2020)
  kIdx = which(all.years%in%year)
  
  message(paste0("Doing cluster-level prediction for survey no. ", surv, " and model ", m) )

  liveSamples = sampling(res.zp = model,
                               samp.dat= cattle_train,
                               data.pred = cattle_test,
                               form.cov = form.cov,
                               nSamp = nSamples,
                               spde= spde_RE1,
                               kIdx=kIdx,
                               year = year,
                               debug = FALSE,
                               inla.seed = inla.seed,
                               loc.pred=pred.locs,
                               prevSample = prevSample)
  
  message("sampling done")
  
  pointPred<-data.frame(matrix(nrow=dim(cattle_test)[1], ncol=6))
  colnames(pointPred)=c("Median", "Mean", "SD", "Min", "Max", "Obs")
  pointPred$Median = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = median)
  pointPred$SD = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = sd)
  pointPred$Mean = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = mean)
  pointPred$Max = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = max)
  pointPred$Min = apply(X = liveSamples$livestock.samples, MARGIN = 2, FUN = min)
  dens<-cattle_test@data$density
  pointPred$Obs=dens
  MSE_point=mean((pointPred$Obs - pointPred$Median)^2, na.rm=T)
  message(paste0("Point level prediction done for model ", m, ", survey ", surv))
  
  resVal=list(survey=surveys[surv], model=paste0("Model_", m),pointPred=pointPred, MSE_point=MSE_point)
  Validation=c(Validation, list(resVal))
  message(paste0("Validation done for model ", m, ", survey ", surv))
  }
  }
  }
  message(paste0("Finished validation for survey number ", surv))
  save(file=paste0("results/livestock_mapping/", country, "/validation", surv, "_", country,".Rdata"), Validation)
  return(Validation)
}
