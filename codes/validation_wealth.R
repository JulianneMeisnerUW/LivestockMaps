#-----------------Code adapted from Geir-Arne Fuglstad----------------#

validFn <- function(surv){
  pcprior<-list(theta = list(prior = "pc.prec", params = c(3, 0.05)))
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
  library(uwIntroStats)
  library(gpclib)
  library(ggplot2)
  library(splancs)
  library(fields)
  library(maptools)
  library(parallel)
  maptools::gpclibPermit()
  
  country<-"DRC"
  if(country=="Uganda"){
  shapefile.name<-"gadm36_UGA_shp"
  layer.name<-"gadm36_UGA_4"
  UTM<-CRS("+proj=utm +zone=35N +ellps=WGS84")
  }
  if(country=="Malawi"){
  shapefile.name<-"mwi_admbnda_adm2_nso_20181016"
  layer.name<-"mwi_admbnda_adm2_nso_20181016"
  UTM<-CRS("+proj=utm +zone=36 +ellps=WGS84")
  }

  if(country=="DRC"){
  shapefile.name<-"gadm36_COD_shp"
  layer.name<-"gadm36_COD_2"
  UTM<-CRS("+proj=utm +zone=33S +ellps=WGS84")
  }
  
  source("codes/functions_prediction.R")
  
  geographic<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84-0,0,0")
  
  #--------------------#
  # Load files
  #--------------------#
  wealth<-readOGR(paste("Data/Confounder_data/", country, "/wealth", sep=""), layer="wealth")
  
   if(country=="DRC"){
    colnames(wealth@data)<-c("cluster", "survey", "coords.x1", "coords.x2", "comb.score", "year", "urban","district", "nightlights","province", "prop_u2", "time.unstruct", "time")
    }else{
    colnames(wealth@data)<-c("cluster", "survey", "coords.x1", "coords.x2", "comb.score", "year", "urban","district", "nightlights","prop_u2", "time.unstruct", "time")
    }
  
  wealth@data$survnum<-as.numeric(wealth@data$survey)
  proj=proj4string(wealth)
  
  districts<-readOGR(paste0("Data/Exposure_data/", country, "/shapefiles/", shapefile.name), layer=layer.name) 
  message("Districts shapefile loaded")
  formedMap<- districts
  formedMap<- districts
  formedMap<-spTransform(formedMap, UTM)
  formedMap2<-rgeos::gBuffer(formedMap, byid=TRUE, width=0)
  formedMap<-spTransform(formedMap2, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84-0,0,0"))
 

  all.years=2000:2020#adjust this later
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
  
  spde_RE<-inla.spde2.pcmatern(mesh, prior.range=c(0.3, 0.05), prior.sigma=c(1,0.05))
  
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
  surveys<-unique(wealth@data$survey)
  surveys<-surveys[!is.na(surveys)]
  survNum=1:length(surveys)
  
  Validation=list()
  
  message(paste0("Preparing data for survey ", surveys[surv]) )

  if(country=="Uganda"){
  region<-over(wealth, districts)
  wealth@data$region<-region$NAME_1
  }
  
  wealth_s<-wealth[which(wealth@data$survnum==surv),]
  message("Data subsetted")

  if(country=="Uganda"){
  wealth_s@data$stratum = interaction(wealth_s@data$region, wealth_s@data$urban)
  }else{
  wealth_s@data$stratum = interaction(wealth_s@data$district, wealth_s@data$urban)
  }

  message("Stratum variable created")
  numInStr = vector('numeric', nlevels(as.factor(wealth_s@data$stratum)))
  for(i in 1:nlevels(as.factor(wealth_s@data$stratum))){
    numInStr[i] = sum(as.numeric(wealth_s@data$stratum) == i, na.rm = TRUE)
  }
  dataSet = list(a = c(), b = c(), c = c(), d = c())
  for(i in 1:length(numInStr)){
    # Choose which cluster to put in training and testing data
    idx.all = which(as.numeric(wealth_s@data$stratum) == i)
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
  message("Loop 1 done")
  for(i in 1:4){
    dataSet[[i]] = wealth_s@data$cluster[dataSet[[i]]]
    dataSet[[i]] = sort(dataSet[[i]])
  }
  message("Loop 2 done")
  # Select dataset to use
  idxKeep = sample(1:4,1)
  head(dataSet[[1]])
  message("idxKeep=",idxKeep)
  
  for(i in 1:4){
    if(i == idxKeep)
      next
    idxRM = which(wealth_s@data$cluster %in% dataSet[[i]])
    wealth_test = wealth_s[-idxRM, ]
  }
  message("Loop 3 done")
  wealth_train<-wealth[which(wealth@data$cluster%ni%wealth_test@data$cluster),]
  
  coords_train<-wealth_train@coords
  
  iset <- inla.spde.make.index('i', n.spde=spde_RE$n.spde)
  message("index made")
  
  A<-inla.spde.make.A(mesh, loc=coords_train)
  message("Projection matrix made")

  hyper<-list(prec=list(prior="pc.prec", param=c(0.5, 0.01)))
  
  N=nrow(wealth_train@data)

  stack<-inla.stack(
    data=list(resp=wealth_train@data$comb.score), 
    A=list(A,1, diag(N), diag(N), diag(N), diag(N)), 
    effects=list(iset, list(intercept=rep(1, N), urban=wealth_train@data$urban, nightlights=wealth_train@data$nightlights), time=wealth_train@data$time, cluster=wealth_train@data$cluster, survey=wealth_train@data$survey, time.unstruct=wealth_train@data$time), 
    tag='est')
  
  message("Stack built")

  year=wealth_test@data$year[1]
  kIdx = which(all.years%in%year)

  n.models<-length(list.files(paste0("results/wealth/", country, "/models")))
  for(m in 1:n.models){
  file.name<-list.files(paste0("results/wealth/", country, "/models"))[m]
  load(paste0("results/wealth/", country, "/models", file.name))
  model<-res
  message(paste0("Running model for survey ", surveys[surv], " and model ", m) )
  formula<-as.formula(model$.args$formula)
  res <- try(inla (formula, data=inla.stack.data(stack), family="gaussian",  control.predictor=list(A=inla.stack.A(stack), link=1), control.inla=list(strategy='adaptive', int.strategy='eb', cmin=0),
                 control.compute=list(config = TRUE), num.threads=inla.getOption("num.threads"), control.fixed=list(prec=1, prec.intercept=1)),TRUE)#, openmp.strategy='pardiso.parallel'
      
  if(class(res)=="try-error"){
    message(paste0("Model ", m, " for survey ", surveys[surv]," exited with error"))
    }else{
    message(paste0("Completed model ", m, " for survey ", surveys[surv]))
    #######################
    #Point-level prediction
    #######################
    pred.locs=wealth_test@coords
    predictors<-paste(res$names.fixed, collapse=" + ")
    form.cov = as.formula(paste("~ 0 +", predictors))
    
    message(paste0("Doing cluster-level prediction for survey no. ", surveys[surv], " and model ", m) )

    wealth_test@data$intercept<-1:nrow(wealth_test@data)
    
    wealthSamples = samplingGauss(res.zp = res,
                           samp.dat= wealth_train,
                           data.pred = wealth_test,
                           form.cov = form.cov,
                           nSamp = nSamples,
                           spde= spde_RE,
                           kIdx=kIdx,
                           year = year,
                           debug = FALSE,
                           inla.seed = inla.seed,
                           loc.pred=pred.locs,
                           prevSample = prevSample) 
    
    pointPred<-data.frame(matrix(nrow=dim(wealth_test)[1], ncol=6))
    colnames(pointPred)=c("Median", "Mean", "SD", "Min", "Max", "Obs")
    pointPred$Median = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = median)
    pointPred$SD = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = sd)
    pointPred$Mean = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = mean)
    pointPred$Max = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = max)
    pointPred$Min = apply(X = wealthSamples$wealth.samples, MARGIN = 2, FUN = min)
    score<-wealth_test@data$comb.score
    pointPred$Obs=score
    MSE_point=mean((pointPred$Obs - pointPred$Median)^2, na.rm=T)
    message("Point level prediction done")
    
    resVal=list(survey=surveys[surv], model=paste0("Model_", m), formula=res$.args$formula, pointPred=pointPred, MSE_point=MSE_point)
    message("resVal object created")
    
    Validation=c(Validation, list(resVal))
    message("validation list updated")
  }
}
  save(file=paste0("results/wealth/", country, "/validation_wealth", surv, "_", country, ".Rdata"), Validation)
  return(Validation)
}
