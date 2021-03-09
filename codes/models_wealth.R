#-----------------Code adapted from Geir-Arne Fuglstad----------------#

modelFn <- function(m){
  
  #--------------------#
  # Load packages
  #--------------------#
  library(dplyr)
  library(tidyr)
  library(rgdal)
  library(sp)
  library(raster)
  library(foreign)
  library(SpatialEpi)
  library(survey)
  library(RColorBrewer)
  library(INLA)
  inla.setOption(pardiso.license= "~/pardiso.lic")
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
  
  results=list()
  
  #--------------------#
  # Load files
  #--------------------#
  set.seed(1154)
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

  folderData<-paste0("results/wealth/", country, "/models/")
  
  if(file.exists(folderData)==FALSE){
    dir.create(file.path(folderData))
  }
  
  wealth_data<-readOGR(paste("Data/Confounder_data/", country, "/wealth", sep=""), layer="wealth")
  if(country=="DRC"){
  colnames(wealth@data)<-c("cluster", "survey", "coords.x1", "coords.x2", "comb.score", "year", "urban","district", "nightlights","province", "prop_u2", "time.unstruct", "time")
  }else{
  colnames(wealth@data)<-c("cluster", "survey", "coords.x1", "coords.x2", "comb.score", "year", "urban","district", "nightlights","prop_u2", "time.unstruct", "time")
  }
  proj<-proj4string(wealth)
  maptools::gpclibPermit()
  
  wealth@data$prop_g2<-1-wealth@data$prop_u2 
  
  districts<-readOGR(paste0("Data/Exposure_data/", country, "/shapefiles/", shapefile.name), layer=layer.name) 
  formedMap<- districts
  formedMap<- districts
  formedMap<-spTransform(formedMap, UTM)
  formedMap2<-rgeos::gBuffer(formedMap, byid=TRUE, width=0)
  formedMap<-spTransform(formedMap2, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84-0,0,0"))
 
  spUnion = maptools::unionSpatialPolygons(formedMap, rep(1, dim(formedMap)[1]))
  loc.mesh = matrix(0, nrow = 0, ncol = 2)
  for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
    loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
  }
  
  mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                      max.edge = c(0.3,0.6), 
                      offset = c(1, 1))
  
  spde_RE<-inla.spde2.pcmatern(mesh, prior.range=c(0.3, 0.05), prior.sigma=c(1,0.05))

  hyper<-list(prec=list(prior="pc.prec", param=c(0.5, 0.01)))
  
  coords_wealth<-wealth@coords
  
  A<-inla.spde.make.A(mesh, loc=coords_wealth)
  iset <- inla.spde.make.index('i', n.spde=spde_RE$n.spde)
  
  #Set up stacks
  N=nrow(wealth@data)

  stack<-inla.stack(
      data=list(resp=wealth@data$comb.score), 
      A=list(A,1, diag(N), diag(N), diag(N), diag(N)), 
      effects=list(iset, list(intercept=rep(1, N), urban=wealth@data$urban, nightlights=wealth@data$nightlights), time=wealth@data$time, cluster=wealth@data$cluster, survey=wealth@data$survey, time.unstruct=wealth@data$time), 
      tag='est')

  stack_comb<-inla.stack(
      data=list(resp=wealth@data$prop_g2), 
      A=list(A,1, diag(N), diag(N), diag(N), diag(N)), 
      effects=list(iset, list(intercept=rep(1, N), comb.score=wealth@data$comb.score), time=wealth@data$time, cluster=wealth@data$cluster, survey=wealth@data$survey, time.unstruct=wealth@data$time), 
      tag='est')

  formulae<-list()
  
  f1 <- resp ~ 0 + intercept + f(i, model=spde_RE, group=i.group, control.group=list(model="ar1", hyper=list(rho=list(prior="pc.cor1", param=c(0.5, 0.8))))) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) + f(survey, model='iid', hyper=hyper) + f(time.unstruct, model = 'iid', hyper=hyper) 
  
  formulae<-c(formulae, f1)
  
  f2 <- resp ~ 0 + intercept + urban + f(i, model=spde_RE, group=i.group, control.group=list(model="ar1", hyper=list(rho=list(prior="pc.cor1", param=c(0.5, 0.8))))) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) + f(survey, model='iid', hyper=hyper) + f(time.unstruct, model = 'iid', hyper=hyper) 
  
  formulae<-c(formulae, f2)
  
  f3 <- resp ~ 0 + intercept + urban + nightlights + f(i, model=spde_RE, group=i.group, control.group=list(model="ar1", hyper=list(rho=list(prior="pc.cor1", param=c(0.5, 0.8))))) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) + f(survey, model='iid', hyper=hyper) + f(time.unstruct, model = 'iid', hyper=hyper) 
  
  formulae<-c(formulae, f3)

  f4 <- resp ~ 0 + intercept + urban + nightlights + f(i, model=spde_RE, group=i.group, control.group=list(model="ar1", hyper=list(rho=list(prior="pc.cor1", param=c(0.5, 0.8))))) + f(time, model = 'rw1', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) + f(survey, model='iid', hyper=hyper) + f(time.unstruct, model = 'iid', hyper=hyper) 
  
  formulae<-c(formulae, f4)

  f5 <- resp ~ 0 + intercept + comb.score + f(i, model=spde_RE) + f(time, model = 'rw2') + f(cluster, model="iid", param=c(0.5, 5e-4)) + f(survey, model='iid', param=c(0.5, 5e-4)) + f(time.unstruct, model = 'iid', param=c(0.5, 5e-4))
  
  formulae<-c(formulae, f5)
  
  formula<-formulae[[m]]
  
  if (m==5){
   res <- inla (formula, data=inla.stack.data(stack_comb), family="gaussian",  control.predictor=list(A=inla.stack.A(stack_comb), link=1), control.inla=list(strategy='adaptive', int.strategy='eb', cmin=0),
                 control.compute=list(config = TRUE), num.threads=inla.getOption("num.threads"), control.fixed=list(prec=1, prec.intercept=1))#, openmp.strategy='pardiso.parallel'
  
  WPcoef<-res$summary.fixed[2,1]
  save(file=paste0(folderData, "/wealthWP.Rdata"), WPcoef)
  
  }else{
   res <- inla (formula, data=inla.stack.data(stack), family="gaussian",  control.predictor=list(A=inla.stack.A(stack), link=1), control.inla=list(strategy='adaptive', int.strategy='eb', cmin=0),
                 control.compute=list(config = TRUE), num.threads=inla.getOption("num.threads"), control.fixed=list(prec=1, prec.intercept=1))#, openmp.strategy='pardiso.parallel'
  
  save(file=paste0(folderData, "/wealth", m, ".Rdata"), res)
  }
  
  message(paste0("Finished model ", m))
  
  return(res)
}
