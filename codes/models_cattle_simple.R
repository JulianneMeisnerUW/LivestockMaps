#-----------------Code adapted from Geir-Arne Fuglstad----------------#

modelFn <- function(m){
t<-Sys.time()
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
  
  #--------------------#
  # Load files
  #--------------------#
  set.seed(1154)
  country<-"Malawi"
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
  
  folderData<-paste0("results/livestock_mapping/", country, "/cattle_results/models/")
  
  if(file.exists(folderData)==FALSE){
    dir.create(file.path(folderData))
  }
  
  cattle_data<-readOGR(paste("Data/Exposure_data/", country, "/Created_datasets/cattle_data", sep=""), layer="cattle_data")
  message("cattle data read in")
  colnames(cattle_data@data)<-c("n.cattle", "mem.cattle", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")
  cattle_data@data<-cattle_data@data[,c("n.cattle", "mem.cattle", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")]
  proj<-proj4string(cattle_data)
  maptools::gpclibPermit()
  coords_cattle<-cattle_data@coords

  if(class(cattle_data@data$mem.cattle)!="numeric"){
  cattle_data@data$mem.cattle<-as.numeric(as.character(cattle_data@data$mem.cattle))
  }

  if(class(cattle_data@data$n.cattle)!="numeric"){
  cattle_data@data$n.cattle<-as.numeric(as.character(cattle_data@data$n.cattle))
  }

  cattle_data<-cattle_data[which(cattle_data@data$year>=2000),]
  table(cattle_data@data$n.cattle)

  districts<-readOGR(paste0("Data/Exposure_data/", country, "/shapefiles/", shapefile.name), layer=layer.name) 
  
  formedMap<- districts
  formedMap<-spTransform(formedMap, UTM)
  formedMap2<-rgeos::gBuffer(formedMap, byid=TRUE, width=0)
  formedMap<-spTransform(formedMap2, CRS("+proj=longlat +datum=WGS84"))
 
  spUnion = maptools::unionSpatialPolygons(formedMap, rep(1, dim(formedMap)[1]))
  loc.mesh = matrix(0, nrow = 0, ncol = 2)
  for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
    loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
  }
 
  mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                      max.edge = c(0.3,0.6), 
                      offset = c(1, 1))

  spde_RE1<-inla.spde2.pcmatern(mesh, prior.range=c(0.3,0.05), prior.sigma=c(1,0.05)) 

  spde_RE2<-inla.spde2.matern(mesh, alpha=2) 

  pcprior<-list(theta=list(prior="pc.prec", params=c(3, 0.05)))
  
  k<-length(unique(cattle_data@data$time))
  iset<-inla.spde.make.index('i', n.spde=spde_RE1$n.spde, n.group=k)
  A<-inla.spde.make.A(mesh, loc=coords_cattle, group=cattle_data@data$time)

  N<-nrow(cattle_data@data)
  
stack<-inla.stack(
    data=list(resp=cattle_data@data$n.cattle), 
    A=list(A,1, diag(N), diag(N), diag(N), diag(N)), 
    effects=list(iset, list(intercept=rep(1, N), urban=cattle_data@data$urban, protected=cattle_data@data$protected, water.body=cattle_data@data$water.body, elevation = cattle_data@data$elevation), time=cattle_data@data$time, cluster=cattle_data@data$cluster, survey=cattle_data@data$survey, time.unstruct=cattle_data@data$time), 
    tag='est')

  hyper<-list(prec=list(prior="pc.prec", param=c(0.5, 0.01)))
  
  formulae<-list()
  
  f1 <- resp ~ 0 + intercept + f(i, model=spde_RE1) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) 
  
  formulae<-c(formulae, f1)
  
  f2 <- resp ~ 0 + intercept + urban + f(i, model=spde_RE1) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) 
  
  formulae<-c(formulae, f2)
  
  f3 <- resp ~ 0 + intercept + urban + protected + f(i, model=spde_RE1) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) 
  
  formulae<-c(formulae, f3)
  
  f4 <- resp ~ 0 + intercept + urban + protected + water.body + f(i, model=spde_RE1) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) 
  
  formulae<-c(formulae, f4)
  
  f5 <- resp ~ 0 + intercept + urban + protected + water.body + elevation + f(i, model=spde_RE1) + f(time, model = 'rw2', diagonal=1e-3, scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) 
  
  formulae<-c(formulae, f5)
  
  f6 <- resp ~ 0 + intercept + urban + protected + water.body + elevation + f(i, model=spde_RE1) + f(time, model = 'rw1', scale.model=TRUE, hyper=hyper) + f(cluster, model="iid", hyper=hyper) 
  
  formulae<-c(formulae, f6)
  
  formula<-formulae[[m]]
  
  message(paste0("Running cattle model ", m))

  res <- try(inla(formula, data=inla.stack.data(stack), family="zeroinflatedpoisson1", offset=log(cattle_data@data$mem.cattle), control.predictor=list(A=inla.stack.A(stack), link=1), control.inla=list(strategy='adaptive', int.strategy='eb', cmin=0),
               control.compute=list(config = TRUE), num.threads=inla.getOption("num.threads"), control.fixed=list(prec=1, prec.intercept=1)), TRUE)#, openmp.strategy='pardiso.parallel'

if(class(res)=="try-error"){
message(paste0("Model ",m," exited with error"))
}else{
message(paste0("Completed model ", m) )
 }

  message(paste0("Finished cattle model ", m))
  save(file=paste0(folderData, "/cattle", m, ".Rdata"), res)
  
  runtime<-(Sys.time()-t)/60
  message(paste0("Runtime for model ", m, "in minutes : ", runtime))
  return(res)
}
