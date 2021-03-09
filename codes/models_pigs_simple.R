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
  folderData<-paste0("results/livestock_mapping/", country, "/pig_results/models/")
  
  if(file.exists(folderData)==FALSE){
    dir.create(file.path(folderData))
  }
  
  pig_data<-readOGR(paste("Data/Exposure_data/", country, "/Created_datasets/pig_data", sep=""), layer="pig_data")
  colnames(pig_data@data)<-c("n.pigs", "mem.pigs", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")
  pig_data@data<-pig_data@data[,c("n.pigs", "mem.pigs", "intercept", "urban", "protected", "water.body", "elevation", "survey", "cluster", "district", "year", "lat", "lon", "time", "time.unstruct")]
  proj<-proj4string(pig_data)
  maptools::gpclibPermit()
  coords_pig<-pig_data@coords

  pig_data<-pig_data[which(pig_data@data$year>=2000),]
  
  if(class(pig_data@data$mem.pigs)!="numeric"){
  pig_data@data$mem.pigs<-as.numeric(as.character(pig_data@data$mem.pigs))
  }

  if(class(pig_data@data$n.pigs)!="numeric"){
  pig_data@data$n.pigs<-as.numeric(as.character(pig_data@data$n.pigs))
  }
  
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
  
  k<-length(unique(pig_data@data$time))
  iset<-inla.spde.make.index('i', n.spde=spde_RE1$n.spde, n.group=k)
  A<-inla.spde.make.A(mesh, loc=coords_pig, group=pig_data@data$time)

  pig_data@data$n.pigs[which(pig_data@data$n.pigs==9509)]<-NA

  N<-nrow(pig_data@data)
  
  stack<-inla.stack(
    data=list(resp=pig_data@data$n.pigs), 
    A=list(A,1, diag(N), diag(N), diag(N), diag(N)), 
    effects=list(iset, list(intercept=rep(1, N), urban=pig_data@data$urban, protected=pig_data@data$protected, water.body=pig_data@data$water.body, elevation = pig_data@data$elevation), time=pig_data@data$time, cluster=pig_data@data$cluster, survey=pig_data@data$survey, time.unstruct=pig_data@data$time), 
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
  
  message(paste0("Running pigs model ", m))

  res <- try(inla (formula, data=inla.stack.data(stack), family="zeroinflatedpoisson1", offset=log(pig_data@data$mem.pigs), control.predictor=list(A=inla.stack.A(stack), link=1), control.inla=list(strategy='adaptive', int.strategy='eb', cmin=0),
               control.compute=list(config = TRUE), num.threads=inla.getOption("num.threads"), control.fixed=list(prec=1, prec.intercept=1)),TRUE)#, openmp.strategy='pardiso.parallel'

  if(class(res)=="try-error"){
  message(paste0("Model ",m," exited with error"))
  }else{
  message(paste0("Completed model ", m) )
   }

  save(file=paste0(folderData, "/pigs", m, ".Rdata"), res)
  
  runtime<-(Sys.time()-t)/60
  message(paste0("Runtime for model ", m, "in minutes : ", runtime))
  return(res)
}
