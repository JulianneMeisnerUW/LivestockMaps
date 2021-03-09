  #--------------------#
  # Load packages
  #--------------------#
  library(rgdal)
  library(INLA)
  library(SpatialEpi)
  library(rmeta)
  library(classInt)
  library(RColorBrewer)
  inla.setOption(mkl=TRUE)
  inla.setOption("num.threads", 8)
  
  #--------------------#
  # Load files
  #--------------------#
  set.seed(1154)
  folderData<-"Data/Exposure_data/South_Sudan/Created_datasets/"
  census<-readOGR(paste(folderData, "census_data", sep=""), layer="census_data")
  
  pc.u=1 
  pc.alpha = 0.01
  phi.u = 0.5
  phi.alpha = 0.7
  
  bym_hyper = list(
    phi = list( prior="pc", param=c(phi.u, phi.alpha)), 
    prec = list( prior="pc.prec", param = c(pc.u, pc.alpha))
  )
  
  hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u, 
                                                            pc.alpha)))
  
  f_c<-ld_c~ 1 +f(struct, graph="census.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="census.graph", model="iid", param=hyperpc1)
  
  f_p<-ld_p~ 1 +f(struct, graph="census.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="census.graph", model="iid", param=hyperpc1)
  
  CI=0.95 
  
  cattle.smooth <- inla(f_c, family = "gaussian", data = census@data, 
                     control.predictor = list(compute = TRUE), 
                     control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                     scale = census@data$prec_c, 
                     lincomb = NULL, 
                     quantiles = c((1 - CI)/2, 0.5, 1 - (1 - CI)/2))
  
  pigs.smooth <- inla(f_p, family = "gaussian", data = census@data, 
                        control.predictor = list(compute = TRUE), 
                        control.family = list( hyper = list(prec = list( initial = log(1), fixed=TRUE))), 
                        scale = census@data$prec_p, 
                        lincomb = NULL, 
                        quantiles = c((1 - CI)/2, 0.5, 1 - (1 - CI)/2))
  
  
  #Results
  cattle.smooth$summary.fixed[,c('mean', '0.5quant', 'sd')]
  cattle.smooth$summary.hyperpar[,c('mean', '0.5quant')]
  cattle.smooth$summary.hyperpar[,c('sd')]
  
  fixed.med.c <- rep(cattle.smooth$summary.fixed[,4],dim(census@data)[1])
  random.iid.c <- cattle.smooth$summary.random$unstruct[,5]
  random.smooth.c <- cattle.smooth$summary.random$struct[,5]
  linpred.c <- cattle.smooth$summary.fitted.values[,'0.5quant']
  pred.c <- exp(linpred.c)
  sd.c <- exp(cattle.smooth$summary.fitted.values[,'sd'])
  res.c <- cbind(census@data,fixed.med.c,random.iid.c,random.smooth.c,linpred.c,pred.c)
  census@data$pred.c<-pred.c
  census@data$pred.c.sd<-sd.c
  
  med.palette <- brewer.pal(n = 7, name = "Purples")
  med.int <- classIntervals(round(res.c$pred.c, 3),
                            n = 7, style = 'jenks')
  med.col <- findColours(med.int, med.palette)
  
  pigs.smooth$summary.fixed[,c('mean', '0.5quant', 'sd')]
  pigs.smooth$summary.hyperpar[,c('mean', '0.5quant')]
  pigs.smooth$summary.hyperpar[,c('sd')]
  
  fixed.med.p <- rep(pigs.smooth$summary.fixed[,4],dim(census@data)[1])
  random.iid.p <- pigs.smooth$summary.random$unstruct[,5]
  random.smooth.p <- pigs.smooth$summary.random$struct[,5]
  linpred.p <- pigs.smooth$summary.fitted.values[,'0.5quant']
  pred.p <- exp(linpred.p)
  sd.p <- exp(pigs.smooth$summary.fitted.values[,'sd'])
  res.p <- cbind(census@data,fixed.med.p,random.iid.p,random.smooth.p,linpred.p,pred.p)
  res.p$log_pred<-log(res.p$pred.p)
  census@data$pred.p<-pred.p
  census@data$pred.p.sd<-sd.p
  
  writeOGR(obj=census, dsn=paste(folderData, "census_data", sep=""), layer="SS_HAT", driver="ESRI Shapefile", overwrite_layer = TRUE)
  
  distribute<-census[,c("new_adm2", "Shape_area", "ADM1_EN", "ADM1_PCODE", "pred_c", "pred_c_sd", "pred_p", "pred_p_sd")]
  
  distribute@data$pigs_median2008<-distribute@data$pred_p
  distribute@data$pigs_sd2008<-distribute@data$pred_p_sd
  
  distribute@data$cattle_median2008<-distribute@data$pred_c
  distribute@data$cattle_sd2008<-distribute@data$pred_c_sd
  
  distribute<-distribute[,c("new_adm2", "Shape_area", "ADM1_EN", "ADM1_PCODE", "cattle_median2008", "cattle_sd2008", "pigs_median2008", "pigs_sd2008")]
  
  writeOGR(obj=distribute, dsn=paste(folderData, "census_data", sep=""), layer="SS_HAT_distribute", driver="ESRI Shapefile", overwrite_layer = TRUE)
  