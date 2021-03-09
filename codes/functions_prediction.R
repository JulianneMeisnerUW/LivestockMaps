#-----------------Code adapted from Geir-Arne Fuglstad----------------#

library(raster)
getUrbanicity = function(loc, regMap, year, country){
  
  if(country=="Malawi"){
    s.frames<-c(1998, 2008)
  }
  
  if(country=="Uganda"){
    s.frames<-c(1991, 2002, 2014)
  }
  
  if(country=="DRC"){
    s.frames<-c(2003, 2010)
  }
  
  if(length(which(s.frames==year))>0){
    urbyear=year
  }else{
    urbyear=s.frames[tail(which(s.frames<year),1)]
  }
  
  # Load raster with global urbanicity for year "year"
  urbRaster = raster(paste("Data/Exposure_data/", country, "/Created_datasets/Urbanicity/urbanicity_", urbyear, ".tif", sep=""))
  # Crop raster to region specifced by regMap
  regMoll = spTransform(regMap, urbRaster@crs)
  ext1 = extent(regMoll)
  ext1@xmin = ext1@xmin-0.5
  ext1@xmax = ext1@xmax+0.5
  ext1@ymin = ext1@ymin-0.5
  ext1@ymax = ext1@ymax+0.5
  crRaster = crop(urbRaster, ext1)

  # Transform locations to correct projection
  gridP = data.frame(Longitude = loc[, 1],
                     Latitude =  loc[, 2])
  coordinates(gridP) = ~ Longitude + Latitude
  proj4string(gridP) = proj4string(regMap)
  gridP_moll = spTransform(gridP, urbRaster@crs)

  urban<-extract(crRaster, gridP_moll)

  if(country!="Malawi"){
   urban<-ifelse(urban>=20,1,0)
  }
  
  # Change to required definitien
  #    1: rural
  #    2: urban
  urbMap = urban
  
  return(urbMap)
}

getPopulation = function(loc, regMap, year, country){
  if(country=="Uganda"){
  c.short="uga"
  }
 if(country=="Malawi"){
  c.short="mwi"
  }
 if(country=="DRC"){
  c.short="cod"
  }
  # Load raster with global population for year "year"
  popRaster = raster(paste('Data/Denominator_data/WorldPop/', c.short, '_ppp_', year, '.tif', sep = "")) 
  
  # Crop raster to region specifced by regMap
  regMoll = spTransform(regMap, popRaster@crs)
  ext1 = extent(regMoll)
  ext1@xmin = ext1@xmin-0.5
  ext1@xmax = ext1@xmax+0.5
  ext1@ymin = ext1@ymin-0.5
  ext1@ymax = ext1@ymax+0.5
  crRaster = crop(popRaster, ext1)

  # Transform locations to correct projection
  gridP = data.frame(Longitude = loc[, 1],
                     Latitude =  loc[, 2])
  coordinates(gridP) = ~ Longitude + Latitude
  proj4string(gridP) = proj4string(regMap)
  gridP_moll = spTransform(gridP, popRaster@crs)
  loc.proj = coordinates(gridP_moll)
  
  # Extract grid from raster
  locRaster = coordinates(crRaster)
  xCor = matrix(locRaster[, 1], ncol = dim(crRaster)[1])
  yCor = matrix(locRaster[, 2], ncol = dim(crRaster)[1])
  val = matrix(values(crRaster), ncol = dim(crRaster)[1])
  
  # Nearest neighbour interpolation
  dx = xCor[2,1]-xCor[1,1]
  dy = yCor[1,2]-yCor[1,1]
  xIdx = round((loc.proj[,1]-xCor[1,1])/dx+1)
  xIdx[(xIdx < 1) | (xIdx > dim(crRaster)[2])] = NA
  yIdx = round((loc.proj[,2]-yCor[1,1])/dy+1)
  yIdx[(yIdx < 1) | (yIdx > dim(crRaster)[1])] = NA
  
  # Get values on grid
  val.int = val[dim(val)[1]*(yIdx-1) + xIdx]
  
  return(val.int)
}
samplingGauss = function(res.zp, data.pred, form.cov, year, kIdx=NULL, debug = FALSE, inla.seed = 0L, nSamp, prevSample = NULL, spde, loc.pred, samp.dat){
    
    # Extract covariates
    previous_na_action <- options('na.action')
    options(na.action='na.pass')
    X.pred = model.frame(form.cov, data.pred)
    X.pred<-data.matrix(X.pred)
    options(na.action=previous_na_action$na.action)
    
    # Check if we are using special configs without predictors
    if(!is.null(res.zp$misc$configs$MODIFIED)){
        iCorr = 1
    } else{
        iCorr = 2
    }
    
    # Get information about latent components in fitted model
    sIdx = c()
    len  = c()
    for(i in 1:length(res.zp$model.random)){
        sIdx = c(sIdx, res.zp$misc$configs$contents$start[iCorr+i])
        len  = c( len, res.zp$misc$configs$contents$length[iCorr+i])
    }
    
    # Covariate information
    cStart = tail(sIdx, n = 1) + tail(len, n = 1)
    cLen   = dim(X.pred)[2]
    
    # Initialize storage for samples
    wealth.samples = matrix(0, nrow = nSamp, ncol = dim(X.pred)[1])
    
    mesh<-spde$mesh
    Aprd<-inla.spde.make.A(mesh, loc.pred)
    
    # Draw all samples at once
    if(is.null(prevSample)){
        sampAll = inla.posterior.sample(n = nSamp, result = res.zp, intern = TRUE, seed = inla.seed)
    } else{
        sampAll = prevSample
    }
    
    # Prediction done by sampling
    for(iSamp in 1:nSamp){
        # Draw sample from posterior
        samp = sampAll[[iSamp]]
        
        # Initialize predictor
        eta = matrix(0, nrow = dim(data.pred)[1], ncol = 1)
        
        #Add random effects to predictor
        #1 = space
        #2 = time
        #6 = space-time
        
        #Spatial RE
        spaceIdx<-which(grepl("SPDE", res.zp$model.random)==TRUE)
        timeIdx<-which(grepl("RW", res.zp$model.random)==TRUE)
        res=c(spaceIdx, timeIdx)
        for(i in res){
          if(i==spaceIdx& is.null(kIdx)){
            idx = 1+res.zp$summary.random[[i]]$ID
          }
          if(i==spaceIdx & !is.null(kIdx)){
            start = spde$n.spde*(kIdx-1)+1 
            end = kIdx*spde$n.spde
            id <- c(start:end)
            idx = 1+res.zp$summary.random[[i]]$ID[id]
          }
          if(i==timeIdx){
            id = year - (min(samp.dat$year)-1)
            idx = res.zp$summary.random[[i]]$ID[id]
          }
          
          xRan = samp$latent[sIdx[i]:(sIdx[i]+len[i]-1)] 
          xRan = xRan[idx]
          xRan[is.na(xRan)] = 0
          
          if (i==spaceIdx){
            upd = as.matrix(Aprd%*%xRan)
          }else{
            upd = as.matrix(rep(xRan, nrow(eta)))
          }
          eta = eta + upd
          
        }
        xCov = samp$latent[cStart:(cStart+cLen-1)]
        eta = eta + X.pred%*%as.vector(xCov)
        wealth.samples[iSamp, ] = eta
        }
    return(list(wealth.samples = wealth.samples, sample = sampAll))
}

sampling = function(res.zp, data.pred, form.cov, year, kIdx=NULL, debug = FALSE, inla.seed = 0L, nSamp, prevSample = NULL, spde, loc.pred, samp.dat){
  
  # Extract covariates
  previous_na_action <- options('na.action')
  options(na.action='na.pass')
  X.pred = model.frame(form.cov, data.pred)
  X.pred<-data.matrix(X.pred)
  options(na.action=previous_na_action$na.action)
  
  # Check if we are using special configs without predictors
  if(!is.null(res.zp$misc$configs$MODIFIED)){
    iCorr = 1
  } else{
    iCorr = 2
  }
  
  # Get information about latent components in fitted model
  sIdx = c()
  len  = c()
  for(i in 1:length(res.zp$model.random)){
    sIdx = c(sIdx, res.zp$misc$configs$contents$start[iCorr+i])
    len  = c( len, res.zp$misc$configs$contents$length[iCorr+i])
  }
  
  # Covariate information
  cStart = tail(sIdx, n = 1) + tail(len, n = 1)
  cLen   = dim(X.pred)[2]
  
  # Initialize storage for samples
  livestock.samples = matrix(0, nrow = nSamp, ncol = dim(X.pred)[1])
  
  mesh<-spde$mesh
  Aprd<-inla.spde.make.A(mesh, loc.pred)
  
  # Draw all samples at once
  if(is.null(prevSample)){
    sampAll = inla.posterior.sample(n = nSamp, result = res.zp, intern = TRUE, seed = inla.seed)
  } else{
    sampAll = prevSample
  }
  
  # Prediction done by sampling
  for(iSamp in 1:nSamp){
    # Draw sample from posterior
    samp = sampAll[[iSamp]]
    
    # Initialize predictor
    eta = matrix(0, nrow = dim(data.pred)[1], ncol = 1)
    
    #Add random effects to predictor
    #1 = space
    #2 = time
    #6 = space-time
    
    #Spatial RE
    spaceIdx<-which(grepl("SPDE", res.zp$model.random)==TRUE)
    timeIdx<-which(grepl("RW", res.zp$model.random)==TRUE)
    res=c(spaceIdx, timeIdx)
      for(i in res){
        if(i==spaceIdx& is.null(kIdx)){
          idx = 1+res.zp$summary.random[[i]]$ID
        }
        if(i==spaceIdx & !is.null(kIdx)){
          start = spde$n.spde*(kIdx-1)+1 
          end = kIdx*spde$n.spde
          id <- c(start:end)
          idx = 1+res.zp$summary.random[[i]]$ID[id]
        }
        if(i==timeIdx){
          id = year - (min(samp.dat$year)-1)
          idx = res.zp$summary.random[[i]]$ID[id]
        }
          
          xRan = samp$latent[sIdx[i]:(sIdx[i]+len[i]-1)] 
          xRan = xRan[idx]
          xRan[is.na(xRan)] = 0
          
          if (i==spaceIdx){
            upd = as.matrix(Aprd%*%xRan)
          }else{
            upd = as.matrix(rep(xRan, nrow(eta)))
          }
          eta = eta + upd
      }
      xCov = samp$latent[cStart:(cStart+cLen-1)]
      eta = exp(eta + X.pred%*%as.vector(xCov))
      livestock.samples[iSamp, ] = eta
    }
  return(list(livestock.samples = livestock.samples, sample = sampAll))
}

aggregateSamples = function(pSamp, loc, crs, regNum, year, regMap, debug = FALSE, level=NULL, country){
  #Get region indices
  loc2<-SpatialPoints(loc, proj4string=crs)
  merge<-over(loc2, regMap)
  if(level=="region"){
   if(country=="Malawi"){
   regIdx = as.numeric(merge$ADM1_EN)
   }
   if(country=="Uganda"){
   regIdx = as.numeric(merge$NAME_4)
   }
  }
  if(level=="district"){
   if(country=="Malawi"){
   regIdx = as.numeric(merge$ADM2_EN)
   }
   if(country=="Uganda"){
    regIdx = as.numeric(merge$NAME_3)
   }
  }
  if(level=="national"){
   if(country=="Malawi"){
   regIdx = as.numeric(merge$ADM0_EN)
   }
   if(country=="Uganda"){
   regIdx = as.numeric(merge$NAME_0)
   }
  }

  mask = regNum
  
  # Initialize storage for each region
  regVal = list()
  for(i in 1:max(regNum, na.rm = TRUE)){
    regVal =  c(regVal, list(idx = vector('numeric', length = 0)))
  }
  
  gridPop = getPopulation(loc = loc, regMap = regMap, year = year, country=country)
  #gridPop = interp.surface(list(x = pop.data$latG[,1], y = pop.data$lonG[1,], z = pop.data$dens), cbind(loc[,2], loc[,1]))
  
  # Generate aggregation for each sample
  nSamp = dim(pSamp)[1]
  for(iSamp in 1:nSamp){
    
    # Use simple integration scheme for each region, one at a time (regIdx= ids for each region)
    for(i in regNum){
      idxes<-which(regIdx==i)
      tmpVal = pSamp[iSamp, idxes]
      tmpDensVal = gridPop[idxes]
      tmpIntVal = sum(tmpVal*tmpDensVal, na.rm = TRUE)/sum(tmpDensVal, na.rm = TRUE)
      regVal[[i]] = c(regVal[[i]], tmpIntVal)
    }    
  }
  
  return(list(regVal = regVal, regIdx = regIdx))
}

