maptools::gpclibPermit()
getNightlights=function(data, year){
  if(year<=2013){
    nightYear=year
    }else{
  nightYear=2013
  }
  nightRaster = raster(paste('Data/Predictor_data/Nighttime_lights/nightlights_',nightYear,'.tif', sep=""))
  nightlights.crop<-crop(nightRaster, data)
  merged.trans <- spTransform(data, crs(nightlights.crop))
  nightlights.extract<-raster::extract(nightlights.crop, merged.trans)
  nightlights.res<-ifelse(is.na(nightlights.extract)==TRUE, round(mean(nightlights.extract, na.rm=TRUE)), nightlights.extract)
  return(nightlights.res)
}

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
  #crRaster = mask(crRaster, regMoll)
  
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

getPopulation = function(loc, regMap, year){
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
  #crRaster = mask(crRaster, regMoll)
  
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

createSPDE_RE = function(data, idxName, max.edge, offset, prior.range, prior.sigma, regMap = NULL, prevSPDE = NULL, repl = NULL, group = NULL, knots = NULL, grpModel = NULL, constr = FALSE, weightedIntToZero = FALSE, nx = 300, ny = 300){
  # Extract data locations
  loc = as.matrix(data[c("lon", "lat")])
  
  # Extract locations to use for mesh
  if(!is.null(regMap)){
    # Extract boundary from spatial map
    spUnion = unionSpatialPolygons(regMap, rep(1, dim(regMap)[1]))
    loc.mesh = matrix(0, nrow = 0, ncol = 2)
    for(i in 1:length(spUnion@polygons[[1]]@Polygons)){
      loc.mesh = rbind(loc.mesh, spUnion@polygons[[1]]@Polygons[[i]]@coords)
    }
  } else{
    # Extract data locations
    loc.mesh = loc
  }
  
  # Create mesh 
  if(is.null(prevSPDE)){
    mesh = inla.mesh.2d(loc.domain = loc.mesh, 
                        max.edge = max.edge, 
                        offset = offset)
    
    # Add a weighted-by-population-integrate-to-zero
    if(weightedIntToZero){
      # Get extent of region
      ext = extent(regMap)
      xlim = c(ext@xmin, ext@xmax)
      ylim = c(ext@ymin, ext@ymax)
      
      # Make grid on Kenya
      Nx = nx
      Ny = ny
      xs = seq(xlim[1], xlim[2], length.out = Nx)
      ys = seq(ylim[1], ylim[2], length.out = Ny)
      Xs = rep(xs, each = Ny)
      Ys = rep(ys, Nx)
      loc.grid = cbind(Xs, Ys)
      
      # Create mask for removing area outside border
      gridP = data.frame(Longitude = Xs,
                         Latitude = Ys)
      coordinates(gridP) = ~ Longitude + Latitude
      proj4string(gridP) = proj4string(regMap)
      mask1 = as.numeric(over(gridP, regMap)[,1])
      mask1[!is.na(mask1)] = 1
      
      # Get population density using year 2014
      vals = getPopulation(loc.grid, regMap, 2014, kenya)
      vals = vals*mask1
      vals[is.na(vals)] = 0
      weights = vals/sum(vals)
      
      # Make A matrix to approximate integral
      Alocs = inla.spde.make.A(mesh, loc = loc.grid)
      Aweight = inla.spde.make.block.A(Alocs, block = rep(1, Nx*Ny), weights = weights)
      
      # Constraint
      e = matrix(0, nrow = 1, ncol = 1)
      
      # Make constraints
      extraconstr = list(A = Aweight, e = e)
    } else{
      extraconstr = NULL
    }
    
    # Create SPDE object
    spde = inla.spde2.pcmatern(mesh = mesh, 
                               prior.range = prior.range, 
                               prior.sigma = prior.sigma,
                               constr = constr,
                               extraconstr = extraconstr)
  } else{
    mesh = prevSPDE$mesh
    spde = prevSPDE$spde
  }
  
  # Pre-parse effects and ex data objects
  effects = list(eval(parse(text = paste("list(", idxName, " = 1:mesh$n)", sep = ""))))
  exData = eval(parse(text = paste("list(spde.", idxName, " = spde)", sep = "")))
  
  # Construct random effect list
  if(is.null(repl) && is.null(group) && is.null(knots)){
    formula = paste("f(", idxName, ", model = spde.", idxName, ")", sep = "")
  } else{
    tmpForm = paste("f(", idxName, ", model = spde.", idxName, sep = "")
    if(is.null(knots)){
      if(!is.null(repl)){
        tmpForm = paste(tmpForm, ", replicate = ", idxName, ".repl", sep = "")
      }
      if(!is.null(group)){
        tmpForm = paste(tmpForm, ", group = ", idxName, ".group, control.group = list(model = \"", grpModel, "\")", sep = "")
      }
      formula = paste(tmpForm, ")", sep = "")
      A = inla.spde.make.A(mesh, loc = loc, repl = repl, group = group)
      effects = list(inla.spde.make.index(idxName, mesh$n, n.repl = max(c(1, repl), na.rm = TRUE), n.group = max(c(1, group), na.rm = TRUE)))
    } else{
      tmpForm = paste(tmpForm, ", group = ", idxName, ".group, control.group = list(model = \"", grpModel, "\"))", sep = "")
      nGroup = length(knots)
      A = NULL
      idxAll = NULL
      # Iterate through each interval and do linear combination
      knots = c(knots, tail(knots, 1))
      for(pIdx in 1:nGroup){
        if(knots[pIdx] != knots[pIdx+1]){
          idx = which((knots[pIdx] <= data$year) & (data$year < knots[pIdx+1]))
          alpha = (data$year[idx]-knots[pIdx])/(knots[pIdx+1]-knots[pIdx])
        } else{
          idx = which(knots[pIdx] == data$year)
          alpha = 0
        }
        Dalpha = .sparseDiagonal(n = length(idx), x = alpha)
        DI = .sparseDiagonal(n = length(idx), x = 1)
        A1 = inla.spde.make.A(mesh, loc = as.matrix(data[idx, c("lon", "lat")]), group = rep(pIdx, length(idx)), n.group = length(knots)-1)
        pEnd = pIdx+1
        if(pEnd > nGroup)
          pEnd = pIdx
        A2 = inla.spde.make.A(mesh, loc = as.matrix(data[idx, c("lon", "lat")]), group = rep(pEnd, length(idx)), n.group = length(knots)-1)
        if(is.null(A) && (dim(A1)[1] > 0)){
          A = (DI-Dalpha)%*%A1+Dalpha%*%A2
          idxAll = idx
        } else if(dim(A1)[1] > 0){
          A = rbind(A, (DI-Dalpha)%*%A1+Dalpha%*%A2)
          idxAll = c(idxAll, idx)
        }
      }
      knots = knots[1:(length(knots)-1)]
      idxPerm = 1:length(idxAll)
      idxPerm[idxAll] = idxPerm
      A = A[idxPerm,, drop = FALSE]
      
      # Create A matrix
      A=inla.spde.make.A(mesh, loc=loc, n.group=length(knots))
      effects = list(inla.spde.make.index(idxName, mesh$n, n.repl = 1, n.group = length(knots)))
      formula = tmpForm
    }
  }
  rEffect = list(formula = formula,
                 effects = effects,
                 A = A,
                 exData = exData,
                 dontUse = list(mesh = mesh, spde = spde),
                 idx = 1:(mesh$n*max(1, repl, na.rm = TRUE)*max(1, length(knots))))
  
  return(rEffect)
}


fitZIPoisson = function(y, X, rEffects = NULL, fix.prior = NULL, verbose = FALSE, debug = FALSE, withConfig = TRUE, eb = FALSE, offset = NULL){
  
  formula = y ~ -1 + X
  
  # Update formula for each random effect
  if(!is.null(rEffects)){
    for(i in 1:length(rEffects)){
      updForm = paste('~.+', rEffects[[i]]$formula)
      formula = update(formula, as.formula(updForm))
    }
  }
  
  # Initialize A matrices and effects list to go into stack
  A = list(1)
  effects = list(list(X = X))
  
  # Update with each random effect
  if(!is.null(rEffects)){
    for(i in 1:length(rEffects)){
      A = c(A, list(rEffects[[i]]$A))
      effects = c(effects, rEffects[[i]]$effects)
    }
  }
  
  # Update extra data that needs to go into data
  # in INLA call
  exData = list()
  if(!is.null(rEffects)){
    for(i in 1:length(rEffects)){
      exData = c(exData, rEffects[[i]]$exData)
    }
  }
  
  # Create stack
  stk = inla.stack(data = list(y = y),
                   A = A,
                   effects = effects, 
                   tag="est")
  
  # Empirical Bayes?
  if(eb){
    intS = "eb"
  } else{
    intS = "auto"
  }
  
  # Fit model with INLA family=, offset=log(merge2010@data$members), control.predictor=list(A=inla.stack.A(stk.dat10c)))
  #--res.spde10c<-inla (f.spde10, data=inla.stack.data(stk), control.predictor=list(A=inla.stack.A(stk.dat10c)))
  #res.spde10c.poisson<-inla (f.spde10, data=inla.stack.data(stk.dat10c), family="poisson", offset=log(merge2010@data$members), control.predictor=list(A=inla.stack.A(stk.dat10c)))
  
  res = inla(formula = formula,
             data = do.call(what = inla.stack.data, 
                            args = c(list(stack = stk), exData)),
             family = "zeroinflatedpoisson1",
             verbose = verbose,
             control.compute = list(config = withConfig),
             control.predictor = list(A = inla.stack.A(stk), link=1),
             control.fixed = list(mean = fix.prior$mean,
                                  prec = fix.prior$prec),
             control.inla = list(int.strategy = intS), #Warning meessage told me to set h value to this (Eigenvalue 6 of the Hessian is -1748.22 < 0)
             offset = offset, 
             num.threads = inla.getOption("num.threads"))
  
  # Debug
  if(debug){
    print(formula)
  }
  
  # return result
  return(list(res=res, stack=stk))
}

fitModel = function(y, form.cov, data, rEffects = NULL, verbose = FALSE, withConfig = TRUE, eb = FALSE, offset = NULL){
  # Extract number of events among the number at risk
  # Extract covariates
  X = model.frame(form.cov, data)
  #b0 = rep(1, nrow(X))
  #X= cbind(X, b0)
  
  # Fit using binomial model in INLA
  #    this only contains estimates of effects on
  #    hazards and need post-processing
  res = fitZIPoisson(y = y, 
                     X = X,
                     rEffects = rEffects,
                     verbose = verbose,
                     withConfig = withConfig,
                     eb = eb,
                     offset = offset)
  
  # Replace fixed effect names with correct name
  rownames(res$res$summary.fixed) = colnames(X)
  
  return(list(res = res, rEffects = rEffects, covNames = colnames(X)))
}

createPredData = function(year, data.pred, loc.pred, map){
  # Get stratification variables on prediction grid for correct year
  uYear = year
  if(uYear > 2016)
    uYear = 2016
  urbanNum = getUrbanicity(loc = loc.pred, map, uYear)            
  urbanFac = factor(x = urbanNum,
                    levels = c(0, 1))
  
  # Create structure
  data.pred$urban = urbanFac
  data.pred$year = rep(year, nrow(data.pred))
  return(data.pred)
}

projection = function(data.pred,  kIdx= kIdx, form.cov, res, loc.pred, spde){
  #kIdx= knot index
  idx<-1:spde$n.spde + (kIdx -1)*spde$n.spde
  r<-res$res$res$summary.random$i$mean[idx]
  mesh<-spde$mesh
  Aprd<-inla.spde.make.A(mesh, loc.pred)
  r.proj<-Aprd%*%r
  
  # Add covariates to linear predictor
  X.pred = model.frame(form.cov, data.pred)
  X.pred[is.na(X.pred)] = 0
  X.pred =X.pred[,-1]
  coef = res$res$res$summary.fixed[,1]
  
  eta = exp(r.proj + X.pred%*%as.vector(coef))
  
  return(eta)
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
    sIdx = c(sIdx, res.zp$misc$configs$contents$start[iCorr+i])#seems to be where each one starts? Looks like the latent field is just one long vector
    len  = c( len, res.zp$misc$configs$contents$length[iCorr+i])#time REs are each of length 15 (2005-2020), survey is 10, space is 5700, cluster is 4430. Length is 1 for the fixed effects
  }
  
  # Covariate information
  cStart = tail(sIdx, n = 1) + tail(len, n = 1)#sum of last entry in sIx and last entry in length, so covar information must start after latent effects end
  cLen   = dim(X.pred)[2]
  
  # Initialize storage for samples
  wealth.samples = matrix(0, nrow = nSamp, ncol = dim(X.pred)[1])
  
  mesh<-spde$mesh
  Aprd<-inla.spde.make.A(mesh, loc.pred)
  
  # Draw all samples at once
  if(is.null(prevSample)){
    sampAll = inla.posterior.sample(n = nSamp, result = res.zp, intern = TRUE, seed = inla.seed)#draws n samples
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
    sIdx = c(sIdx, res.zp$misc$configs$contents$start[iCorr+i])#seems to be where each one starts? Looks like the latent field is just one long vector
    len  = c( len, res.zp$misc$configs$contents$length[iCorr+i])#time REs are each of length 15 (2005-2020), survey is 10, space is 5700, cluster is 4430. Length is 1 for the fixed effects
  }
  
  # Covariate information
  cStart = tail(sIdx, n = 1) + tail(len, n = 1)#sum of last entry in sIx and last entry in length, so covar information must start after latent effects end
  cLen   = dim(X.pred)[2]
  
  # Initialize storage for samples
  livestock.samples = matrix(0, nrow = nSamp, ncol = dim(X.pred)[1])
  
  mesh<-spde$mesh
  Aprd<-inla.spde.make.A(mesh, loc.pred)
  
  # Draw all samples at once
  if(is.null(prevSample)){
    sampAll = inla.posterior.sample(n = nSamp, result = res.zp, intern = TRUE, seed = inla.seed)#draws n samples
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

aggregateSamples = function(pSamp, loc, regNum, year, regMap, debug = FALSE){
  mask = regNum
  
  # Get indices of each region
  regIdx = list()
  for(i in 1:max(regNum, na.rm = TRUE)){
    regIdx = c(regIdx, list(idx = which(mask == i)))
  }
  
  # Initialize storage for each region
  regVal = list()
  for(i in 1:max(regNum, na.rm = TRUE)){
    regVal =  c(regVal, list(idx = vector('numeric', length = 0)))
  }
  
  gridPop = getPopulation(loc = loc, regMap = regMap, year = year)
  #gridPop = interp.surface(list(x = pop.data$latG[,1], y = pop.data$lonG[1,], z = pop.data$dens), cbind(loc[,2], loc[,1]))
  
  # Generate aggregation for each sample
  nSamp = dim(pSamp)[1]
  for(iSamp in 1:nSamp){
    
    # Use simple integration scheme for each region, one at a time (regIdx= ids for each region)
    for(i in 1:length(regIdx)){
      tmpVal = pSamp[iSamp, regIdx[[i]]]
      tmpDensVal = gridPop[regIdx[[i]]]
      tmpIntVal = sum(tmpVal*tmpDensVal, na.rm = TRUE)/sum(tmpDensVal, na.rm = TRUE)
      regVal[[i]] = c(regVal[[i]], tmpIntVal)
    }    
  }
  
  return(list(regVal = regVal, regIdx = as.vector(mask)))
}