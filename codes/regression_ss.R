library(rgdal)
library(spdep)
library(raster)
library(dplyr)
library(INLA)
library(xtable)
library(classInt)
library(RColorBrewer)
library(uwIntroStats)
library(robustHD)
outcomeFolderData<-"Data/Outcome Data/"

source("codes/functions_prediction.R")

folderData<-paste0("results/regression/South_Sudan/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

HAT<-readOGR("Outcome_data/South_Sudan/HAT_final_SS", layer="HAT_final_SS")
HAT.active<-readOGR("Outcome_data/South_Sudan/HAT_final_SS", layer="HAT_final_SS_active")
HAT.active<-HAT.active[which(HAT.active@data$n2008>0),]

nb.map <- poly2nb(HAT)
nb2INLA("HAT.graph",nb.map)

nb.map.active <- poly2nb(HAT.active)
nb2INLA("HAT.graph.active",nb.map.active)

#-------------------#
# Main model
#-------------------#
#Put exposure on log scale, and transform so that a 1 unit increase = 50% increase
logtrans<-log(1.5)
density_log<-log(HAT@data$pred_c)
HAT@data$c_trans<-density_log/logtrans
density_log<-log(HAT@data$pred_p)
HAT@data$p_trans<-density_log/logtrans

density_log<-log(HAT.active@data$pred_c)
HAT.active@data$c_trans<-density_log/logtrans
density_log<-log(HAT.active@data$pred_p)
HAT.active@data$p_trans<-density_log/logtrans

pc.u=1 
pc.alpha = 0.01

hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u, 
                                                          pc.alpha)))

HAT@data$struct<-HAT@data$unstruct<-1:nrow(HAT@data)
HAT.active@data$struct<-HAT.active@data$unstruct<-1:nrow(HAT.active@data)

#-----------------------------------#
#Naive fit ignoring measurement error
#-----------------------------------#

formula_c<-cs2008~ c_trans + elevation + wealth + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cWP <- inla(formula_c, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                      control.predictor = list(compute = TRUE))

results_cLS <- inla(formula_c, family = "poisson", data = HAT@data, offset=log(HAT@data$LS2008),
                    control.predictor = list(compute = TRUE))

formula_p<-cs2008~ p_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pWP <- inla(formula_p, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                  control.predictor = list(compute = TRUE))

h.value=2.59364
results_pLS <- inla(formula_p, family = "poisson", data = HAT@data, offset=log(HAT@data$LS2008),
                    control.predictor = list(compute = TRUE), control.inla=list(h=h.value))

#-----------------------------------#
#Priors for MEC model
#-----------------------------------#
n=nrow(HAT@data)
d<-runif(n, min=0.5, max=1.5)
lgamma.prior=c(10,10) 

init.beta_cWP<-results_cWP$summary.fixed[2,1]
beta_cWP.var<-(results_cWP$summary.fixed[2,2])^2
init.beta_cLS<-results_cLS$summary.fixed[2,1]
beta_cLS.var<-(results_cLS$summary.fixed[2,2])^2
init.beta_pWP<-results_pWP$summary.fixed[2,1]
beta_pWP.var<-(results_pWP$summary.fixed[2,2])^2
init.beta_pLS<-results_pLS$summary.fixed[2,1]
beta_pLS.var<-(results_pLS$summary.fixed[2,2])^2

#Wealth
init.beta_wcWP<-results_cWP$summary.fixed[4,1]
beta_wcWP.var<-(results_cWP$summary.fixed[4,2])^2
init.beta_wcLS<-results_cLS$summary.fixed[4,1]
beta_wcLS.var<-(results_cLS$summary.fixed[4,2])^2
init.beta_wpWP<-results_pWP$summary.fixed[4,1]
beta_wpWP.var<-(results_pWP$summary.fixed[4,2])^2
init.beta_wpLS<-results_pLS$summary.fixed[4,1]
beta_wpLS.var<-(results_pLS$summary.fixed[4,2])^2

beta.prior_cWP=c(init.beta_cWP, 1/beta_cWP.var)
beta.prior_cLS=c(init.beta_cLS, 1/beta_cLS.var)
beta.prior_pWP=c(init.beta_pWP, 1/beta_pWP.var)
beta.prior_pLS=c(init.beta_pLS, 1/beta_pLS.var)

#Wealth
beta.prior_wcWP=c(init.beta_wcWP, 1/beta_wcWP.var)
beta.prior_wcLS=c(init.beta_wcLS, 1/beta_wcLS.var)
beta.prior_wpWP=c(init.beta_wpWP, 1/beta_wpWP.var)
beta.prior_wpLS=c(init.beta_wpLS, 1/beta_wpLS.var)

#Loggamma
cvar<-HAT@data$pred_c_sd^2
pvar<-HAT@data$pred_p_sd^2
wvar<-HAT@data$wealth_se^2
init.prec.u.c<-1/var(cvar)
init.prec.u.p<-1/var(pvar)
init.prec.u.w<-1/var(wvar)

#Do this as fixed, so no prior distribution
init.mean.c<-mean(HAT@data$c_trans)
init.mean.p<-mean(HAT@data$p_trans)
init.mean.w<-mean(HAT@data$wealth)

#Loggamma
init.prec.x.c<-1/var(HAT@data$c_trans)
init.prec.x.p<-1/var(HAT@data$p_trans)
init.prec.x.w<-1/var(HAT@data$wealth)

formula_cWP<-cs2008~ f(c_trans, model="mec", scale=d, values=c_trans,
                     hyper=list(
                       beta=list(prior="gaussian", param=beta.prior_cWP, fixed=FALSE),
                       prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                       mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                       prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + elevation + wealth +  war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cWP_mec <- inla(formula_cWP, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                    control.predictor = list(compute = TRUE))

exp(results_cWP_mec$summary.hyperpar)

formula_cLS<-cs2008~ f(c_trans, model="mec", scale=d, values=c_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_cLS, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + elevation + wealth + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cLS_mec <- inla(formula_cLS, family = "poisson", data = HAT@data, offset=log(HAT@data$LS2008),
                        control.predictor = list(compute = TRUE))

exp(results_cLS_mec$summary.hyperpar)

formula_pWP<-cs2008~ f(p_trans, model="mec", scale=d, values=p_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_pWP, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + elevation + wealth + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pWP_mec <- inla(formula_pWP, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                        control.predictor = list(compute = TRUE))

exp(results_pWP_mec$summary.hyperpar)

formula_pLS<-cs2008~ f(p_trans, model="mec", scale=d, values=p_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_pLS, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

h=1.16118
results_pLS_mec <- inla(formula_pLS, family = "poisson", data = HAT@data, offset=log(HAT@data$LS2008),
                        control.predictor = list(compute = TRUE), control.inla=list(h=h.value))

exp(results_pLS_mec$summary.hyperpar)

#---------------------#
# With MEC for wealth
#---------------------#
formula_cWP<-cs2008~ f(c_trans, model="mec", scale=d, values=c_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_cWP, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                          hyper=list(
                                                                                                                            beta=list(prior="gaussian", param=beta.prior_wcWP, fixed=FALSE),
                                                                                                                            prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                            mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                            prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

h.value=0.010814
results_cWP_mec2 <- inla(formula_cWP, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                        control.predictor = list(compute = TRUE), control.inla=list(h=h.value))

exp(results_cWP_mec2$summary.hyperpar)

formula_cLS<-cs2008~ f(c_trans, model="mec", scale=d, values=c_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_cLS, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                          hyper=list(
                                                                                                                            beta=list(prior="gaussian", param=beta.prior_wcLS, fixed=FALSE),
                                                                                                                            prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                            mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                            prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

h.value=  0.0164338
results_cLS_mec2 <- inla(formula_cLS, family = "poisson", data = HAT@data, offset=log(HAT@data$LS2008),
                        control.predictor = list(compute = TRUE), control.inla=list(h=h.value))

exp(results_cLS_mec2$summary.hyperpar)

formula_pWP<-cs2008~ f(p_trans, model="mec", scale=d, values=p_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_pWP, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                          hyper=list(
                                                                                                                            beta=list(prior="gaussian", param=beta.prior_wpWP, fixed=FALSE),
                                                                                                                            prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                            mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                            prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

h.value=0.0111998
results_pWP_mec2 <- inla(formula_pWP, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                        control.predictor = list(compute = TRUE), control.inla=list(h=h.value))

exp(results_pWP_mec2$summary.hyperpar)

formula_pLS<-cs2008~ f(p_trans, model="mec", scale=d, values=p_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_pLS, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                          hyper=list(
                                                                                                                            beta=list(prior="gaussian", param=beta.prior_wpLS, fixed=FALSE),
                                                                                                                            prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                            mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                            prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pLS_mec2 <- inla(formula_pLS, family = "poisson", data = HAT@data, offset=log(HAT@data$LS2008),
                        control.predictor = list(compute = TRUE))

exp(results_pLS_mec2$summary.hyperpar)

#-----------------------------------#
#Mediation analysis, naive fit
#-----------------------------------#
#Outcome model
formula_cLST<-cs2012~ c_trans + elevation  + wealth  + war2007 + LST_2010 + LST_2007 + NDVI2007 + war2009+ f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cLST <- inla(formula_cLST, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                    control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_cNDVI<-cs2012~ c_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + NDVI2010 + war2009 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cNDVI <- inla(formula_cNDVI, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                     control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_pLST<-cs2012~ p_trans + elevation + wealth  + war2007 + LST_2010 + LST_2007 + NDVI2007 + war2009+ f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pLST <- inla(formula_pLST, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                     control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_pNDVI<-cs2012~ p_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + NDVI2010 + war2009 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pNDVI <- inla(formula_pNDVI, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                      control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

#Mediator model
formula_cLST2<-LST_2010~ c_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cLST2 <- inla(formula_cLST2, family = "gaussian", data = HAT@data, control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_cNDVI2<-NDVI2010~ c_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cNDVI2 <- inla(formula_cNDVI2, family = "gaussian", data = HAT@data,
                      control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_pLST2<-LST_2010~ p_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pLST2 <- inla(formula_pLST2, family = "gaussian", data = HAT@data,
                     control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_pNDVI2<-NDVI2010~ p_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pNDVI2 <- inla(formula_pNDVI2, family = "gaussian", data = HAT@data,
                      control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))
#---------------------------------------------------#
# Mediation analysis with MEC for density and wealth
#---------------------------------------------------#
init.beta_cLST<-results_cLST$summary.fixed[2,1]
beta_cLST.var<-(results_cLST$summary.fixed[2,2])^2
init.beta_cNDVI<-results_cNDVI$summary.fixed[2,1]
beta_cNDVI.var<-(results_cNDVI$summary.fixed[2,2])^2
init.beta_pLST<-results_pLST$summary.fixed[2,1]
beta_pLST.var<-(results_pLST$summary.fixed[2,2])^2
init.beta_pNDVI<-results_pNDVI$summary.fixed[2,1]
beta_pNDVI.var<-(results_pNDVI$summary.fixed[2,2])^2

init.beta_cLST2<-results_cLST2$summary.fixed[2,1]
beta_cLST.var2<-(results_cLST2$summary.fixed[2,2])^2
init.beta_cNDVI2<-results_cNDVI2$summary.fixed[2,1]
beta_cNDVI.var2<-(results_cNDVI2$summary.fixed[2,2])^2
init.beta_pLST2<-results_pLST2$summary.fixed[2,1]
beta_pLST.var2<-(results_pLST2$summary.fixed[2,2])^2
init.beta_pNDVI2<-results_pNDVI2$summary.fixed[2,1]
beta_pNDVI.var2<-(results_pNDVI2$summary.fixed[2,2])^2

#Wealth
init.beta_wcLST<-results_cLST$summary.fixed[4,1]
beta_wcLST.var<-(results_cLST$summary.fixed[4,2])^2
init.beta_wcNDVI<-results_cNDVI$summary.fixed[4,1]
beta_wcNDVI.var<-(results_cNDVI$summary.fixed[4,2])^2
init.beta_wpLST<-results_pLST$summary.fixed[4,1]
beta_wpLST.var<-(results_pLST$summary.fixed[4,2])^2
init.beta_wpNDVI<-results_pNDVI$summary.fixed[4,1]
beta_wpNDVI.var<-(results_pNDVI$summary.fixed[4,2])^2

init.beta_wcLST2<-results_cLST2$summary.fixed[4,1]
beta_wcLST.var2<-(results_cLST2$summary.fixed[4,2])^2
init.beta_wcNDVI2<-results_cNDVI2$summary.fixed[4,1]
beta_wcNDVI.var2<-(results_cNDVI2$summary.fixed[4,2])^2
init.beta_wpLST2<-results_pLST2$summary.fixed[4,1]
beta_wpLST.var2<-(results_pLST2$summary.fixed[4,2])^2
init.beta_wpNDVI2<-results_pNDVI2$summary.fixed[4,1]
beta_wpNDVI.var2<-(results_pNDVI2$summary.fixed[4,2])^2

beta.prior_cLST=c(init.beta_cLST, 1/beta_cLST.var)
beta.prior_cNDVI=c(init.beta_cNDVI, 1/beta_cNDVI.var)
beta.prior_pLST=c(init.beta_pLST, 1/beta_pLST.var)
beta.prior_pLS=c(init.beta_pLST, 1/beta_pLST.var)

beta.prior_cLST2=c(init.beta_cLST2, 1/beta_cLST.var)
beta.prior_cNDVI2=c(init.beta_cNDVI2, 1/beta_cNDVI.var)
beta.prior_pLST2=c(init.beta_pLST2, 1/beta_pLST.var)
beta.prior_pLS2=c(init.beta_pLST2, 1/beta_pLST.var)

#Wealth
beta.prior_wcLST=c(init.beta_wcLST, 1/beta_wcLST.var)
beta.prior_wcNDVI=c(init.beta_wcNDVI, 1/beta_wcNDVI.var)
beta.prior_wpLST=c(init.beta_wpLST, 1/beta_wpLST.var)
beta.prior_wpNDVI=c(init.beta_wpNDVI, 1/beta_wpNDVI.var)

beta.prior_wcLST2=c(init.beta_wcLST2, 1/beta_wcLST.var)
beta.prior_wcNDVI2=c(init.beta_wcNDVI2, 1/beta_wcNDVI.var)
beta.prior_wpLST2=c(init.beta_wpLST2, 1/beta_wpLST.var)
beta.prior_wpNDVI2=c(init.beta_wpNDVI2, 1/beta_wpNDVI.var)

formula_cLST<-cs2012~ f(c_trans, model="mec", scale=d, values=c_trans,
                       hyper=list(
                         beta=list(prior="gaussian", param=beta.prior_cLST, fixed=FALSE),
                         prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                         mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                         prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                          hyper=list(
                                                                                                                            beta=list(prior="gaussian", param=beta.prior_wcLST, fixed=FALSE),
                                                                                                                            prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                            mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                            prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + LST_2010 + war2009 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)
results_cLST_mec <- inla(formula_cLST, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                         control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

exp(results_cLST_mec$summary.hyperpar)

formula_cNDVI<-cs2012~ f(c_trans, model="mec", scale=d, values=c_trans,
                        hyper=list(
                          beta=list(prior="gaussian", param=beta.prior_cLST, fixed=FALSE),
                          prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                          mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                          prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                           hyper=list(
                                                                                                                             beta=list(prior="gaussian", param=beta.prior_wcLST, fixed=FALSE),
                                                                                                                             prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                             mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                             prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + NDVI2010 + war2009 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cNDVI_mec <- inla(formula_cNDVI, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                         control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

exp(results_cNDVI_mec$summary.hyperpar)

formula_pLST<-cs2012~ f(p_trans, model="mec", scale=d, values=p_trans,
                        hyper=list(
                          beta=list(prior="gaussian", param=beta.prior_pLST, fixed=FALSE),
                          prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                          mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                          prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                           hyper=list(
                                                                                                                             beta=list(prior="gaussian", param=beta.prior_wpLST, fixed=FALSE),
                                                                                                                             prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                             mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                             prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + LST_2010 + war2009 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)
h.value=0.0116792
results_pLST_mec <- inla(formula_pLST, family = "poisson", data = HAT@data, offset=log(HAT@data$WP2008),
                         control.predictor = list(compute = TRUE), control.inla=list(h=h.value),control.compute=list(config = TRUE))

exp(results_pLST_mec$summary.hyperpar)

formula_pNDVI<-cs2012~ f(p_trans, model="mec", scale=d, values=p_trans,
                         hyper=list(
                           beta=list(prior="gaussian", param=beta.prior_pLST2, fixed=FALSE),
                           prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                           mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                           prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                            hyper=list(
                                                                                                                              beta=list(prior="gaussian", param=beta.prior_wpLST2, fixed=FALSE),
                                                                                                                              prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                              mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                              prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + NDVI2010 + war2009 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pNDVI_mec <- inla(formula_pNDVI, family = "poisson", data=HAT@data, offset=log(HAT@data$WP2008),
                          control.predictor = list(compute = TRUE), control.inla=list(h=h.value),control.compute=list(config = TRUE))

formula_cLST2<-LST_2010~ f(c_trans, model="mec", scale=d, values=c_trans,
                        hyper=list(
                          beta=list(prior="gaussian", param=beta.prior_cLST2, fixed=FALSE),
                          prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                          mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                          prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                           hyper=list(
                                                                                                                             beta=list(prior="gaussian", param=beta.prior_wcLST2, fixed=FALSE),
                                                                                                                             prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                             mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                             prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)
h.value=0.0419701
results_cLST_mec2 <- inla(formula_cLST2, family = "gaussian", data = HAT@data,
                         control.predictor = list(compute = TRUE), control.inla=list(h=h.value),control.compute=list(config = TRUE))


formula_cNDVI2<-NDVI2010~ f(c_trans, model="mec", scale=d, values=c_trans,
                         hyper=list(
                           beta=list(prior="gaussian", param=beta.prior_cLST2, fixed=FALSE),
                           prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.c), fixed=FALSE),
                           mean.x=list(prior="gaussian", initial=init.mean.c, fixed=TRUE),
                           prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.c), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                            hyper=list(
                                                                                                                              beta=list(prior="gaussian", param=beta.prior_wcLST2, fixed=FALSE),
                                                                                                                              prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                              mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                              prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 +f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_cNDVI_mec2 <- inla(formula_cNDVI2, family = "gaussian", data = HAT@data, 
                          control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_pLST2<-LST_2010~ f(p_trans, model="mec", scale=d, values=p_trans,
                        hyper=list(
                          beta=list(prior="gaussian", param=beta.prior_pLST2, fixed=FALSE),
                          prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                          mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                          prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                           hyper=list(
                                                                                                                             beta=list(prior="gaussian", param=beta.prior_wpLST2, fixed=FALSE),
                                                                                                                             prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                             mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                             prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation  + war2007 + LST_2007 + NDVI2007 +  f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)
results_pLST_mec2 <- inla(formula_pLST2, family = "gaussian", data = HAT@data,
                         control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

formula_pNDVI2<-NDVI2010~ f(p_trans, model="mec", scale=d, values=p_trans,
                         hyper=list(
                           beta=list(prior="gaussian", param=beta.prior_pLST2, fixed=FALSE),
                           prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.p), fixed=FALSE),
                           mean.x=list(prior="gaussian", initial=init.mean.p, fixed=TRUE),
                           prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.p), fixed=FALSE))) + f(wealth, model="mec", scale=d, values=wealth,
                                                                                                                            hyper=list(
                                                                                                                              beta=list(prior="gaussian", param=beta.prior_wpLST2, fixed=FALSE),
                                                                                                                              prec.u=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.u.w), fixed=FALSE),
                                                                                                                              mean.x=list(prior="gaussian", initial=init.mean.w, fixed=TRUE),
                                                                                                                              prec.x=list(prior="loggamma", param=lgamma.prior, initial=log(init.prec.x.w), fixed=FALSE))) + elevation + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph", model="iid", param=hyperpc1)

results_pNDVI_mec2 <- inla(formula_pNDVI2, family = "gaussian", data = HAT@data,
                          control.predictor = list(compute = TRUE),control.compute=list(config = TRUE))

#-------------------------#
# Active surveillance only#
#-------------------------#
formula_c<-cs2008~ c_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph.active", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph.active", model="iid", param=hyperpc1)

results_c <- inla(formula_c, family = "poisson", data = HAT.active@data, offset=log(HAT.active@data$n2008),
                    control.predictor = list(compute = TRUE))

formula_p<-cs2008~ p_trans + elevation + wealth  + war2007 + LST_2007 + NDVI2007 + f(struct, graph="HAT.graph.active", model="besag", param=hyperpc1, scale.model=TRUE)+f(unstruct, graph="HAT.graph.active", model="iid", param=hyperpc1)

results_p <- inla(formula_p, family = "poisson", data = HAT.active@data, offset=log(HAT.active@data$n2008),
                    control.predictor = list(compute = TRUE))

results=list("WP_main_cattle"=results_cWP, "LS_main_cattle"=results_cLS,
             "WP_MEC2_cattle"=results_cWP_mec2, "LS_MEC2_cattle"=results_cLS_mec2,
             "active_cattle"=results_c, "WP_main_pigs"=results_pWP, 
             "LS_main_pigs"=results_pLS, "WP_MEC2_pigs"=results_pLS_mec2, 
             "LS_MEC2_pigs"=results_pLS_mec2, "active_pigs"=results_p)

save(results, file="results/regression/South_Sudan/SS_regression_results.Rdata")

med_results=list("cattle_naive_NDVI"= results_cNDVI, "cattle_naive_LST"=results_cLST,
                 "pigs_naive_NDVI" = results_pNDVI, "pigs_naive_LST"=results_pLST,
                 "cattle_naive_NDVI2"= results_cNDVI2, "cattle_naive_LST2"=results_cLST2,
                 "pigs_naive_NDVI2" = results_pNDVI2, "pigs_naive_LST2"=results_pLST2,
                 "cattle_NDVI_MEC"= results_cNDVI_mec, "cattle_LST_MEC"=results_cLST_mec,
                 "pigs_NDVI_MEC" = results_pNDVI_mec, "pigs_LST_MEC"=results_pLST_mec,
                 "cattle_NDVI_MEC2"= results_cNDVI_mec2, "cattle_LST_MEC2"=results_cLST_mec2,
                 "pigs_NDVI_MEC2" = results_pNDVI_mec2, "pigs_LST2_MEC2"=results_pLST_mec2)

save(med_results, file="results/regression/South_Sudan/SS_mediation_results.Rdata")

results_table<-as.data.frame(matrix(nrow=length(results), ncol=2))
colnames(results_table)<-c("Model", "Posterior Mean (95% CI)")
results_table$Model<-names(results)

results_table_plot<-as.data.frame(matrix(nrow=length(results), ncol=4))
colnames(results_table_plot)<-c("Model", "Mean", "Lower", "Upper")
results_table_plot$Model<-names(results)

for(i in 1:length(results)){
  res<-results[[i]]
  if(rownames(results[[i]]$summary.fixed)[2]=="elevation"){
    print(paste(rownames(res$summary.hyperpar)[1]),i)
    mean<-round(exp(res$summary.hyperpar)[1,1],2)
    lower<-round(exp(res$summary.hyperpar)[1,3],2)
    upper<-round(exp(res$summary.hyperpar)[1,5],2)
    results_table[i,2]<-paste(mean, " (", lower, ", ", upper, ")", sep="")
    
    results_table_plot$Mean[i]<-mean
    results_table_plot$Lower[i]<-lower
    results_table_plot$Upper[i]<-upper
  }else{
      print(paste(rownames(res$summary.fixed)[2]),i)
      mean<-round(exp(res$summary.fixed)[2,1],2)
      lower<-round(exp(res$summary.fixed)[2,3],2)
      upper<-round(exp(res$summary.fixed)[2,5],2)
      results_table[i,2]<-paste(mean, " (", lower, ", ", upper, ")", sep="")
      
      results_table_plot$Mean[i]<-mean
      results_table_plot$Lower[i]<-lower
      results_table_plot$Upper[i]<-upper
  }
}

restab<-xtable(results_table)

save(restab, file="results/regression/South_Sudan/SS_results_table.Rdata")

#-----------------------#
# Mediation analysis
#-----------------------#
#Draw samples from posterior
outcomeIdx<-c(1:4,9:12)
med_samps<-as.data.frame(matrix(NA, nrow=1000, ncol=24))
spp<-c("cattle", "pigs")
meds<-c("NDVI", "LST")
coefs<-c("theta1", "theta2", "beta1")
mec<-c("", "_MEC")
c1<-paste(rep(spp, each = length(meds)), meds, sep = "_")
c2<-paste(rep(c1, each=length(coefs)), coefs, sep="_")
ns<-paste0(rep(c2, each=length(mec)), mec)
colnames(med_samps)<-ns

for(i in 1:length(med_results)){
  sampAll = inla.posterior.sample(n = 1000, result = med_results[[i]], intern = TRUE, seed = 0L)
  if(grepl("NDVI", names(med_results)[i])==TRUE){
    mediator<-paste0("NDVI", 2010)
    colIdx_1<-which(grepl("NDVI", names(med_samps)))
  }else{
    mediator<-paste("LST", 2010, sep="_")
    colIdx_1<-which(grepl("LST", names(med_samps)))
  }
  if(grepl("pig", names(med_results)[i])){
    colIdx_2<-which(grepl("pigs", names(med_samps)))
  }else{
    colIdx_2<-which(grepl("cattle", names(med_samps)))
  }
  
  if(grepl("MEC", names(med_results)[i])==TRUE){
    MEC<-TRUE
    colIdx_3<-which(grepl("MEC", names(med_samps)))
  }else{
    MEC<-FALSE
    colIdx_3<-which(grepl("MEC", names(med_samps))==FALSE)
  }
  
  cIdx<-colIdx_2[which(colIdx_2%in%colIdx_1)]
  
  cIdx<-cIdx[which(cIdx%in%colIdx_3)]
  
  if(length(cIdx)!=3){
    warning("Problem with assigning output column")
  }
  if(i %in% outcomeIdx){
    collapse<-sum(med_samps[1, cIdx[1]],med_samps[1, cIdx[2]], na.rm=T)
  }else{
    collapse<-sum(med_samps[1, cIdx[3]], na.rm=T)
  }
  if(collapse!=0){
    warning(paste0("Something is already here! Index ", i))
  }
  
  if(i %in% outcomeIdx){
    if(MEC==TRUE){
      theta1_idx<-which(grepl("_trans", names(sampAll[[1]]$hyperpar)))[1]
    }else{
      theta1_idx<-which(grepl("_trans", rownames(sampAll[[1]]$latent)))
    }
    theta2_idx<-which(grepl(mediator, rownames(sampAll[[1]]$latent)))
    
    for(j in 1:length(sampAll)){
      if(MEC==TRUE){
        med_samps[j,cIdx[1]]<-as.numeric(sampAll[[j]]$hyperpar[theta1_idx])
        med_samps[j,cIdx[2]]<-as.numeric(sampAll[[j]]$latent[theta2_idx])
      }else{
        med_samps[j,cIdx[1]]<-as.numeric(sampAll[[j]]$latent[theta1_idx])
        med_samps[j,cIdx[2]]<-as.numeric(sampAll[[j]]$latent[theta2_idx])
      }
    }
  }else{
    if(MEC==TRUE){
      beta1_idx<-which(grepl("_trans", names(sampAll[[1]]$hyperpar)))[1]
    }else{
      beta1_idx<-which(grepl("_trans", rownames(sampAll[[1]]$latent)))
    }
    for(j in 1:length(sampAll)){
      if(MEC==TRUE){
        med_samps[j,cIdx[3]]<-as.numeric(sampAll[[j]]$hyperpar[beta1_idx])
      }else{
        med_samps[j,cIdx[3]]<-as.numeric(sampAll[[j]]$latent[beta1_idx])
      }
    }
  }
  print(i)
}

Mediation<-as.data.frame(matrix(NA, ncol=(1+2*3), nrow=length(outcomeIdx)))
colnames(Mediation)<-c("Model", "CDE_med", "CDE_lower", "CDE_upper","NIE_med", "NIE_lower",  "NIE_upper") 
Mediation$Model<-names(med_results)[outcomeIdx]

for(i in 1:nrow(Mediation)){
  if(grepl("cattle", Mediation$Model[i])){
    spp<-"cattle"
  }else{
    spp<-"pigs"
  }
  if(grepl("NDVI",  Mediation$Model[i])){
    mediator<-"NDVI"
  }else{
    mediator<-"LST"
  }
  if(grepl("MEC",  Mediation$Model[i])){
    MEC<-TRUE
  }else{
    MEC<-FALSE
  }
  colIdx_1<-which(grepl(spp, names(med_samps)))
  colIdx_2<-which(grepl(mediator, names(med_samps)))
  if(MEC=="TRUE"){
    colIdx_3<-which(grepl("MEC", names(med_samps))==TRUE)
  }else{
    colIdx_3<-which(grepl("MEC", names(med_samps))==FALSE)
  }
  
  cIdx_a<-colIdx_1[which(colIdx_1%in%colIdx_2)]
  cIdx<-cIdx_a[which(cIdx_a%in%colIdx_3)]
  
  CDE<-exp(median(med_samps[,cIdx[1]]))
  CDE_l<-exp(quantile(probs=0.025, med_samps[,cIdx[1]]))
  CDE_u<-exp(quantile(probs=0.975, med_samps[,cIdx[1]]))
  
  prod<-med_samps[,cIdx[2]] * med_samps[,cIdx[3]]
  NIE<-exp(median(prod))
  NIE_l<-exp(quantile(probs=0.025, prod))
  NIE_u<-exp(quantile(probs=0.975, prod))
  
  Mediation$CDE_med[i]<-CDE
  Mediation$CDE_lower[i]<-CDE_l
  Mediation$CDE_upper[i]<-CDE_u
  Mediation$NIE_med[i]<-NIE
  Mediation$NIE_lower[i]<-NIE_l
  Mediation$NIE_upper[i]<-NIE_u
}

medtable<-xtable(Mediation)

save(Mediation, file="results/regression/South_Sudan/SS_mediation_table.Rdata")



