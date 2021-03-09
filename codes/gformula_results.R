rm(list = ls())
'%ni%' <- Negate('%in%')

############################
# Mean across years  
############################
Uganda.years<-c(2006:2018)
Uganda.results<-as.data.frame(matrix(nrow=20, ncol=4))
colnames(Uganda.results)<-c("Model", "Mean", "Lower", "Upper")
Uganda.results$Model<-c("Total_cattle_rhod", "Total_pigs_rhod", "Total_cattle_gamb", "Total_pigs_gamb", 
                        "NDE_cattle_NDVI_rhod", "NDE_pigs_NDVI_rhod", "NDE_cattle_LST_rhod", "NDE_pigs_LST_rhod",
                        "NDE_cattle_NDVI_gamb", "NDE_pigs_NDVI_gamb", "NDE_cattle_LST_gamb", "NDE_pigs_LST_gamb",
                        "NIE_cattle_NDVI_rhod", "NIE_pigs_NDVI_rhod", "NIE_cattle_LST_rhod", "NIE_pigs_LST_rhod",
                        "NIE_cattle_NDVI_gamb", "NIE_pigs_NDVI_gamb", "NIE_cattle_LST_gamb", "NIE_pigs_LST_gamb")

Malawi.years<-c(2003:2014)
Malawi.results<-as.data.frame(matrix(nrow=6, ncol=4))
colnames(Malawi.results)<-c("Model", "Mean", "Lower", "Upper")
Malawi.results$Model<-c("Total_cattle_rhod", "Total_pigs_rhod", 
                        "NDE_cattle_NDVI_rhod", "NDE_pigs_NDVI_rhod", 
                        "NIE_cattle_NDVI_rhod", "NIE_pigs_NDVI_rhod")

DRC.years<-c(2010:2013)
DRC.results<-as.data.frame(matrix(nrow=10, ncol=4))
colnames(DRC.results)<-c("Model", "Mean", "Lower", "Upper")
DRC.results$Model<-c("Total_cattle_gamb", "Total_pigs_gamb", 
                        "NDE_cattle_NDVI_gamb", "NDE_pigs_NDVI_gamb", "NDE_cattle_LST_gamb", "NDE_pigs_LST_gamb",
                        "NIE_cattle_NDVI_gamb", "NIE_pigs_NDVI_gamb", "NIE_cattle_LST_gamb", "NIE_pigs_LST_gamb")


#-----Total effect------#
#Uganda, rhodesiense
load("gform_results/Uganda/MSEs_ug_rhod.Rdata")
MSEs.years<-MSEs[which(MSEs$Year%in%Uganda.years),]

#-Cattle
load("gform_results/Uganda/gform_cattle_inc_rhod.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
Uganda.results$Mean[1]<-mean(RRvec)
Uganda.results$Lower[1]<-quantile(RRvec, prob=0.025)
Uganda.results$Upper[1]<-quantile(RRvec, prob=0.975)

#-Pigs
load("gform_results/Uganda/gform_pigs_inc_rhod.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
Uganda.results$Mean[2]<-mean(RRvec, na.rm=T)
Uganda.results$Lower[2]<-quantile(RRvec, prob=0.025, na.rm=T)
Uganda.results$Upper[2]<-quantile(RRvec, prob=0.975, na.rm=T)

#Uganda, gambiense
load("gform_results/Uganda/MSEs_ug_gamb.Rdata")
MSEs.years<-MSEs[which(MSEs$Year%in%Uganda.years),]
Uganda.pigs.gamb.years<-Uganda.years[Uganda.years%ni%c(2010, 2012)]

#-Cattle
load("gform_results/Uganda/gform_cattle_inc_gamb.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
Uganda.results$Mean[3]<-mean(RRvec)
Uganda.results$Lower[3]<-quantile(RRvec, prob=0.025)
Uganda.results$Upper[3]<-quantile(RRvec, prob=0.975)

#-Pigs
load("gform_results/Uganda/gform_pigs_inc_gamb.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.pigs.gamb.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
Uganda.results$Mean[4]<-mean(RRvec)
Uganda.results$Lower[4]<-quantile(RRvec, prob=0.025)
Uganda.results$Upper[4]<-quantile(RRvec, prob=0.975)

#Malawi
load("gform_results/Malawi/MSEs_mw_rhod.Rdata")
MSEs.years<-MSEs[which(MSEs$Year%in%Malawi.years),]

#-Cattle
load("gform_results/Malawi/gform_cattle_inc_rhod.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Malawi.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
Malawi.results$Mean[1]<-mean(RRvec)
Malawi.results$Lower[1]<-quantile(RRvec, prob=0.025)
Malawi.results$Upper[1]<-quantile(RRvec, prob=0.975)

#-Pigs
load("gform_results/Malawi/gform_pigs_inc_rhod.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Malawi.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
Malawi.results$Mean[2]<-mean(RRvec)
Malawi.results$Lower[2]<-quantile(RRvec, prob=0.025)
Malawi.results$Upper[2]<-quantile(RRvec, prob=0.975)

#DRC
load("gform_results/DRC/MSEs_drc_gamb.Rdata")
MSEs.years<-MSEs[which(MSEs$Year%in%DRC.years),]

#-Cattle
load("gform_results/DRC/gform_cattle_inc_gamb.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%DRC.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
DRC.results$Mean[1]<-mean(RRvec)
DRC.results$Lower[1]<-quantile(RRvec, prob=0.025)
DRC.results$Upper[1]<-quantile(RRvec, prob=0.975)

#-Pigs
load("gform_results/DRC/gform_pigs_inc_gamb.Rdata")
RRvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]$RR<-results.inc[[i]]$a1/results.inc[[i]]$a0
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%DRC.years),]
  results.inc[[i]]$RR[is.infinite(results.inc[[i]]$RR)] <- NA  
  
  med<-median(results.inc[[i]]$RR, na.rm=T)
  RRvec<-c(RRvec, med)
}
print(summary(RRvec))
DRC.results$Mean[2]<-mean(RRvec)
DRC.results$Lower[2]<-quantile(RRvec, prob=0.025)
DRC.results$Upper[2]<-quantile(RRvec, prob=0.975)


#-----Mediation------#

###########
# Uganda
###########

########
#-Rhod-#
########
#-NDVI
#--cattle
load("gform_results/Uganda/gform_cattle_inc_med_NDVI_rhod.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA  
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[5]<-mean(NDEvec)
Uganda.results$Lower[5]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[5]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[13]<-mean(NIEvec)
Uganda.results$Lower[13]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[13]<-quantile(NIEvec, prob=0.975)

#--pigs
load("gform_results/Uganda/gform_pigs_inc_med_NDVI_rhod.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA   
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[6]<-mean(NDEvec)
Uganda.results$Lower[6]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[6]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[14]<-mean(NIEvec)
Uganda.results$Lower[14]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[14]<-quantile(NIEvec, prob=0.975)

#-LST
#--cattle
load("gform_results/Uganda/gform_cattle_inc_med_LST_rhod.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA   
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[7]<-mean(NDEvec)
Uganda.results$Lower[7]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[7]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[15]<-mean(NIEvec)
Uganda.results$Lower[15]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[15]<-quantile(NIEvec, prob=0.975)

#--pigs
load("gform_results/Uganda/gform_pigs_inc_med_LST_rhod.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA   
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[8]<-mean(NDEvec)
Uganda.results$Lower[8]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[8]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[16]<-mean(NIEvec)
Uganda.results$Lower[16]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[16]<-quantile(NIEvec, prob=0.975)

########
#-Gamb-#
########
#-NDVI
#--cattle
load("gform_results/Uganda/gform_cattle_inc_med_NDVI_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA   
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[9]<-mean(NDEvec)
Uganda.results$Lower[9]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[9]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[17]<-mean(NIEvec)
Uganda.results$Lower[17]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[17]<-quantile(NIEvec, prob=0.975)

#--pigs
load("gform_results/Uganda/gform_pigs_inc_med_NDVI_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA   
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[10]<-mean(NDEvec)
Uganda.results$Lower[10]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[10]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[18]<-mean(NIEvec)
Uganda.results$Lower[18]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[18]<-quantile(NIEvec, prob=0.975)

#-LST
#--cattle
load("gform_results/Uganda/gform_cattle_inc_med_LST_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA   
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[11]<-mean(NDEvec)
Uganda.results$Lower[11]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[11]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[19]<-mean(NIEvec)
Uganda.results$Lower[19]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[19]<-quantile(NIEvec, prob=0.975)

#--pigs
load("gform_results/Uganda/gform_pigs_inc_med_LST_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Uganda.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA  
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Uganda.results$Mean[12]<-mean(NDEvec)
Uganda.results$Lower[12]<-quantile(NDEvec, prob=0.025)
Uganda.results$Upper[12]<-quantile(NDEvec, prob=0.975)
Uganda.results$Mean[20]<-mean(NIEvec)
Uganda.results$Lower[20]<-quantile(NIEvec, prob=0.025)
Uganda.results$Upper[20]<-quantile(NIEvec, prob=0.975)

###########
# Malawi
###########
#NDVI only: too much missingness in LST

#-Cattle
load("gform_results/Malawi/gform_cattle_inc_med_NDVI_rhod.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Malawi.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA  
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Malawi.results$Mean[3]<-mean(NDEvec)
Malawi.results$Lower[3]<-quantile(NDEvec, prob=0.025)
Malawi.results$Upper[3]<-quantile(NDEvec, prob=0.975)
Malawi.results$Mean[5]<-mean(NIEvec)
Malawi.results$Lower[5]<-quantile(NIEvec, prob=0.025)
Malawi.results$Upper[5]<-quantile(NIEvec, prob=0.975)

#-Pigs
load("gform_results/Malawi/gform_pigs_inc_med_NDVI_rhod.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%Malawi.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA   
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
Malawi.results$Mean[4]<-mean(NDEvec)
Malawi.results$Lower[4]<-quantile(NDEvec, prob=0.025)
Malawi.results$Upper[4]<-quantile(NDEvec, prob=0.975)
Malawi.results$Mean[6]<-mean(NIEvec)
Malawi.results$Lower[6]<-quantile(NIEvec, prob=0.025)
Malawi.results$Upper[6]<-quantile(NIEvec, prob=0.975)

###########
# DRC
###########

#-NDVI
#--cattle
load("gform_results/DRC/gform_cattle_inc_med_NDVI_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%DRC.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA  
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
DRC.results$Mean[3]<-mean(NDEvec)
DRC.results$Lower[3]<-quantile(NDEvec, prob=0.025)
DRC.results$Upper[3]<-quantile(NDEvec, prob=0.975)
DRC.results$Mean[7]<-mean(NIEvec)
DRC.results$Lower[7]<-quantile(NIEvec, prob=0.025)
DRC.results$Upper[7]<-quantile(NIEvec, prob=0.975)

#--pigs
load("gform_results/DRC/gform_pigs_inc_med_NDVI_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%DRC.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA  
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
DRC.results$Mean[4]<-mean(NDEvec)
DRC.results$Lower[4]<-quantile(NDEvec, prob=0.025)
DRC.results$Upper[4]<-quantile(NDEvec, prob=0.975)
DRC.results$Mean[8]<-mean(NIEvec)
DRC.results$Lower[8]<-quantile(NIEvec, prob=0.025)
DRC.results$Upper[8]<-quantile(NIEvec, prob=0.975)

#-LST
#--cattle
load("gform_results/DRC/gform_cattle_inc_med_LST_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%DRC.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA  
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
DRC.results$Mean[5]<-mean(NDEvec)
DRC.results$Lower[5]<-quantile(NDEvec, prob=0.025)
DRC.results$Upper[5]<-quantile(NDEvec, prob=0.975)
DRC.results$Mean[9]<-mean(NIEvec)
DRC.results$Lower[9]<-quantile(NIEvec, prob=0.025)
DRC.results$Upper[9]<-quantile(NIEvec, prob=0.975)

#--pigs
load("gform_results/DRC/gform_pigs_inc_med_LST_gamb.Rdata")
NIEvec<-c()
NDEvec<-c()
for(i in 1:length(results.inc)){
  results.inc[[i]]<-results.inc[[i]][which(results.inc[[i]]$year%in%DRC.years),]
  results.inc[[i]]$NDE<-results.inc[[i]]$a1m0/results.inc[[i]]$a0m0
  results.inc[[i]]$NIE<-results.inc[[i]]$a1m1/results.inc[[i]]$a1m0
  results.inc[[i]]$NIE[is.infinite(results.inc[[i]]$NIE)] <- NA   
  results.inc[[i]]$NDE[is.infinite(results.inc[[i]]$NDE)] <- NA  
  
  med.NDE<-median(results.inc[[i]]$NDE, na.rm=T)
  med.NIE<-median(results.inc[[i]]$NIE, na.rm=T)
  
  NIEvec<-c(NIEvec, med.NIE) #Median across years for all iterations
  NDEvec<-c(NDEvec, med.NDE) #Median across years for all iterations
}
print(summary(NIEvec))
print(summary(NDEvec))
DRC.results$Mean[6]<-mean(NDEvec)
DRC.results$Lower[6]<-quantile(NDEvec, prob=0.025)
DRC.results$Upper[6]<-quantile(NDEvec, prob=0.975)
DRC.results$Mean[10]<-mean(NIEvec)
DRC.results$Lower[10]<-quantile(NIEvec, prob=0.025)
DRC.results$Upper[10]<-quantile(NIEvec, prob=0.975)

