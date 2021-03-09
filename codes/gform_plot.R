rm(list = ls())
library(ggplot2)
library(dplyr)

country<-"DRC"
type<-"gamb"
sim="false"

folderData<-paste0("plots/", country, "/")

if(file.exists(folderData)==FALSE){
  dir.create(file.path(folderData))
}

'%ni%'<-Negate('%in%')

if(country=="Uganda"){
  country.short<-"ug"
}
if(country=="Malawi"){
  country.short<-"mw"
}
if(country=="DRC"){
  country.short<-"drc"
}

if(sim=="true"){
  load(paste0("gform_results/", country, "/natural_course_pigs_sim_", country.short,"_", type, ".Rdata"))
}else{
  load(paste0("gform_results/", country, "/natural_course_pigs_", country.short, "_", type, ".Rdata"))
}
length(validation)

res.df<-as.data.frame(matrix(NA, nrow=length(validation)*dim(validation[[1]])[1], ncol=3))
names(res.df)<-c("obs.outcome", "pred.outcome", "Year")

counter=1
for(i in 1:length(validation)){
  valid<-validation[[i]]
  idx<-c(counter:(counter+dim(valid)[1]-1))
  res.df$obs.outcome[idx]<-valid$observed.outcome
  res.df$pred.outcome[idx]<-valid$predicted.outcome
  res.df$Year[idx]<-rownames(valid)
  counter=counter+dim(valid)[1]
}

res.df<-res.df[which(res.df$Year>2001),]

pigs<-res.df

v.test<-res.df%>%group_by(Year)%>%summarise(med.obs=median(obs.outcome, na.rm=T), med.pred=median(pred.outcome, na.rm=T),
                                            sd.pred=sd(pred.outcome, na.rm=T))%>% as.data.frame()

if(sim=="true"){
  png(paste0("plots/", country, "/natural_course_pigs_sim_", country.short, "_", type, ".png"), height=600, width=1000)
}else{
  png(paste0("plots/", country, "/natural_course_pigs_", country.short, "_", type, ".png"), height=600, width=1000)
}
ggplot(v.test, aes(x=med.obs, y=med.pred, color=Year)) + geom_point()  + geom_pointrange(aes(ymin = med.pred-sd.pred, ymax = med.pred+sd.pred)) + geom_abline(intercept=0,slope=1, lwd=0.25)+ geom_text(label=v.test$Year, nudge_x=0.00045, nudge_y=0.01, size=6)+ xlab("Observed") + ylab("Predicted")+ggtitle("Observed versus predicted (100 simulations) median cases by year, pig model")+theme_bw(base_size = 14)
dev.off()

if(sim=="true"){
  load(paste0("gform_results/", country, "/natural_course_cattle_sim_", country.short,"_", type, ".Rdata"))
}else{
  load(paste0("gform_results/", country, "/natural_course_cattle_", country.short, "_", type, ".Rdata"))
}

length(validation)

res.df<-as.data.frame(matrix(NA, nrow=length(validation)*dim(validation[[1]])[1], ncol=3))
names(res.df)<-c("obs.outcome", "pred.outcome", "Year")

counter=1
for(i in 1:length(validation)){
  valid<-validation[[i]]
  idx<-c(counter:(counter+dim(valid)[1]-1))
  res.df$obs.outcome[idx]<-valid$observed.outcome
  res.df$pred.outcome[idx]<-valid$predicted.outcome
  res.df$Year[idx]<-rownames(valid)
  counter=counter+dim(valid)[1]
}

res.df<-res.df[which(res.df$Year>2001),]

cattle<-res.df

v.test<-res.df%>%group_by(Year)%>%summarise(med.obs=median(obs.outcome, na.rm=T), med.pred=median(pred.outcome, na.rm=T),
                                            sd.pred=median(pred.outcome, na.rm=T))%>% as.data.frame()

if(sim=="true"){
  png(paste0("plots/", country, "/natural_course_cattle_sim_", country.short, "_", type, ".png"), height=600, width=1000)
}else{
  png(paste0("plots/", country, "/natural_course_cattle_", country.short, "_", type, ".png"), height=600, width=1000)

}
ggplot(v.test, aes(x=med.obs, y=med.pred, color=Year)) + geom_point()  + geom_pointrange(aes(ymin = med.pred-sd.pred, ymax = med.pred+sd.pred)) + geom_abline(intercept=0,slope=1, lwd=0.25)+ geom_text(label=v.test$Year, nudge_x=0.00075, nudge_y=0.0005, size=6)+ xlab("Observed") + ylab("Predicted")+ggtitle("Observed versus predicted (100 simulations) median cases by year, cattle model")+theme_bw(base_size = 14)
dev.off()

#-------MSEs-------#
MSEs<-as.data.frame(matrix(NA, nrow=length(unique(res.df$Year)), ncol=3))
colnames(MSEs)<-c("Year", "MSE_cattle", "MSE_pigs")
MSEs$Year<-unique(res.df$Year)

for(y in unique(res.df$Year)){
  sub<-cattle[which(cattle$Year==y),]
  mse<-median((sub$pred.outcome-sub$obs.outcome)^2, na.rm=T)
  idx<-which(MSEs$Year==y)
  MSEs$MSE_cattle[idx]<-mse
  
  sub<-pigs[which(pigs$Year==y),]
  mse<-median((sub$pred.outcome-sub$obs.outcome)^2, na.rm=T)
  idx<-which(MSEs$Year==y)
  MSEs$MSE_pigs[idx]<-mse
}

save(file=paste0("gform_results/", country, "/MSEs_", country.short, "_", type, ".Rdata"), MSEs)
