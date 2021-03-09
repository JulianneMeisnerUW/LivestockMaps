rm(list=ls())

library(ggplot2)
library(dplyr)
library(xtable)

country<-"Uganda"
type<-"rhod"
intervention<-"int_1tt_1itc"

if(type=="gamb"|country=="DRC"){
  underreport<-0.65
}
if(type=="rhod"){
  underreport<-0.083
}

finalyear<-2030
startyear<-2005
  
if(country=="Uganda"|country=="DRC"){
  if(country=="Uganda"){
    country.short<-"ug"
  }else{
    country.short<-"drc"
  }
  if(intervention=="base_case"){
    output<-read.csv(paste0("trajectories_", country.short, "_", type, ".csv"), header=F)
  }else{
    output<-read.csv(paste0("trajectories_", country.short, "_", type, "_", intervention, ".csv"), header=F)
  }
  data.years<-c(2005:2018)
  if(country=="Uganda"&type=="rhod"){
    observed<-read.csv("HAT_ug_Tbr.csv")
  }
  if(country=="Uganda"&type=="gamb"){
    observed<-read.csv("HAT_ug_Tbg.csv")
  }
  if(country=="DRC"){
    observed<-read.csv("HAT_drc.csv")
  }
}

if(country=="South_Sudan"){
  country.short<-"ss"
  if(intervention=="base_case"){
    output<-read.csv(paste0("trajectories_", country.short, ".csv"), header=F)
  }else{
    output<-read.csv(paste0("trajectories_", country.short, "_", intervention, ".csv"), header=F)
  }
  data.years<-c(2005:2014)
  observed<-read.csv("HAT_ss.csv")
}

if(country=="Malawi"){
  country.short<-"mw"
  if(intervention=="base_case"){
    output<-read.csv(paste0("trajectories_", country.short, ".csv"), header=F)
  }else{
    output<-read.csv(paste0("trajectories_", country.short, "_", intervention, ".csv"), header=F)
  }
  data.years<-c(2005:2014)
  observed<-read.csv("HAT_mw.csv")
}
  
rownames(output)<-output[,1]
output<-output[,-1]
colnames(output)<-paste("t", 1:ncol(output), sep=".")
data<-as.data.frame(t(as.matrix(output)))

reps<-10 # number of simulations run
if(type=="rhod"){
  #n.compartments<-20 # number of observable compartments
  species<-c("human", "tsetse", "reservoir-sd", "reservoir-rd", "reservoir-sw", "reservoir-rw")
}else{
  #n.compartments<-14 # number of observable compartments
  species<-c("human", "tsetse", "reservoir")
}
n.compartments<-ncol(data)/reps

results<-list()
for(i in 1:reps){
  idx<-seq(i, (n.compartments)*reps, reps)
  r<-data[,idx]
  r$time<-1:nrow(r)
  r$new_infections<-rep(NA, nrow(r))
  for(j in 2:nrow(r)){
    last<-r$`human-infection-cumulative`[j-1]
    current<-r$`human-infection-cumulative`[j]
    r$new_infections[j]<-current-last
  }
  results<-c(results, list(r))
}
#Time is in days, so run for 10 years (3,650 time points)

#Get median over all runs
median<-as.data.frame(matrix(NA, nrow=nrow(results[[1]]), ncol=ncol(results[[1]])))
colnames(median)<-colnames(results[[1]])
rownames(median)<-rownames(results[[1]])

for(i in 1:ncol(median)){
  dat<-matrix(NA, nrow=nrow(median), ncol=length(results))
  for(j in 1:length(results)){
    dat[,j]<-results[[j]][,i]
  }
  med<-apply(X=dat, MARGIN=1, FUN=median, na.rm=TRUE)
  median[,i]<-med
}

#Get 2.5%ile over all runs
lower<-as.data.frame(matrix(NA, nrow=nrow(results[[1]]), ncol=ncol(results[[1]])))
colnames(lower)<-colnames(results[[1]])
rownames(lower)<-rownames(results[[1]])

for(i in 1:ncol(lower)){
  dat<-matrix(NA, nrow=nrow(lower), ncol=length(results))
  for(j in 1:length(results)){
    dat[,j]<-results[[j]][,i]
  }
  low<-apply(X=dat, MARGIN=1, FUN=quantile, probs=0.025, na.rm=TRUE)
  lower[,i]<-low
}

#Get 97.5%ile over all runs
upper<-as.data.frame(matrix(NA, nrow=nrow(results[[1]]), ncol=ncol(results[[1]])))
colnames(upper)<-colnames(results[[1]])
rownames(upper)<-rownames(results[[1]])

for(i in 1:ncol(upper)){
  dat<-matrix(NA, nrow=nrow(upper), ncol=length(results))
  for(j in 1:length(results)){
    dat[,j]<-results[[j]][,i]
  }
  up<-apply(X=dat, MARGIN=1, FUN=quantile, probs=0.975, na.rm=TRUE)
  upper[,i]<-up
}

years<-nrow(median)/365
#----Plots----#

#----By species----#
lims<-list()
for(spp in species){
  idx<-which(grepl(spp, names(results[[1]]))==TRUE)
  nidx<-which(grepl("cumulative", names(results[[1]]))==TRUE)
  idx<-idx[! idx %in% nidx]
  agg<-as.vector(as.matrix(results[[1]][,idx]))
  lims<-c(lims, list(agg))
}

compartments<-colnames(results[[1]])
compartments<-compartments[!compartments%in%c("time")]
rescol<-c("grey", "pink", "lightblue", "seagreen3", "lightgoldenrod2")
medcol<-c("black", "red", "blue", "seagreen4", "lightgoldenrod4")
compnames<-c("Susceptible" , "Exposed" ,"Infectious, stage 1", "Infectious, stage 2", "Recovered", "Infectious", "Non-susceptible")
#infection 2 = seagreen3/seagreen4, recovered = lightgoldenrod2/lightgoldenrod4

xtick<-seq(0, nrow(results[[1]]), 365)

for(spp in 1:length(species)){
  compartments.spp<-which(grepl(species[spp], compartments))
  comps<-compartments[compartments.spp]
  
  if(species[spp]=="human"){
    cnames<-compnames[1:5]
    rc<-rescol
    mc<-medcol
  }
  if(species[spp]=="tsetse"){
    rc<-rescol[c(1:3,5)]
    mc<-medcol[c(1:3,5)]
    cnames<-compnames[c(1,2,6,7)]
  }
  if(grepl("reservoir", species[spp])==TRUE){
    rc<-rescol[c(1:3)]
    mc<-medcol[c(1:3)]
    cnames<-compnames[c(1,2,6)]
  }
  
  if(country!="Uganda"&country!="DRC"){
    png(paste0("plots/", country, "/", intervention, "/", species[spp], "_", country.short, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
  }else{
    png(paste0("plots/", country, "/", type, "/", intervention, "/",species[spp], "_", country.short, "_", type, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
  }
  
  plot(results[[1]]$time, results[[1]][,compartments.spp[1]], type="l", xlab="Year", ylab="N", col=rc[1], lwd=0.5, xaxt='n', ylim=c(min(lims[[spp]]), max(lims[[spp]])), xaxs="i")
  axis(side=1, at=xtick, labels=c(2014:(2014+years)))
  for(i in 2:length(results)){
    lines(results[[i]]$time, results[[i]][,compartments.spp[1]], type="l", xlab="Day", ylab="N", col=rc[1], lwd=0.5, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
  }
  #lines(median$time, median[,compartments.spp[1]], type="l", xlab="Day", ylab="N", col=mc[1], lwd=1, xaxs="i")
  
  for(j in 2:length(comps)){
    for(i in 1:length(results)){
      lines(results[[i]]$time, results[[i]][,compartments.spp[j]], type="l", xlab="Day", ylab="N", col=rc[j], lwd=0.5, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
    }
    #lines(median$time, median[,compartments.spp[j]], type="l", xlab="Day", ylab="N", col=mc[j], lwd=1, xaxs="i")
  }
  abline(v=c(365*1:(years-1)), lty=2, lwd=0.25)
  
  legend("topright", c(cnames), col=mc, lty=1, lwd=2, seg.len=3)
  
  dev.off()
}

#----By compartment x species----#
for(c in 1:length(compartments)){
  idx<-which(names(results[[1]])==compartments[c])
  
  if(country!="Uganda"&country!="DRC"){
    png(paste0("plots/", country, "/", intervention, "/", compartments[c], "_", country.short, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
  }else{
    png(paste0("plots/", country, "/", type, "/", intervention, "/",compartments[c], "_", country.short, "_", type, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
  }
  
  plot(results[[1]]$time, results[[1]][,idx], type="l", xlab="Year", ylab="N", xaxt='n', col="gray", lwd=0.5, xaxs="i")
  axis(side=1, at=xtick, labels=c(2014:(2014+years)))
  for(i in 2:length(results)){
    lines(results[[i]]$time, results[[i]][,idx], type="l", xlab="Year", ylab="N", col="gray", lwd=0.5, xaxs="i")
  }
  lines(median$time, median[,idx], type="l", xlab="Year", ylab="N", col="black", lwd=1, xaxs="i")
  abline(v=c(365*1:(years-1)), lty=2, lwd=0.25)
  dev.off()
}

#----Humans, susceptibles removed----#
ticks<-seq(0, nrow(results[[1]]), 365*5)
labs<-seq(2005, 2005+years, 5)

colIdx<-which(names(results[[1]])%in%c("human-infectious-one", "human-infectious-two"))
for(i in 1:length(results)){
  results[[i]]$`human-infec`<-rowSums(results[[i]][,colIdx])
}

median$`human-infec`<-rowSums(median[,colIdx])
lower$`human-infec`<-rowSums(lower[,colIdx])
upper$`human-infec`<-rowSums(upper[,colIdx])

idx<-c(which(names(median)%in%c("human-exposed", "human-infec")))
ns_lims<-c(as.vector(as.matrix(median[,idx])),as.vector(as.matrix(lower[,idx])),as.vector(as.matrix(upper[,idx])))

if(country!="Uganda"&country!="DRC"){
  png(paste0("plots/", country, "/", intervention, "/humans_nosuscep_", country.short, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}else{
  png(paste0("plots/", country, "/", type, "/", intervention, "/humans_nosuscep_", country.short, "_", type, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}

colors <- c("Exposed" = "red", "Infected" = "blue")

p1 = ggplot(median, aes(time)) + 
  geom_line(aes(x=time, y=`human-exposed`, colour="Exposed")) + geom_ribbon(aes(ymin=lower$`human-exposed`, ymax=upper$`human-exposed`), alpha=0.2, fill="red") +theme(text = element_text(size=20)) 
p2 = p1 + geom_line(aes(x=time, y=`human-infec`, colour="Infected")) + geom_ribbon(aes(ymin=lower$`human-infec`, ymax=upper$`human-infec`), alpha=0.2, fill="blue")
p3 = p2 + ylab("N")+ theme_bw()+theme(text = element_text(size=20))  + xlab("Year")
p4 = p3 + scale_x_continuous(breaks=ticks,labels=labs, expand = c(0, 0))
p5 = p4 + theme(legend.title = element_blank())
plot(p5)

#plot(median$time, median$`human-exposed`, type="l", xlab="Year", ylab="N", col="red", xaxt='n', ylim=c(min(ns_lims), max(ns_lims)), xaxs="i")
#axis(side=1, at=xtick, labels=c(2014:(2014+years)))
#lines(lower$time, lower$`human-exposed`, type="l", xlab="", ylab="", col="pink", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
#lines(upper$time, upper$`human-exposed`, type="l", xlab="", ylab="", col="pink", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))

#lines(median$time, median$`human-infec`, type="l", xlab="", ylab="", col="blue", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
#lines(lower$time, lower$`human-infec`, type="l", xlab="", ylab="", col="lightblue", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
#lines(upper$time, upper$`human-infec`, type="l", xlab="", ylab="", col="lightblue", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))

#legend("topright", c("Exposed", "Infectious (sum stages 1 and 2)"), col=c("red", "blue"), lty=1, lwd=2, seg.len=3)

dev.off()

#----Tsetse, proportion----#
colIdx<-which(grepl("tsetse", names(results[[1]]))==TRUE)
for(i in 1:length(results)){
  results[[i]]$`tsetse-total`<-rowSums(results[[i]][,colIdx])
  results[[i]]$`tsetse-exposed-prop`<-results[[i]]$`tsetse-exposed`/results[[i]]$`tsetse-total`
  results[[i]]$`tsetse-infec-prop`<-results[[i]]$`tsetse-infectious`/results[[i]]$`tsetse-total`
}

median$`tsetse-total`<-rowSums(median[,colIdx])
median$`tsetse-exposed-prop`<-median$`tsetse-exposed`/median$`tsetse-total`
median$`tsetse-infec-prop`<-median$`tsetse-infectious`/median$`tsetse-total`

lower$`tsetse-total`<-rowSums(lower[,colIdx])
lower$`tsetse-exposed-prop`<-lower$`tsetse-exposed`/lower$`tsetse-total`
lower$`tsetse-infec-prop`<-lower$`tsetse-infectious`/lower$`tsetse-total`

upper$`tsetse-total`<-rowSums(upper[,colIdx])
upper$`tsetse-exposed-prop`<-upper$`tsetse-exposed`/upper$`tsetse-total`
upper$`tsetse-infec-prop`<-upper$`tsetse-infectious`/upper$`tsetse-total`

idx<-which(names(results[[1]])%in% c("tsetse-exposed-prop", "tsetse-infec-prop"))
ts_lims<-c(as.vector(as.matrix(upper[,idx])), as.vector(as.matrix(lower[,idx])))

if(country!="Uganda"&country!="DRC"){
  png(paste0("plots/", country, "/", intervention, "/tsetse_prop_", country.short, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}else{
  png(paste0("plots/", country, "/", type, "/", intervention, "/tsetse_prop_", country.short, "_", type, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}

colors <- c("Exposed" = "red", "Infected" = "blue")

p1 = ggplot(median, aes(time)) + 
  geom_line(aes(x=time, y=`tsetse-exposed-prop`, colour="Exposed")) + geom_ribbon(aes(ymin=lower$`tsetse-exposed-prop`, ymax=upper$`tsetse-exposed-prop`), alpha=0.2, fill="red")+theme(text = element_text(size=20)) 
p2 = p1 + geom_line(aes(x=time, y=`tsetse-infec-prop`, colour="Infected")) + geom_ribbon(aes(ymin=lower$`tsetse-infec-prop`, ymax=upper$`tsetse-infec-prop`), alpha=0.2, fill="blue")
p3 = p2 + ylab("Proportion of total population")+ theme_bw()+theme(text = element_text(size=20))  + xlab("Year")
p4 = p3 + scale_x_continuous(breaks=ticks,labels=labs, expand = c(0, 0))
p5 = p4 + theme(legend.title = element_blank())
plot(p5)

#plot(median$time, median$`tsetse-exposed-prop`, type="l", xlab="Year", ylab="Proportion of total population", col="red", xaxt='n', ylim=c(min(ts_lims), max(ts_lims)), xaxs="i")
#axis(side=1, at=xtick, labels=c(2014:(2014+years)))
#lines(lower$time, lower$`tsetse-exposed-prop`, type="l", xlab="", ylab="", col="pink", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
#lines(upper$time, upper$`tsetse-exposed-prop`, type="l", xlab="", ylab="", col="pink", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))

#lines(median$time, median$`tsetse-infec-prop`, type="l", xlab="", ylab="", col="blue", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
#lines(lower$time, lower$`tsetse-infec-prop`, type="l", xlab="", ylab="", col="lightblue", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))
#lines(upper$time, upper$`tsetse-infec-prop`, type="l", xlab="", ylab="", col="lightblue", lty=2, xaxs="i", xaxt="none",axis(1, at=seq(0, nrow(results[[1]]), 365), labels=c(2014:(2014+years))))

#legend("topright", c("Exposed", "Infectious"), col=c("red", "blue"), lty=1, lwd=2, seg.len=3)

dev.off()

#----Distribution of time to elimination----#
elim.year=c()
elim.bin=c()#1 = elimination achieved by 2030, 0 = no elimination achieved
infec.vec<-c()
results.year<-list()

observed.collapse<-observed%>%group_by(Year)%>%summarize(new_infections=sum(New_HAT_cases, na.rm=T))%>%as.data.frame()
observed.collapse<-observed.collapse[which(observed.collapse$Year%in%data.years),]

day<-(finalyear-startyear)*365

for(i in 1:length(results)){
  infections<-results[[i]]$new_infections
  results[[i]]$year<-as.numeric(substr(((results[[i]]$time/365)+startyear), 1,4))
  
  if(sum(infections, na.rm=T)==0){
    t<-startyear
  }else{
    t<-max(which(infections>0))
  }
  
  if(t<=day){
    eb<-1
  }else{
    eb<-0
  }
  
  t.year<-round(t/365)+startyear
  if(eb==1){
    elim.year=c(elim.year, t.year)
  }else{
    elim.year=c(elim.year, 3000)
  }
  
  elim.bin<-c(elim.bin, eb)

  results.collapse<-results[[i]]%>%group_by(year)%>%summarize(predicted=sum(new_infections, na.rm=T))%>%as.data.frame
  
  idx<-which(results.collapse$year==finalyear)
  inf<-results.collapse$predicted[idx]
  infec.vec<-c(infec.vec, inf)
  
  results.collapse<-results.collapse[which(results.collapse$year%in%data.years),]
  results.collapse$observed<-observed.collapse$new_infections
  results.collapse$predicted.adjusted<-round(results.collapse$predicted*underreport,0)
  results.collapse$diff<-results.collapse$predicted.adjusted-results.collapse$observed
  results.year<-c(results.year, list(results.collapse))
}

lastyear<-max(results[[1]]$year)
prop<-(sum(elim.bin)/length(elim.bin)*100)

if(country!="Uganda"&country!="DRC"){
  png(paste0("plots/", country, "/", intervention, "/time_elimination_", country.short, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}else{
  png(paste0("plots/", country, "/", type, "/", intervention, "/time_elimination_", country.short, "_", type, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}

qplot(elim.year, geom="histogram")+coord_cartesian(xlim=c(startyear, lastyear)) +geom_vline(xintercept=2030,linetype="dashed")+ theme_bw()+theme(text = element_text(size=20))+xlab(paste0("Elimination year (EOT observed in ", prop, "% of runs)"))
dev.off()

#----Histogram of new transmissions by 2030----#

if(country!="Uganda"&country!="DRC"){
  png(paste0("plots/", country, "/", intervention, "/histogram_", country.short, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}else{
  png(paste0("plots/", country, "/", type, "/", intervention, "/histogram_", country.short, "_", type, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
}

qplot(infec.vec, geom="histogram") + theme_bw()+xlab(paste0("Total human infections (", finalyear, "), ", reps, " runs")) + ylab("Runs")+theme(text = element_text(size=20)) 
dev.off()

print(paste0("Zero cases by 2030 in ", ((length(which(elim.bin==1))/length(elim.bin))*100), "% of runs"))

#----Compare with observed data----#
if(intervention=="base_case"){
  years<-results.year[[1]]$year
  obs.results<-as.data.frame(matrix(NA, nrow=reps*length(years), ncol=3))
  colnames(obs.results)<-c("iteration", "year", "difference")
  
  idx=1
  for(i in 1:length(results.year)){
    results<-results.year[[i]]
    results<-results[,c("year", "diff")]
    obs.results$year[idx:(idx+(length(years))-1)]<-results$year
    obs.results$iteration[idx:(idx+(length(years))-1)]<-i
    obs.results$difference[idx:(idx+(length(years))-1)]<-results$diff
    
    idx=idx+length(years)
  }
  
  breaks<-seq(min(obs.results$year),max(obs.results$year), 2)
  labs<-breaks
  obs.results$Year<-as.factor(obs.results$year)
  
  if(country!="Uganda"&country!="DRC"){
    png(paste0("plots/", country, "/", intervention, "/validation_", country.short, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
    print(ggplot(obs.results, aes(x=Year, y=difference))+ theme_bw()+theme(text = element_text(size=20)) + geom_jitter(width=0.1, alpha=0.25, size=3)+ ylab("Difference in reported cases")+scale_x_discrete(name=paste0("Year"), breaks=breaks, labels=labs))
    dev.off()
  }else{
    png(paste0("plots/", country, "/", type, "/", intervention, "/validation_", country.short, "_", type, "_", intervention, ".png"), height=8, width=13, units="in", res=300)
    print(ggplot(obs.results, aes(x=Year, y=difference)) + theme_bw()+theme(text = element_text(size=20)) + geom_jitter(width=0.1, alpha=0.25, size=3)+ ylab("Difference in reported cases")+scale_x_discrete(name=paste0("Year"), breaks=breaks, labels=labs)) 
    dev.off()
 }
}

#----Save results----#
results.list<-list("results.year"=results.year, "infec.vec"=infec.vec, "elim.year"=elim.year, "elim.bin"=elim.bin)
if(country=="Uganda"|country=="DRC"){
  save(file=paste0("plots/", country, "/", type, "/", intervention, "/results_", intervention, ".Rdata"), results.list)
}else{
  save(file=paste0("plots/", country, "/", intervention, "/results_", intervention, ".Rdata"), results.list)
}


#---------------#
# Overall table#
#---------------#
folders<-list.files("plots")
folderslev2<-c()
for(i in 1:length(folders)){
  if(folders[i]=="DRC"|folders[i]=="Uganda"){
    folders_new<-paste(folders[i], list.files(paste0("plots/", folders[i])), sep="/")
  }else{
    folders_new<-folders[i]
  }
  folderslev2<-c(folderslev2, folders_new)
}
folderslev3<-c()
tests<-c()
for(i in 1:length(folderslev2)){
  folds<-list.files(paste0("plots/", folderslev2[i]))
  keep<-which(grepl("redundant", folds)==FALSE)
  folds<-folds[keep]
  keep<-which(grepl("old", folds)==FALSE)
  folds<-folds[keep]
  res<-paste(folderslev2[i], folds, sep="/")
  folderslev3<-c(folderslev3, res)
  tests<-c(tests, folds)
}

tests2<-paste0("results_", tests, ".Rdata")

filenames<-paste(folderslev3, tests2, sep="/")
names.nice<-c(rep("drcf_", 5), rep("drcq_", 5), rep("mw_", 6), rep("ss_", 5), rep("ugg_", 5),  rep("ugr_", 6))
names.nice<-paste0(names.nice, tests)

overall.results<-as.data.frame(matrix(NA, nrow=length(filenames), ncol=5))
colnames(overall.results)<-c("Model", "EOT_yr_mean", "EOT_yr_sd", "Cases2030_mean", "Cases2030_sd")
overall.results$Model<-names.nice

for(i in 1:length(filenames)){
  load(paste0("plots/", filenames[i]))
  overall.results$Cases2030_mean[i]<-mean(results.list$infec.vec, na.rm=T)
  overall.results$Cases2030_sd[i]<-sd(results.list$infec.vec, na.rm=T)
  
  for(j in 1:length(results.list$elim.bin)){
    if(results.list$elim.bin[j]==0){
      results.list$elim.year[j]<-NA
    }
  }
  
  overall.results$EOT_yr_mean[i]<-mean(results.list$elim.year, na.rm=T)
  overall.results$EOT_yr_sd[i]<-sd(results.list$elim.year, na.rm=T)
}

overall.results<-xtable(overall.results)

