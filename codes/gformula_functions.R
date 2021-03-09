gformMed = function(tm.vec, HAT, years, type, species, mediator, Y.vars, X.vars, family, lag1vars, lag2vars, lag1names, lag2names, cat.vars){
  
  Y<-Y.vars
  X<-X.vars
  
  res.mean=data.frame(matrix(nrow=length(years), ncol=4))
  names(res.mean)<-c("year",tm.vec)
  res.mean$year<-c(years)
  
  res.total=data.frame(matrix(nrow=length(years), ncol=4))
  names(res.total)<-c("year",tm.vec)
  res.total$year<-c(years)
  
  res.inc=data.frame(matrix(nrow=length(years), ncol=4))
  names(res.inc)<-c("year",tm.vec)
  res.inc$year<-c(years)
  
  pos<-HAT[which(HAT$cases>0),]
  keep<-unique(pos$cluster)
  
  #bootstrap for confidence intervals
  index<-sample(unique(HAT$cluster), size=length(unique(HAT$cluster)), replace=T)
  boot.dat<-as.data.frame(matrix(nrow=1, ncol=ncol(HAT)))
  names(boot.dat)<-names(HAT)
  for(year in years){
    sub<-HAT[which(HAT$year==year),]
    sub<-sub[index,]
    #sub<-sub[which(!is.na(sub$cluster)),]
    boot.dat<-rbind(boot.dat, sub)
  }
  boot.dat<-boot.dat[-1,]
  boot.dat$clus.id<-1:nrow(boot.dat)
  
  #Monte Carlo samples of baseline (not downstream of exposure) confounders
  index<-sample(unique(boot.dat$cluster), size=10000, replace=T)
  add<-which(keep%ni%index) #Keep all clusters that had cases
  index<-c(index,add)
  MC<-as.data.frame(matrix(nrow=1, ncol=ncol(HAT)+1))
  names(MC)<-c(names(HAT), "clus.id")
  for(year in years){
    sub<-boot.dat[which(boot.dat$year==year),]
    sub<-sub[index,]
    MC<-rbind(MC, sub)
  }
  MC<-MC[-1,]
  
  for(tm in tm.vec){
    gmHAT<-boot.dat
    gmMC<-MC
    gmHAT<-gmHAT[which(gmHAT$clus.id%in%gmMC$clus.id),] #Subset to clusters in the MC sample
    
    if(species=="cattle"){
      gmHAT$livestock<-log(gmHAT$density_med_c)
      gmHAT$PO<-log(gmHAT$cattle_POmean)
      
      gmMC$livestock<-log(gmMC$density_med_c)
      gmMC$PO<-log(gmMC$cattle_POmean)
    }else{
      gmHAT$livestock<-log(gmHAT$density_med_p)
      gmHAT$PO<-log(gmHAT$pigs_POmean)
      
      gmMC$livestock<-log(gmMC$density_med_p)
      gmMC$PO<-log(gmMC$pigs_POmean)
    }
    
    if(tm=="a1m1"){
      if(species=="cattle"){
        exp_med<-exp<-log(a1_c)
      }else{
        exp_med<-exp<-log(a1_p)
      }
    }

    mean<-mean(HAT$density_med_c, na.rm=T)
    sd<-sd(HAT$density_med_c, na.rm=T)
    a0_c<-mean
    a1_c<-mean+0.5*mean
    
    mean<-mean(HAT$density_med_c, na.rm=T)
    sd<-sd(HAT$density_med_c, na.rm=T)
    a0_p<-mean
    a1_p<-mean+0.5*mean
    
    if(tm=="a1m0"){
      if(species=="cattle"){
        exp<-log(a1_c) 
        exp_med<-log(a0_c)
      }else{
        exp<-log(a1_p)
        exp_med<-log(a0_p)
      }
    }
    if(tm=="a0m0"){
      if(species=="cattle"){
        exp_med<-exp<-log(a0_c)
      }else{
        exp_med<-exp<-log(a0_p)
      }        
    }
    
    gmMC[,c("livestock", "PO")]<-exp
    
    for(year in years){
      #message(paste0("Year: ", year))
      #Subset data by year
      dat<-gmHAT[which(gmHAT$year==year),]
      datMC<-gmMC[which(gmMC$year==year),]
      
      if(year==2001){
        #Lag confounders
        dat_prev1<-gmHAT[which(gmHAT$year==year-1),]
        dat[,c(lag1names)]<-dat_prev1[,c(lag1vars)]
        
        dat_prev1MC<-gmMC[which(gmMC$year==year-1),]
        datMC[,c(lag1names)]<-dat_prev1MC[,c(lag1vars)]
        
        problem_vars<-c()
        for(c in 1:length(cat.vars)){
          idx<-which(names(dat)==cat.vars[c])
          t<-table(dat[,idx])
          if(dim(t)==1){
            warning(paste0("Zero cells in ", cat.vars[c], " variable, year ", year))
            problem_vars<-c(problem_vars, cat.vars[c])
          }
        }
        #-------------------------#
        # Fit models for t>0
        #-------------------------#
        r<-which(Y%in%c("cases", "livestock"))
        Y.base<-Y[-r]
        X.base<-X[-r]
        family.base<-family[-r]
        
        for(x in 1:length(X.base)){
          idx<-which(X.base[[x]]%in%problem_vars)
          if(length(idx)>0){
            X.base[[x]]<-X.base[[x]][-idx]
          }
        }
        
        datMC$livestock<-datMC$PO<-exp
        
        res.bind<-as.data.frame(matrix(NA, nrow=nrow(datMC), ncol=length(Y.base)))
        colnames(res.bind)<-Y.base
        
        for(var in 1:length(Y.base)){
          #Fit model
          form<-as.formula(paste(Y.base[var], paste(X.base[[var]], collapse=" + "), sep=" ~ "))
          res<-glm (form, data=dat, family=family.base[var])
          
          #Predict
          if(Y.base[var]!=mediator){
            p<-predict(res, newdata=data.frame(datMC), type="response")
          }else{
            datMC$livestock_lag1<-exp_med
            p<-predict(res, newdata=data.frame(datMC), type="response")
          }
          
          res.bind[,var]<-p
          
          #nms<-str_split(names(res.bind), "_")
          #nms<-unlist(lapply(nms, `[[`, 1))
        }
        #names(res.bind)<-nms
        #gmMC[which(gmMC$year==year),c(nms)]<-res.bind[,c(nms)]
        gmMC[which(gmMC$year==year),c(names(res.bind))]<-res.bind
        
      }
      
      if(year>2001){
        #Lag confounders
        dat_prev1<-gmHAT[which(gmHAT$year==year-1),]
        dat[,c(lag1names)]<-dat_prev1[,c(lag1vars)]
        
        dat_prev2<-gmHAT[which(gmHAT$year==year-2),]
        dat[,c(lag2names)]<-dat_prev2[,c(lag2vars)]
        
        dat_prev1MC<-gmMC[which(gmMC$year==year-1),]
        datMC[,c(lag1names)]<-dat_prev1MC[,c(lag1vars)]
        
        dat_prev2MC<-gmMC[which(gmMC$year==year-2),]
        datMC[,c(lag2names)]<-dat_prev2MC[,c(lag2vars)]
        
        for(c in 1:length(cat.vars)){
          idx<-which(names(dat)==cat.vars[c])
          t<-table(dat[,idx])
          if(dim(t)==1){
            warning(paste0("Zero cells in ", cat.vars[c], " variable, year ", year))
          }
        }
        
        #-------------------------#
        # Fit models for t>0
        #-------------------------#
        r<-which(Y%in%c("livestock"))
        Y.test<-Y[-r]
        X.test<-X[-r]
        family.test<-family[-r]
        
        for(x in 1:length(X.test)){
          idx<-which(X.test[[x]]%in%problem_vars)
          if(length(idx)>0){
            X.test[[x]]<-X.test[[x]][-idx]
          }
        }
        
        datMC$livestock<-datMC$PO<-exp
        
        res.bind<-as.data.frame(matrix(NA, nrow=nrow(datMC), ncol=length(Y.test)))
        colnames(res.bind)<-Y.test
        
        for(var in 1:length(Y.test)){
          #Fit model
          form<-as.formula(paste(Y.test[var], paste(X.test[[var]], collapse=" + "), sep=" ~ "))
          res<-glm (form, data=dat, family=family.test[var])
          
          #Predict
          if(Y.test[var]!="cases"){
            if(Y.test[var]!=mediator){
              p<-predict(res, newdata=data.frame(datMC), type="response")
            }else{
              datMC$livestock_lag1<-exp_med
              p<-predict(res, newdata=data.frame(datMC), type="response")
            }
          }else{
            pRisk<-predict(res, newdata=data.frame(datMC), type="response")
            p<-rpois(lambda=median(pRisk, na.rm=T), n=datMC$WorldPop)
          }
          
          res.bind[,var]<-p
          
          nms<-str_split(names(res.bind), "_")
          nms<-unlist(lapply(nms, `[[`, 1))
        }
        
        #names(res.bind)<-nms
        #gmMC[which(gmMC$year==year),c(nms)]<-res.bind[,c(nms)]
        gmMC[which(gmMC$year==year),c(names(res.bind))]<-res.bind
        
      }
    }
    idx<-which(tm.vec==tm)
    res.collapse<-gmMC%>%group_by(year)%>%summarise(mean.cases=mean(cases, na.rm=T), total.cases=sum(cases, na.rm=T),
                                                    pop=sum(WorldPop, na.rm=T))%>% as.data.frame()
    res.collapse$inc<-res.collapse$total.cases/res.collapse$pop
    res.collapse<-res.collapse[which(!is.na(res.collapse$year)),]
    
    res.mean[,idx+1]<-res.collapse$mean.cases
    res.total[,idx+1]<-res.collapse$total.cases
    res.inc[,idx+1]<-res.collapse$inc
    
  }
  
  return(list(incidence=res.inc, mean=res.mean, total=res.total))
}


gform = function(tx.vec, HAT, years, type, species, Y.vars, X.vars, family, test, lag1vars, lag2vars, lag1names, lag2names, cat.vars){
  Y<-Y.vars
  X<-X.vars
  
  mean<-mean(HAT$density_med_c, na.rm=T)
  sd<-sd(HAT$density_med_c, na.rm=T)
  a0_c<-mean
  a1_c<-mean+0.5*mean
  
  mean<-mean(HAT$density_med_c, na.rm=T)
  sd<-sd(HAT$density_med_c, na.rm=T)
  a0_p<-mean
  a1_p<-mean+0.5*mean
  
  valid<-data.frame(matrix(nrow=length(years), ncol=4))
  rownames(valid)<-years
  names(valid)<-c("observed.exposure", "predicted.expoosure", "observed.outcome", "predicted.outcome")
  
  pos<-HAT[which(HAT$cases>0),]
  keep<-unique(pos$cluster)
  
  #bootstrap for confidence intervals
  index<-sample(unique(HAT$cluster), size=length(unique(HAT$cluster)), replace=T)
  boot.dat<-as.data.frame(matrix(nrow=1, ncol=ncol(HAT)))
  names(boot.dat)<-names(HAT)
  for(year in years){
    sub<-HAT[which(HAT$year==year),]
    sub<-sub[index,]
    #sub<-sub[which(!is.na(sub$cluster)),]
    boot.dat<-rbind(boot.dat, sub)
  }
  boot.dat<-boot.dat[-1,]
  boot.dat$clus.id<-1:nrow(boot.dat)
  
  #Monte Carlo samples of baseline (not downstream of exposure) confounders
  index<-sample(unique(boot.dat$cluster), size=10000, replace=T)
  add<-which(keep%ni%index) #Keep all clusters that had cases
  index<-c(index,add)
  MC<-as.data.frame(matrix(nrow=1, ncol=ncol(HAT)+1))
  names(MC)<-c(names(HAT), "clus.id")
  for(year in years){
    sub<-boot.dat[which(boot.dat$year==year),]
    sub<-sub[index,]
    MC<-rbind(MC, sub)
  }
  MC<-MC[-1,]
  
  if(test=="yes"){
    for(tx in tx.vec){
      gmHAT<-boot.dat
      gmMC<-MC
      gmHAT<-gmHAT[which(gmHAT$clus.id%in%gmMC$clus.id),] #Subset to clusters in the MC sample
      
      if(species=="cattle"){
        gmHAT$livestock<-log(gmHAT$density_med_c)
        gmHAT$PO<-log(gmHAT$cattle_POmean)
        
        gmMC$livestock<-log(gmMC$density_med_c)
        gmMC$PO<-log(gmMC$cattle_POmean)
      }else{
        gmHAT$livestock<-log(gmHAT$density_med_p)
        gmHAT$PO<-log(gmHAT$pigs_POmean)
        
        gmMC$livestock<-log(gmMC$density_med_p)
        gmMC$PO<-log(gmMC$pigs_POmean)
      }
      
      if(tx=="a1"){
        if(species=="cattle"){
          exp<-log(a1_c) 
        }else{
          exp<-log(a1_p)
        }
      }
      if(tx=="a0"){
        if(species=="cattle"){
          exp<-log(a0_c) 
        }else{
          exp<-log(a0_p)
        }
      }
      
      gmMC[,c("livestock", "PO")]<-exp
      
      for(year in years){
        #message(paste0("Year: ", year))
        #Subset data by year
        dat<-gmHAT[which(gmHAT$year==year),]
        datMC<-gmMC[which(gmMC$year==year),]
        
        if(year==2001){
          #Lag confounders
          dat_prev1<-gmHAT[which(gmHAT$year==year-1),]
          dat[,c(lag1names)]<-dat_prev1[,c(lag1vars)]
          
          dat_prev1MC<-gmMC[which(gmMC$year==year-1),]
          datMC[,c(lag1names)]<-dat_prev1MC[,c(lag1vars)]
          
          problem_vars<-c()
          for(c in 1:length(cat.vars)){
            idx<-which(names(dat)==cat.vars[c])
            t<-table(dat[,idx])
            if(dim(t)==1){
              warning(paste0("Zero cells in ", cat.vars[c], " variable, year ", year))
              problem_vars<-c(problem_vars, cat.vars[c])
            }
          }
          #-------------------------#
          # Fit models for t>0
          #-------------------------#
          r<-which(Y%in%c("cases", "livestock"))
          Y.base<-Y[-r]
          X.base<-X[-r]
          family.base<-family[-r]
          
          for(x in 1:length(X.base)){
            idx<-which(X.base[[x]]%in%problem_vars)
            if(length(idx)>0){
              X.base[[x]]<-X.base[[x]][-idx]
            }
          }
          
          datMC$livestock<-datMC$PO<-exp
          
          res.bind<-as.data.frame(matrix(NA, nrow=nrow(datMC), ncol=length(Y.base)))
          colnames(res.bind)<-Y.base
          
          for(var in 1:length(Y.base)){
            #Fit model
            form<-as.formula(paste(Y.base[var], paste(X.base[[var]], collapse=" + "), sep=" ~ "))
            res<-glm (form, data=dat, family=family.base[var])
            
            #Predict
            p<-predict(res, newdata=data.frame(datMC), type="response")
            
            res.bind[,var]<-p
            
            #nms<-str_split(names(res.bind), "_")
            #nms<-unlist(lapply(nms, `[[`, 1))
            
            
          }
          #names(res.bind)<-nms
          #gmMC[which(gmMC$year==year),c(nms)]<-res.bind[,c(nms)]
          gmMC[which(gmMC$year==year),c(names(res.bind))]<-res.bind
          
        }
        
        if(year>2001){
          #Lag confounders
          dat_prev1<-gmHAT[which(gmHAT$year==year-1),]
          dat[,c(lag1names)]<-dat_prev1[,c(lag1vars)]
          
          dat_prev2<-gmHAT[which(gmHAT$year==year-2),]
          dat[,c(lag2names)]<-dat_prev2[,c(lag2vars)]
          
          dat_prev1MC<-gmMC[which(gmMC$year==year-1),]
          datMC[,c(lag1names)]<-dat_prev1MC[,c(lag1vars)]
          
          dat_prev2MC<-gmMC[which(gmMC$year==year-2),]
          datMC[,c(lag2names)]<-dat_prev2MC[,c(lag2vars)]
          
          for(c in 1:length(cat.vars)){
            idx<-which(names(dat)==cat.vars[c])
            t<-table(dat[,idx])
            if(dim(t)==1){
              warning(paste0("Zero cells in ", cat.vars[c], " variable, year ", year))
            }
          }
          
          #-------------------------#
          # Fit models for t>0
          #-------------------------#
          r<-which(Y%in%c("livestock"))
          Y.test<-Y[-r]
          X.test<-X[-r]
          family.test<-family[-r]
          
          for(x in 1:length(X.test)){
            idx<-which(X.test[[x]]%in%problem_vars)
            if(length(idx)>0){
              X.test[[x]]<-X.test[[x]][-idx]
            }
          }
          
          datMC$livestock<-datMC$PO<-exp
          
          res.bind<-as.data.frame(matrix(NA, nrow=nrow(datMC), ncol=length(Y.test)))
          colnames(res.bind)<-Y.test
          
          for(var in 1:length(Y.test)){
            
            #Fit model
            form<-as.formula(paste(Y.test[var], paste(X.test[[var]], collapse=" + "), sep=" ~ "))
            res<-glm (form, data=dat, family=family.test[var])
            
            #Predict
            if(Y.test[var]!="cases"){
              p<-predict(res, newdata=data.frame(datMC), type="response")
            }else{
              pRisk<-predict(res, newdata=data.frame(datMC), type="response")
              p<-rpois(lambda=median(pRisk, na.rm=T), n=datMC$WorldPop)
            }
            
            res.bind[,var]<-p
            
            nms<-str_split(names(res.bind), "_")
            nms<-unlist(lapply(nms, `[[`, 1))
          }
          
          #names(res.bind)<-nms
          #gmMC[which(gmMC$year==year),c(nms)]<-res.bind[,c(nms)]
          gmMC[which(gmMC$year==year),c(names(res.bind))]<-res.bind
          
        }
      }
      idx<-which(tx.vec==tx)
      res.collapse<-gmMC%>%group_by(year)%>%summarise(mean.cases=mean(cases, na.rm=T), total.cases=sum(cases, na.rm=T),
                                                      pop=sum(WorldPop, na.rm=T))%>% as.data.frame()
      res.collapse$inc<-res.collapse$total.cases/res.collapse$pop
      res.collapse<-res.collapse[which(!is.na(res.collapse$year)),]
      
      res.mean[,idx+1]<-res.collapse$mean.cases
      res.total[,idx+1]<-res.collapse$total.cases
      res.inc[,idx+1]<-res.collapse$inc
    }
    return(list(incidence=res.inc, mean=res.mean, total=res.total))
  }else{
    gmHAT<-boot.dat
    gmMC<-MC
    gmHAT<-gmHAT[which(gmHAT$clus.id%in%gmMC$clus.id),] #Subset to clusters in the MC sample
    
    if(species=="cattle"){
      gmHAT$livestock<-log(gmHAT$density_med_c)
      gmHAT$PO<-log(gmHAT$cattle_POmean)
      
      gmMC$livestock<-log(gmMC$density_med_c)
      gmMC$PO<-log(gmMC$cattle_POmean)
    }else{
      gmHAT$livestock<-log(gmHAT$density_med_p)
      gmHAT$PO<-log(gmHAT$pigs_POmean)
      
      gmMC$livestock<-log(gmMC$density_med_p)
      gmMC$PO<-log(gmMC$pigs_POmean)
    }
    
    gmHAT.test<-gmHAT.train<-gmHAT
    
    for(year in years){
      #message(paste0("Year: ", year))
      #Subset data by year
      dat<-gmHAT.train[which(gmHAT.train$year==year),]
      datMC<-gmHAT.test[which(gmHAT.test$year==year),]
      
      problem_vars<-c()
      for(c in 1:length(cat.vars)){
        idx<-which(names(dat)==cat.vars[c])
        if(length(idx)>0){
          t<-table(dat[,idx])
          if(dim(t)==1){
            warning(paste0("Zero cells in ", cat.vars[c], " variable, year ", year))
            problem_vars<-c(problem_vars, cat.vars[c])
          }
        }
      }
        
        if(year==2000){
          #--------------------------------#
          # Exposure model (for validation)
          #--------------------------------#
          #Predictors: wealth, , NDVI (removed protected areas and disaster as couldn't fit outcome models d/t zero cells, and d/t overfitting)
          idx<-which(Y=="livestock")
          
          predictors<-X[[idx]]
          r.x<-which(predictors%in%problem_vars)
          if(length(r.x)>0){
            predictors<-predictors[-r.x]
          }
          
          rm.lag<-which(grepl("lag", predictors)==TRUE)
          predictors<-predictors[-rm.lag]
          
          form<-as.formula(paste(Y[idx], paste(predictors, collapse=" + "), sep=" ~ "))
          
          res.livestock<-glm (form, data=dat, family=family[idx])
          
          #-----------------------------------#
          # Predict on MC sample
          #-----------------------------------#
          pLivestock<-predict(res.livestock, newdata=datMC, type="response")
          #pLivestock<-predict(res.livestock, type="response")
          og<-mean(dat$livestock, na.rm=T)
          pred<-mean(pLivestock, na.rm=T)
          rel.diff<-100*((abs(og-pred))/abs(pred))
          #print(paste0(year, ": Prediction is different from observed livestock by ", round(rel.diff,2), "%"))
          
          res.bind<-as.data.frame(cbind(pLivestock, clus.id=datMC$clus.id))
          
          gmHAT.test[which(gmHAT.test$year==year),c("livestock")]<-res.bind[,c("pLivestock")]
          
        }
        if(year==2001){
          #Lag confounders
          dat_prev1<-gmHAT.train[which(gmHAT.train$year==year-1),]
          dat[,c(lag1names)]<-dat_prev1[,c(lag1vars)]
          
          dat_prev1MC<-gmHAT.test[which(gmHAT.test$year==year-1),]
          datMC[,c(lag1names)]<-dat_prev1MC[,c(lag1vars)]
          
          r<-which(Y%in%c("cases"))
          Y.base<-Y[-r]
          X.base<-X[-r]
          family.base<-family[-r]
          
          for(x in 1:length(X.base)){
            idx<-which(X.base[[x]]%in%problem_vars)
            if(length(idx)>0){
              X.base[[x]]<-X.base[[x]][-idx]
            }
          }
          
          res.bind<-as.data.frame(matrix(NA, nrow=nrow(datMC), ncol=length(Y.base)))
          colnames(res.bind)<-Y.base
          
          for(var in 1:length(Y.base)){
            #Fit model
            form<-as.formula(paste(Y.base[var], paste(X.base[[var]], collapse=" + "), sep=" ~ "))
            res<-glm (form, data=dat, family=family.base[var])
            
            #Predict
            if(Y.base[var]!="livestock"){
              p<-predict(res, newdata=data.frame(datMC), type="response")
            }else{
              m<-which(X.base[[var]]%in%Y.base)
              datMC.live<-datMC
              mc.idx<-which(names(datMC.live)%in%X.base[[var]][m])
              rb.idx<-which(names(res.bind)%in%X.base[[var]][m])
              datMC.live[,mc.idx]<-res.bind[,rb.idx]
              p<-predict(res, newdata=data.frame(datMC.live), type="response")
            }
            res.bind[,var]<-p
            
          }
          #nms<-str_split(names(res.bind), "_")
          #nms<-unlist(lapply(nms, `[[`, 1))
          #names(res.bind)<-nms
          #gmMC[which(gmMC$year==year),c(nms)]<-res.bind[,c(nms)]
          gmHAT.test[which(gmHAT.test$year==year),c(names(res.bind))]<-res.bind
        }
        if(year>2001){
          #Lag confounders
          dat_prev1<-gmHAT.train[which(gmHAT.train$year==year-1),]
          dat[,c(lag1names)]<-dat_prev1[,c(lag1vars)]
          
          dat_prev2<-gmHAT.train[which(gmHAT.train$year==year-2),]
          dat[,c(lag2names)]<-dat_prev2[,c(lag2vars)]
          
          dat_prev1MC<-gmHAT.test[which(gmHAT.test$year==year-1),]
          datMC[,c(lag1names)]<-dat_prev1MC[,c(lag1vars)]
          
          dat_prev2MC<-gmHAT.test[which(gmHAT.test$year==year-2),]
          datMC[,c(lag2names)]<-dat_prev2MC[,c(lag2vars)]
          
          res.bind<-as.data.frame(matrix(NA, nrow=nrow(datMC), ncol=length(Y)))
          colnames(res.bind)<-Y
          
          for(var in 1:length(Y)){
            X.test<-X[[var]]
            r.x<-which(X.test%in%problem_vars)
            if(length(r.x)>0){
              X.test<-X.test[-r.x]
            }
            
            #Fit model
            form<-as.formula(paste(Y[var], paste(X.test, collapse=" + "), sep=" ~ "))
            res<-glm (form, data=dat, family=family[var])
            
            #Predict
            if(Y[var]%ni%c("livestock","cases")){
              p<-predict(res, newdata=data.frame(datMC), type="response")
            }
            if(Y[var]=="livestock"){
              m<-which(X.test%in%Y)
              datMC.live<-datMC
              mc.idx<-which(names(datMC.live)%in%X.test[m])
              rb.idx<-which(names(res.bind)%in%X.test[m])
              datMC.live[,mc.idx]<-res.bind[,rb.idx]
              p<-predict(res, newdata=data.frame(datMC.live), type="response")
            }
            if(Y[var]=="cases"){
              pRisk<-predict(res, newdata=data.frame(datMC), type="response")
              p<-rpois(lambda=median(pRisk, na.rm=T), n=datMC$WorldPop)
            }
            res.bind[,var]<-p
          }
          #nms<-str_split(names(res.bind), "_")
          #nms<-unlist(lapply(nms, `[[`, 1))
          #names(res.bind)<-nms
          #gmMC[which(gmMC$year==year),c(nms)]<-res.bind[,c(nms)]
          
          gmHAT.test[which(gmHAT.test$year==year),c(names(res.bind))]<-res.bind
          
        }
      }
      message(paste0("finished all models, natural course iteration ", i))
      
      exp.og<-aggregate(gmHAT.train$livestock, list(gmHAT.train$year), mean)
      exp.og$x<-exp(exp.og$x)
      
      exp.pred<-aggregate(gmHAT.test$livestock, list(gmHAT.test$year), mean)
      exp.pred$x<-exp(exp.pred$x)
      
      out.og<-aggregate(gmHAT.train$cases, list(gmHAT.train$year), mean)
      
      out.pred<-aggregate(gmHAT.test$cases, list(gmHAT.test$year), mean)
      
      valid$observed.exposure<-exp.og$x
      valid$predicted.expoosure<-exp.pred$x
      valid$observed.outcome<-out.og$x
      valid$predicted.outcome<-out.pred$x
      return(valid)
    }
}
