FAOSTAT<-read.csv("Data/External_validation/FAOSTAT_data_2-17-2021-2.csv")
WB<-read.csv("Data/External_validation/API_SP.POP.TOTL_DS2_en_csv_v2_2055649.csv")

#------------------#
#------Malawi------#
#------------------#

#Livestock
mw<-FAOSTAT[which(FAOSTAT$Area=="Malawi"),]
mw.cattle<-mw[which(mw$Item=="Cattle"),]
mw.pigs<-mw[which(mw$Item=="Pigs"),]
mean(mw.cattle$Value) #1093163
mean(mw.pigs$Value) #2149757

#Humans
mw.h<-WB[which(WB$Country.Name=="Malawi"),]
names<-paste0("X", 2000:2019)
mw.h<-mw.h[,c(names )]
mw.h<-as.numeric(mw.h[1,])
mean(mw.h)

mean(mw.cattle$Value)/mean(mw.h) 
mean(mw.pigs$Value)/mean(mw.h) 

#------------------#
#------Uganda------#
#------------------#

#Livestock
ug<-FAOSTAT[which(FAOSTAT$Area=="Uganda"),]
ug.cattle<-ug[which(ug$Item=="Cattle"),]
ug.pigs<-ug[which(ug$Item=="Pigs"),]

#Humans
ug.h<-WB[which(WB$Country.Name=="Uganda"),]
names<-paste0("X", 2000:2019)
ug.h<-ug.h[,c(names )]
ug.h<-as.numeric(ug.h[1,])

mean(ug.cattle$Value)/mean(ug.h) 
mean(ug.pigs$Value)/mean(ug.h) 

#------------------#
#---South Sudan----#
#------------------#

#Livestock
ss<-FAOSTAT[which(FAOSTAT$Area=="South Sudan"),]
ss.cattle<-ss[which(ss$Item=="Cattle"),]
ss.pigs<-ss[which(ss$Item=="Pigs"),]

#Humans
ss.h<-WB[which(WB$Country.Name=="South Sudan"),]
names<-paste0("X", 2000:2019)
ss.h<-ss.h[,c(names )]
ss.h<-as.numeric(ss.h[1,])

mean(ss.cattle$Value)/mean(ss.h) 
mean(ss.pigs$Value)/mean(ss.h) 

#------------------#
#--------DRC-------#
#------------------#

#Livestock
drc<-FAOSTAT[which(FAOSTAT$Area=="Democratic Republic of the Congo"),]
drc.cattle<-drc[which(drc$Item=="Cattle"),]
drc.pigs<-drc[which(drc$Item=="Pigs"),]

#Humans
drc.h<-WB[which(WB$Country.Name=="Congo, Dem. Rep."),]
names<-paste0("X", 2008:2015)
drc.h<-drc.h[,c(names )]
drc.h<-as.numeric(drc.h[1,])

mean(drc.cattle$Value)/mean(drc.h) 
mean(drc.pigs$Value)/mean(drc.h) 