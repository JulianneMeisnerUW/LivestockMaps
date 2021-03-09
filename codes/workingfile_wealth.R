library(parallel)
folderData<-"For_Brain/"

country<-"DRC"
if(country=="Uganda"){
n<-11
}
if(country=="Malawi"){
n<-15
}
if(country=="DRC"){
n<-3
}

message("Running models")
n.cores<-8
clus<-makeCluster(n.cores, outfile='models_wealth.txt')
source("codes/models_wealth.R")
parLapply(clus, 1:5, modelFn)
stopCluster(clus)

###############################
#   Validation
###############################
n.cores <- 8
clus <- makeCluster(n.cores, outfile = 'validOut_wealth.txt')
source('codes/validation_wealth.R')
message("Validation")
parLapply(clus, 1:n, validFn) 
stopCluster(clus)



