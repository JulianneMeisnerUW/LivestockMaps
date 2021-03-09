library(parallel)
country<-"Malawi"
if(country=="Uganda"){
 country.short<-"ug"
 n<-10
}
if(country=="Malawi"){
 country.short<-"mw"
 n<-12
}
if(country=="DRC"){
 country.short<-"drc"
 n<-2
}

folderData<-paste0("results/livestock_mapping/", country, "/")

if(file.exists(folderData)==FALSE){
        dir.create(file.path(folderData))
}

###############################
#   Models
###############################
message("Running models, cattle")
n.cores<-6
clus<-makeCluster(n.cores, outfile='models_cattle.txt')
if(country!="DRC"){
 source("codes/models_cattle.R")
}else{
 source("codes/models_cattle_simple.R")
}
ST.results.cattle<-parLapply(clus,1:6, modelFn)
save(ST.results.cattle, file= paste0(folderData, "cattle_results/cattle_results_", country.short, ".Rdata"))
stopCluster(clus)

message("Running models, pigs")
n.cores<-6
clus<-makeCluster(n.cores, outfile='models_pigs.txt')
if(country!="DRC"){
        source("codes/models_pigs.R")
}else{
        source("codes/models_cattle_pigs.R")
}
ST.results.pigs<-parLapply(clus,1:6, modelFn)
save(ST.results.pigs, file= paste0(folderData, "pig_results/pig_results_", country.short, ".Rdata"))
stopCluster(clus)


###############################
#   Random field prediction
###############################
message("Predictions, cattle")
n.cores<-6
clus<-makeCluster(n.cores, outfile='predict_cattle.txt')
source("codes/predict_cattle.R")
Predict_nST.results.cattle<-parLapply(clus, 1:6, PredictFn)
save(Predict_nST.results.cattle, file= paste0(folderData, "cattle_results/prediction_cattle_", country.short, ".Rdata"))
stopCluster(clus)

message("Predictions, pigs")
n.cores<-6
clus<-makeCluster(n.cores, outfile='predict_pigs.txt')
source("codes/predict_pigs.R")
Predict_nST.results.pigs<-parLapply(clus, 1:6, PredictFn)
save(Predict_nST.results.pigs, file= paste0(folderData, "pig_results/prediction_pigs_", country.short, ".Rdata"))
stopCluster(clus)

###############################
#   Validation
###############################
n.cores <- 8
clus <- makeCluster(n.cores, outfile = 'validOut.txt')
source('codes/validation.R')
message("Validation, cattle")
parLapply(clus, 1:n, validFn) 

stopCluster(clus)
