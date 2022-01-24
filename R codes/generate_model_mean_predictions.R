#!bin/usr/bin/env Rscript
library("gbm")
library('randomForest')
library('mgcv')
library('nnet')
library("dismo")
library("FactoMineR")
#library("readxl")
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library("CDFt")
library('parallel')

df <- readRDS('Genocenoses_env_parameters_woa_scaled.rds')

best_models <- read.table('best_selected_models.txt', header = T)
best_models <- as.data.frame(best_models)
N <-dim(best_models)[1]
row.names(best_models)<-c(1:N)
variables <- c(6:11,13)
variables1 <- c('thetao','so', 'si','no3','po4', 'dfe')
fractions <- c('180-2000', '20-180', '43952', '0.8-5', '0.22-3', '0-0.2')

model0<-'model-mean'
black_caspien <- expand.grid(seq(28.5, 60.5, 1),seq(37.5, 48.5, 1))
med_sea <- expand.grid(seq(-2.5, 35.5, 1),seq(30.5, 45.5, 1))
closed_sea <- rbind(black_caspien, med_sea)
variables <- c(6:11,13)
variables1 <- c('thetao','so', 'si','no3','po4', 'dfe')

data_w <- readRDS('data_w_scaled.rds')
data_p <- readRDS('data_p_scaled.rds')
data_f <- readRDS('data_f_scaled.rds')
data <- readRDS('model-mean_2006-15_scaled.rds')
data1 <- readRDS('model-mean_2090-99_scaled.rds')

models_bt <- readRDS('models_bt.rds')
models_gam <- readRDS('models_gam.rds')
models_nn <- readRDS('models_nn.rds')
models_rf <- readRDS('models_rf.rds')
good <- readRDS('good_models.rds')

start <- Sys.time()
long<-seq(0.5,359.5,1)
long[181:360] =as.numeric(long[181:360])-360
lat<-seq(-89.5,89.5,1)
new_stations <- expand.grid(lat[seq(1,150,1)], long[seq(1,360,1)])
new_stations_1 <- NULL
new_stations_1_nc <- NULL
predictors_new_stations_w <- NULL
predictors_new_stations <- NULL
predictors_new_stations_nc <- NULL
predictors_noT <- NULL
predictors_f_drivers_al <- rep(list(NULL), length(variables))
for (st in 1:dim(new_stations)[1]){
  lt = new_stations[st,1]
  lg = new_stations[st,2]
  predictors_w <- NULL
  predictors <- NULL
  predictors_f <- NULL
  predictors_nc <- NULL
  predictors_f_nc <- NULL
  predictors_f_noT <- NULL
  for (l in 1:length(variables)){
    lgi = which(long==lg)
    lti = which(lat==lt)
    predictors_w <- append(predictors_w, data_w[[l]][lgi,lti])
    predictors <- append(predictors, data_p[[l]][lgi,lti])
    predictors_f <- append(predictors_f, data_f[[l]][lgi,lti])
    if (l!=1){
      predictors_f_noT <- append(predictors_f_noT, data_f[[l]][lgi,lti])
    } else{
      predictors_f_noT <- append(predictors_f_noT, data_p[[l]][lgi,lti])
    }
    predictors_nc <- append(predictors_nc, data[[l]][lgi,lti])
    predictors_f_nc <- append(predictors_f_nc, data1[[l]][lgi,lti])
  }
  
  predictors_drivers_al <- rep(list(NULL), length(variables))
  for (c in 1:length(variables)){
    pred_d_al <- predictors
    pred_d_al[c] <-  predictors_f[c]
    predictors_drivers_al[[c]] <- append(predictors_drivers_al[[c]], pred_d_al)
  }
  
  condi <- !is.na(sum(predictors_w)) & !is.null(predictors_w) & !is.na(sum(predictors)) & !is.null(predictors) & !is.na(sum(predictors_f)) & !is.null(predictors_f) & !(lg %in% closed_sea[[1]] & lt %in% closed_sea[[2]]) 
  if (condi){
    new_stations_1 <- rbind(new_stations_1, c(lt, lg))
    predictors_new_stations_w <- append(predictors_new_stations_w, predictors_w)
    predictors_new_stations <- append(predictors_new_stations, c(predictors,predictors_f))
    predictors_noT <- append(predictors_noT, predictors_f_noT)
    predictors_new_stations_nc <- append(predictors_new_stations_nc, c(predictors_nc,predictors_f_nc))
  }
  for (j in 1:length(variables)){
    if (!is.na(sum(predictors_drivers_al[[j]])) & !is.null(predictors_drivers_al[[j]]) & !(lg %in% closed_sea[[1]] & lt %in% closed_sea[[2]])){
      predictors_f_drivers_al[[j]] <- append(predictors_f_drivers_al[[j]], predictors_drivers_al[[j]])
    }
  }
}
end <- Sys.time()
duration <- end-start

predictors_new_stations_w <- matrix(predictors_new_stations_w, ncol=length(variables), byrow=T)
predictors_new_stations_w <- as.data.frame(predictors_new_stations_w)
colnames(predictors_new_stations_w)<-names(df[0, variables])
predictors_new_stations <- matrix(predictors_new_stations, ncol=length(variables), byrow=T)
predictors_new_stations <- as.data.frame(predictors_new_stations)
colnames(predictors_new_stations)<-names(df[0, variables])
predictors_new_stations_nc <- matrix(predictors_new_stations_nc, ncol=length(variables), byrow=T)
predictors_new_stations_nc <- as.data.frame(predictors_new_stations_nc)
colnames(predictors_new_stations_nc)<-names(df[0, variables])
predictors_noT<- matrix(predictors_noT, ncol=length(variables), byrow=T)
predictors_noT<- as.data.frame(predictors_noT)
colnames(predictors_noT)<-names(df[0, variables])
for (h in 1:length(variables)){
  predictors_f_drivers_al[[h]] <- matrix(predictors_f_drivers_al[[h]], ncol=length(variables), byrow=T)
  predictors_f_drivers_al[[h]] <- as.data.frame(predictors_f_drivers_al[[h]])
  colnames(predictors_f_drivers_al[[h]])<-names(df[0, variables])
}

u=1
predictions_new_stations_w <- NULL
predictions_new_stations <- NULL
predictions_new_stations_nc <- NULL
predictions_new_stations_noT <- NULL
pred_woa_list_bt <-NULL
pred_2006_list_bt <- NULL
pred_2090_list_bt <- NULL
pred_woa_list_rf<-NULL
pred_2006_list_rf <- NULL
pred_2090_list_rf <- NULL
pred_woa_list_nn<-NULL
pred_2006_list_nn <- NULL
pred_2090_list_nn <- NULL
pred_woa_list_gam<-NULL
pred_2006_list_gam <- NULL
pred_2090_list_gam <- NULL
pred_woa_list_delta_sdm <-NULL
pred_2006_list_delta_sdm <- NULL
pred_2090_list_delta_sdm <- NULL
pred_woa_list_sd_sdm <-NULL
pred_2006_list_sd_sdm <- NULL
pred_2090_list_sd_sdm <- NULL
drivers_al0 <- rep(list(NULL), length(variables))
for (m in good){
  bt_model <- models_bt[[u]]
  rf_model <- models_rf[[u]]
  nn_model <- models_nn[[u]]
  gam_model <- models_gam[[u]]
  
  pred_rf_w <- stats::predict(rf_model, predictors_new_stations_w, type='prob')[,2]
  pred_nn_w <- stats::predict(nn_model, predictors_new_stations_w, type='raw')[,2]
  pred_rf <- stats::predict(rf_model, predictors_new_stations, type='prob')[,2]
  pred_nn <- stats::predict(nn_model, predictors_new_stations, type='raw')[,2]
  pred_rf_nc <- stats::predict(rf_model, predictors_new_stations_nc, type='prob')[,2]
  pred_nn_nc <- stats::predict(nn_model, predictors_new_stations_nc, type='raw')[,2]
  pred_rf_noT <- stats::predict(rf_model, predictors_noT, type='prob')[,2]
  pred_nn_noT <- stats::predict(nn_model, predictors_noT, type='raw')[,2]
  len <- length(pred_rf_w)
  if (!is.na(bt_model)){
    pred_bt_w <- gbm::predict.gbm(bt_model,predictors_new_stations_w, n.trees=bt_model$gbm.call$best.trees, type="response")
    pred_bt <- gbm::predict.gbm(bt_model,predictors_new_stations, n.trees=bt_model$gbm.call$best.trees, type="response")
    pred_bt_nc <- gbm::predict.gbm(bt_model,predictors_new_stations_nc, n.trees=bt_model$gbm.call$best.trees, type="response")
    pred_bt_noT <- gbm::predict.gbm(bt_model,predictors_noT, n.trees=bt_model$gbm.call$best.trees, type="response")
    
    pred_woa_list_bt <-append(pred_woa_list_bt, pred_bt_w)
    pred_2006_list_bt <- append(pred_2006_list_bt, pred_bt[seq(1, length(pred_bt), 2)])
    pred_2090_list_bt <- append(pred_2090_list_bt, pred_bt[seq(2, length(pred_bt), 2)])
  } else{
    pred_woa_list_bt <-append(pred_woa_list_bt,rep(NA, len))
    pred_2006_list_bt <- append(pred_2006_list_bt, rep(NA, len))
    pred_2090_list_bt <- append(pred_2090_list_bt, rep(NA, len))
  }
  if (!is.na(gam_model)){
    pred_gam_w <- mgcv::predict.gam(gam_model, predictors_new_stations_w, type='response')
    pred_gam <- mgcv::predict.gam(gam_model, predictors_new_stations, type='response')
    pred_gam_nc <- mgcv::predict.gam(gam_model, predictors_new_stations_nc, type='response')
    pred_gam_noT <- mgcv::predict.gam(gam_model, predictors_noT, type='response')
    
    pred_woa_list_gam<-append(pred_woa_list_gam,pred_gam_w)
    pred_2006_list_gam <- append(pred_2006_list_gam,pred_gam[seq(1, length(pred_gam), 2)])
    pred_2090_list_gam <- append(pred_2090_list_gam , pred_gam[seq(2, length(pred_gam), 2)])
  } else{
    pred_woa_list_gam<-append(pred_woa_list_gam,rep(NA, len))
    pred_2006_list_gam <- append(pred_2006_list_gam,rep(NA, len))
    pred_2090_list_gam <- append(pred_2090_list_gam , rep(NA, len))
  }
  
  
  cond_06 = seq(1, length(pred_rf), 2)
  cond_90 = seq(2, length(pred_rf), 2)
  pred_woa_list_rf<-append(pred_woa_list_rf,pred_rf_w)
  pred_2006_list_rf <- append(pred_2006_list_rf, pred_rf[cond_06])
  pred_2090_list_rf <- append(pred_2090_list_rf,pred_rf[cond_90])
  
  pred_woa_list_nn<-append(pred_woa_list_nn,pred_nn_w)
  pred_2006_list_nn <- append(pred_2006_list_nn,pred_nn[cond_06])
  pred_2090_list_nn <- append(pred_2090_list_nn, pred_nn[cond_90])
  
  
  if (!is.na(bt_model) & !is.na(gam_model)){
    all_pred_w <- data.frame('bt'=pred_bt_w, 'nn'=pred_nn_w, 'rf'=pred_rf_w, 'gam'=pred_gam_w)
    all_pred_06 <- data.frame('bt'=pred_bt[cond_06], 'nn'=pred_nn[cond_06],'rf'=pred_rf[cond_06], 'gam'=pred_gam[cond_06])
    all_pred_90 <- data.frame('bt'=pred_bt[cond_90], 'nn'=pred_nn[cond_90],'rf'=pred_rf[cond_90], 'gam'=pred_gam[cond_90])
    pred_w <- (pred_gam_w+pred_rf_w+pred_bt_w+pred_nn_w)/4
    pred <- (pred_gam+pred_rf+pred_bt+pred_nn)/4
    pred_nc <- (pred_gam_nc+pred_rf_nc+pred_bt_nc+pred_nn_nc)/4
    pred_noT <- (pred_gam_noT+pred_rf_noT+pred_bt_noT+pred_nn_noT)/4
  } else if (is.na(bt_model) & !is.na(gam_model)){
    all_pred_w <- data.frame( 'nn'=pred_nn_w, 'rf'=pred_rf_w, 'gam'=pred_gam_w)
    all_pred_06 <- data.frame( 'nn'=pred_nn[cond_06],'rf'=pred_rf[cond_06], 'gam'=pred_gam[cond_06])
    all_pred_90 <- data.frame( 'nn'=pred_nn[cond_90],'rf'=pred_rf[cond_90], 'gam'=pred_gam[cond_90])
    pred_w <- (pred_gam_w+pred_rf_w+pred_nn_w)/3
    pred <- (pred_gam+pred_rf+pred_nn)/3
    pred_nc <- (pred_gam_nc+pred_rf_nc+pred_nn_nc)/3
    pred_noT <- (pred_gam_noT+pred_rf_noT+pred_nn_noT)/3
  } else if (!is.na(bt_model) & is.na(gam_model)){
    all_pred_w <- data.frame( 'nn'=pred_nn_w, 'rf'=pred_rf_w, 'bt'=pred_gam_w)
    all_pred_06 <- data.frame( 'nn'=pred_nn[cond_06],'rf'=pred_rf[cond_06], 'bt'=pred_gam[cond_06])
    all_pred_90 <- data.frame( 'nn'=pred_nn[cond_90],'rf'=pred_rf[cond_90], 'bt'=pred_gam[cond_90])
    pred_w <- (pred_bt_w+pred_rf_w+pred_nn_w)/3
    pred <- (pred_bt+pred_rf+pred_nn)/3
    pred_nc <- (pred_bt_nc+pred_rf_nc+pred_nn_nc)/3
    pred_noT <- (pred_bt_noT+pred_rf_noT+pred_nn_noT)/3
  }
  
  delta_max <- function(u){
    return(max(u)-min(u))
  }
  
  sd_w <- apply(all_pred_w, 1, sd)
  sd_06 <- apply(all_pred_06, 1, sd)
  sd_90 <- apply(all_pred_90, 1, sd)
  d_w <- apply(all_pred_w, 1, delta_max)
  d_06 <- apply(all_pred_06, 1, delta_max)
  d_90 <- apply(all_pred_90, 1, delta_max)
  
  pred_woa_list_delta_sdm <-append(pred_woa_list_delta_sdm , d_w)
  pred_2006_list_delta_sdm <- append(pred_2006_list_delta_sdm , d_06)
  pred_2090_list_delta_sdm <- append(pred_2090_list_delta_sdm , d_90)
  pred_woa_list_sd_sdm <-append(pred_woa_list_sd_sdm ,sd_w)
  pred_2006_list_sd_sdm <- append(pred_2006_list_sd_sdm ,sd_06)
  pred_2090_list_sd_sdm <- append(pred_2090_list_sd_sdm , sd_90)
  
  predictions_new_stations_w <- append(predictions_new_stations_w, pred_w)
  predictions_new_stations <- append(predictions_new_stations, pred)
  predictions_new_stations_nc <- append(predictions_new_stations_nc, pred_nc)
  predictions_new_stations_noT <- append(predictions_new_stations_noT, pred_noT)
  
  for (h in 1:length(variables)){
    prdctors <- predictors_f_drivers_al[[h]]
    if (!is.na(bt_model)){
      pred_bt <- gbm::predict.gbm(bt_model,prdctors, n.trees=bt_model$gbm.call$best.trees, type="response")
    }
    pred_rf <- stats::predict(rf_model, prdctors, type='prob')[,2]
    pred_nn <- stats::predict(nn_model, prdctors, type='raw')[,2]
    if (!is.na(gam_model)){
      pred_gam <- mgcv::predict.gam(gam_model, prdctors, type='response')
    }
    if (!is.na(bt_model) & !is.na(gam_model)){
      pred <- (pred_gam+pred_rf+pred_bt+pred_nn)/4
    } else if (is.na(bt_model) & !is.na(gam_model)){
      pred <- (pred_gam+pred_rf+pred_nn)/3
    } else if (!is.na(bt_model) & is.na(gam_model)){
      pred <- (pred_bt+pred_rf+pred_nn)/3
    }
    drivers_al0[[h]] <- append(drivers_al0[[h]], pred) 
  }
  u=u+1
}
predictions_new_stations_w <- matrix(predictions_new_stations_w, ncol=length(good))
saveRDS(predictions_new_stations_w, 'predictions_new_stations_w_model-mean_1deg.rds')
predictions_new_stations <- matrix(predictions_new_stations, ncol=length(good))
saveRDS(predictions_new_stations, 'predictions_new_stations_model-mean_1deg.rds')
saveRDS(new_stations_1, 'new_stations_1deg.rds')

predictions_new_stations_nc <- matrix(predictions_new_stations_nc, ncol=length(good))
saveRDS(predictions_new_stations_nc, 'predictions_new_stations_model-mean_1deg_nc.rds')
saveRDS(new_stations_1_nc, 'new_stations_1deg_nc.rds')

predictions_new_stations_noT <- matrix(predictions_new_stations_noT, ncol=length(good))
saveRDS(predictions_new_stations_noT, 'predictions_new_stations_noT_model-mean_1deg.rds')

for (h in 1:length(variables)){
  drivers_al0[[h]]<- matrix(drivers_al0[[h]], ncol=length(good))
}

count=1
labs <- NULL
for (u in best_models$Fraction){
  gen <- best_models$Gen[count]
  if (u != "43952"){
    labs <-append(labs, paste(u, gen, sep='_'))
  } else{
    labs <- append(labs, paste('5-20', gen, sep='_'))
  }
  count=count+1
}

pred_woa_list <- predictions_new_stations_w
pred_2006_list <- predictions_new_stations[seq(1, dim(predictions_new_stations)[1], 2),]
pred_2090_list <- predictions_new_stations[seq(2, dim(predictions_new_stations)[1], 2),]
pred_2090_list_noT <- predictions_new_stations_noT


d06_90 <- abs(pred_2006_list-pred_2090_list)
s06_90 <- pred_2006_list+pred_2090_list
diff06_90 <- apply(d06_90,1,sum)
sum06_90 <- apply(s06_90,1,sum)
bc <-1-diff06_90/sum06_90
#

bray_curtis <- new_stations_1
bray_curtis <- cbind(new_stations_1, bc)
colnames(bray_curtis) <- c('Lat', 'Long', 'b_c')
rownames(bray_curtis) <- NULL
bray_curtis <- as.data.frame(bray_curtis)
bray_curtis$cell<-as.character(paste(bray_curtis$Lat,bray_curtis$Long, sep='_'))

sum_driv_al <- rep(0, dim(pred_2006_list)[1])
for (e in 1:length(variables)){
  x <- abs(drivers_al0[[e]]-pred_2006_list) 
  y <- y <- apply(x, 1, sum)
  sum_driv_al <- sum_driv_al+y 
}
drivers_al <- rep(list(NULL), length(variables))
for (e in 1:length(variables)){
  drivers_al[[e]] <- new_stations_1
  x <- abs(drivers_al0[[e]]-pred_2006_list)
  y <- apply(x, 1, sum)
  drivers_al[[e]] <- cbind(drivers_al[[e]], y/sum_driv_al)
  drivers_al[[e]] <- cbind(drivers_al[[e]], 1-bray_curtis$b_c)
}

pred_woa_list_sd_sdm <- matrix(pred_woa_list_sd_sdm, ncol=length(good))
pred_2006_list_sd_sdm <- matrix(pred_2006_list_sd_sdm, ncol=length(good))
pred_2090_list_sd_sdm <- matrix(pred_2090_list_sd_sdm, ncol=length(good))
pred_woa_list_delta_sdm <- matrix(pred_woa_list_delta_sdm, ncol=length(good))
pred_2006_list_delta_sdm <- matrix(pred_2006_list_delta_sdm, ncol=length(good))
pred_2090_list_delta_sdm <- matrix(pred_2090_list_delta_sdm, ncol=length(good))
pred_woa_list_gam <- matrix(pred_woa_list_gam, ncol=length(good))
pred_woa_list_nn <- matrix(pred_woa_list_nn, ncol=length(good))
pred_woa_list_rf <- matrix(pred_woa_list_rf, ncol=length(good))
pred_woa_list_bt <- matrix(pred_woa_list_bt, ncol=length(good))
pred_2006_list_gam <- matrix(pred_2006_list_gam, ncol=length(good))
pred_2006_list_nn <- matrix(pred_2006_list_nn, ncol=length(good))
pred_2006_list_rf <- matrix(pred_2006_list_rf, ncol=length(good))
pred_2006_list_bt <- matrix(pred_2006_list_bt, ncol=length(good))
pred_2090_list_gam <- matrix(pred_2090_list_gam, ncol=length(good))
pred_2090_list_nn <- matrix(pred_2090_list_nn, ncol=length(good))
pred_2090_list_rf <- matrix(pred_2090_list_rf, ncol=length(good))
pred_2090_list_bt <- matrix(pred_2090_list_bt, ncol=length(good))

colnames(pred_woa_list)<-labs
colnames(pred_2006_list)<-labs
colnames(pred_2090_list)<-labs
colnames(pred_2090_list_noT)<-labs
row.names(pred_woa_list)<-NULL
row.names(pred_2006_list)<-NULL
row.names(pred_2090_list)<-NULL
row.names(pred_2090_list_noT)<-NULL
colnames(pred_woa_list_delta_sdm)<-labs
colnames(pred_2006_list_delta_sdm)<-labs
colnames(pred_2090_list_delta_sdm)<-labs
row.names(pred_woa_list_delta_sdm)<-NULL
row.names(pred_2006_list_delta_sdm)<-NULL
row.names(pred_2090_list_delta_sdm)<-NULL
colnames(pred_woa_list_sd_sdm)<-labs
colnames(pred_2006_list_sd_sdm)<-labs
colnames(pred_2090_list_sd_sdm)<-labs
row.names(pred_woa_list_sd_sdm)<-NULL
row.names(pred_2006_list_sd_sdm)<-NULL
row.names(pred_2090_list_sd_sdm)<-NULL
colnames(pred_woa_list_gam)<-labs
colnames(pred_2006_list_gam)<-labs
colnames(pred_2090_list_gam)<-labs
row.names(pred_woa_list_gam)<-NULL
row.names(pred_2006_list_gam)<-NULL
row.names(pred_2090_list_gam)<-NULL
colnames(pred_woa_list_rf)<-labs
colnames(pred_2006_list_rf)<-labs
colnames(pred_2090_list_rf)<-labs
row.names(pred_woa_list_rf)<-NULL
row.names(pred_2006_list_rf)<-NULL
row.names(pred_2090_list_rf)<-NULL
colnames(pred_woa_list_nn)<-labs
colnames(pred_2006_list_nn)<-labs
colnames(pred_2090_list_nn)<-labs
row.names(pred_woa_list_nn)<-NULL
row.names(pred_2006_list_nn)<-NULL
row.names(pred_2090_list_nn)<-NULL
colnames(pred_woa_list_bt)<-labs
colnames(pred_2006_list_bt)<-labs
colnames(pred_2090_list_bt)<-labs
row.names(pred_woa_list_bt)<-NULL
row.names(pred_2006_list_bt)<-NULL
row.names(pred_2090_list_bt)<-NULL

saveRDS(bray_curtis, paste(model0, '_bray_curtis.rds', sep=''))
saveRDS(pred_woa_list, paste(model0, '_pred_woa_list.rds', sep=''))
saveRDS(pred_2006_list, paste(model0, '_pred_2006_list.rds', sep=''))
saveRDS(pred_2090_list, paste(model0, '_pred_2090_list.rds', sep=''))
saveRDS(pred_2090_list_noT, paste(model0, '_pred_2090_list_noT.rds', sep=''))
saveRDS(pred_woa_list_delta_sdm, paste(model0, '_pred_woa_list_delta_sdm.rds', sep=''))
saveRDS(pred_2006_list_delta_sdm, paste(model0, '_pred_2006_list_delta_sdm.rds', sep=''))
saveRDS(pred_2090_list_delta_sdm, paste(model0, '_pred_2090_list_delta_sdm.rds', sep=''))
saveRDS(pred_woa_list_sd_sdm, paste(model0, '_pred_woa_list_sd_sdm.rds', sep=''))
saveRDS(pred_2006_list_sd_sdm, paste(model0, '_pred_2006_list_sd_sdm.rds', sep=''))
saveRDS(pred_2090_list_sd_sdm, paste(model0, '_pred_2090_list_sd_sdm.rds', sep=''))
saveRDS(pred_woa_list_gam, paste(model0, '_pred_woa_list_gam.rds', sep=''))
saveRDS(pred_2006_list_gam, paste(model0, '_pred_2006_list_gam.rds', sep=''))
saveRDS(pred_2090_list_gam, paste(model0, '_pred_2090_list_gam.rds', sep=''))
saveRDS(pred_woa_list_rf, paste(model0, '_pred_woa_list_rf.rds', sep=''))
saveRDS(pred_2006_list_rf, paste(model0, '_pred_2006_list_rf.rds', sep=''))
saveRDS(pred_2090_list_rf, paste(model0, '_pred_2090_list_rf.rds', sep=''))
saveRDS(pred_woa_list_nn, paste(model0, '_pred_woa_list_nn.rds', sep=''))
saveRDS(pred_2006_list_nn, paste(model0, '_pred_2006_list_nn.rds', sep=''))
saveRDS(pred_2090_list_nn, paste(model0, '_pred_2090_list_nn.rds', sep=''))
saveRDS(pred_woa_list_bt, paste(model0, '_pred_woa_list_bt.rds', sep=''))
saveRDS(pred_2006_list_bt, paste(model0, '_pred_2006_list_bt.rds', sep=''))
saveRDS(pred_2090_list_bt, paste(model0, '_pred_2090_list_bt.rds', sep=''))
saveRDS(drivers_al,paste(model0, '_drivers_al.rds', sep='') )
weight_cos <- cos(new_stations_1[,1]*pi/180)
saveRDS(weight_cos, paste(model0, '_weights_cos.rds', sep=''))

