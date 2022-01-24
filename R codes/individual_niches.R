#!bin/usr/bin/env Rscript
# projections_cor.R
# Code to project individual environmental niches of each 'metacommunity' of each size fraction on 3 oceanic biogeochemical Earth System Models (IPSL-CM5A-LR,
# GFDL-ESM2G and MPI-ESM-LR) and a mean of 7 models (IPSL-CM5A-LR/MR, GFDL-ESM2G/M and MPI-ESM-LR/MR and CESM1-BGC) in 2006-15 decade and 2090-2099 decade. 
# Environmental niches models are fitted using World Ocean Atlas 2013 (WOA13) extrapolated observations + WOA01 for chlorophyll A and 4 machine learning algorithms
# WOA13 data is also used to correct outputs of ESMs in 2006-15 and 2090-99 using Cumulative Distribution Function transform method (CDFt) 
library("gbm")
library("dismo")
library("FactoMineR")
#library("readxl")
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('maptools')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('CDFt')
library('randomForest')
library('mgcv')
library('nnet')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
source('axis_map0.R')
source('hide_arctic.R')
cambria <-readRDS('cambria.rds')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
# Function to calculate distance between 2 earth points
getDistance <-function(lat1,lon1,lat2,lon2){
  R = 6371; # Radius of the earth in km
  dLat = deg2rad(lat2-lat1)  # deg2rad below
  dLon = deg2rad(lon2-lon1) 
  a = sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2) * sin(dLon/2)
  c = 2 * atan2(sqrt(a), sqrt(1-a)); 
  d = R * c; # Distance in km
  getDistance <- d
}

# Converting degree to radian
deg2rad <- function(deg) {
  deg2ead <- deg * (pi/180)
}

vec2grid <- function(data, j){
  test <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(data)[1]){
    test[which(longi==data$Long[i]),which(lati==data$Lat[i])]= round(data[i,j])+1
  }
  return(test)
}
set = colorRampPalette(c('blue','yellow', 'red'))(101)
set0 = colorRampPalette(c('blue','yellow', 'red'))(99)
set1 = colorRampPalette(c("darkgreen","light blue", "dark violet"))(201)
set2 = colorRampPalette(c("darkgreen","light blue", "dark violet"))(21)
best_models <- read.table('best_selected_models.txt', header = T)
best_models <- as.data.frame(best_models)
N <- dim(best_models)[1]
row.names(best_models)<-c(1:N)
#models <- c('IPSL-CM5A-LR', 'MPI-ESM-LR','model-mean', 'GFDL-ESM2G')
models_bt <- readRDS('models_bt.rds')
models_gam <- readRDS('models_gam.rds')

lati<-seq(-89.5, 60, 1)
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
longi <- sort(longi)

model='model-mean'
fractions <- c('180-2000', '20-180', '43952','0.8-5','0.22-3', '0-0.2')
distances1 <- NULL
cors_all_woa <- NULL
uncertainties_all_woa <- list(NULL,NULL,NULL)
cors_all_2006 <- NULL
uncertainties_all_2006 <- list(NULL,NULL,NULL)
cors_all_2090 <- NULL
uncertainties_all_2090 <- list(NULL,NULL,NULL)
pred_woa_list <- readRDS(paste(model,'_pred_woa_list.rds', sep=''))
pred_2006_list <- readRDS(paste(model,'_pred_2006_list.rds', sep=''))
pred_2090_list <- readRDS(paste(model,'_pred_2090_list.rds', sep=''))
pred_woa_list_delta_sdm <- readRDS(paste(model,'_pred_woa_list_delta_sdm.rds', sep=''))
pred_2006_list_delta_sdm <- readRDS(paste(model,'_pred_2006_list_delta_sdm.rds', sep=''))
pred_2090_list_delta_sdm <- readRDS(paste(model,'_pred_2090_list_delta_sdm.rds', sep=''))
pred_woa_list_sd_sdm <- readRDS(paste(model,'_pred_woa_list_sd_sdm.rds', sep=''))
pred_2006_list_sd_sdm <- readRDS(paste(model,'_pred_2006_list_sd_sdm.rds', sep=''))
pred_2090_list_sd_sdm <- readRDS(paste(model,'_pred_2090_list_sd_sdm.rds', sep=''))
preds_all_nc <- readRDS('predictions_new_stations_model-mean_1deg_nc.rds')
pred_2006_list_nc <- preds_all_nc[seq(1, dim(preds_all_nc)[1], 2),]
pred_2090_list_nc <- preds_all_nc[seq(2, dim(preds_all_nc)[1], 2),]
colnames(pred_2006_list_nc) <- colnames(pred_2006_list)
colnames(pred_2090_list_nc) <- colnames(pred_2006_list)
stats <- readRDS('new_stations_1deg.rds')
colnames(stats)<- c('Lat', 'Long')
stats <- as.data.frame(stats)
df <- readRDS('Genocenoses_env_parameters_woa.rds')
frac <- '180-2000'
clusts <- c(5, 8)
# for(clust in clusts){
#   df1 <-df[df$Fraction==frac,]
#   df1$Genocenose <- as.numeric(df1$Genocenose==clust)+1
#   colors <- c('red','green' )
#   pdf(family="Helvetica",paste(frac,'_', clust,'.pdf', sep=''),width=10,height=4.065)
#   par(mar=c(0,0,0,0))
#   maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
#   plot(coastline,lwd=0.0475, col='black', add=T)
#   points(df1$Long, df1$Lat, col=colors[df1$Genocenose], pch=19, cex=20)
#   axis_map()
#   dev.off()
# }
count=1
for (count in 1:N){
  frac <- best_models$Fraction[count]
  k <- best_models$Gen[count]
  bt_model<-models_bt[[count]]
  gam_model<-models_gam[[count]]
  preds_all_woa <- NULL
  preds_all_06 <- NULL
  preds_all_90 <- NULL
  if (frac=='43952'){
    frac_ok <- '5-20'
  } else{
    frac_ok <- frac
  }
  index <- which(colnames(pred_2006_list)==paste(frac_ok, k, sep='_'))
  pred_bt_woa <- readRDS(paste(model,'_pred_woa_list_bt.rds', sep=''))[, index]
  pred_bt_06 <- readRDS(paste(model,'_pred_2006_list_bt.rds', sep=''))[, index]
  pred_bt_90 <- readRDS(paste(model,'_pred_2090_list_bt.rds', sep=''))[, index]
  preds_all_woa <- append(preds_all_woa, pred_bt_woa)
  preds_all_06 <- append(preds_all_06, pred_bt_06)
  preds_all_90 <- append(preds_all_90, pred_bt_90)
  pred_rf_woa <- readRDS(paste(model,'_pred_woa_list_rf.rds', sep=''))[, index]
  pred_rf_06 <- readRDS(paste(model,'_pred_2006_list_rf.rds', sep=''))[, index]
  pred_rf_90 <- readRDS(paste(model,'_pred_2090_list_rf.rds', sep=''))[, index]
  preds_all_woa <- append(preds_all_woa, pred_rf_woa)
  preds_all_06 <- append(preds_all_06, pred_rf_06)
  preds_all_90 <- append(preds_all_90, pred_rf_90)
  pred_nn_woa <- readRDS(paste(model,'_pred_woa_list_nn.rds', sep=''))[, index]
  pred_nn_06 <- readRDS(paste(model,'_pred_2006_list_nn.rds', sep=''))[, index]
  pred_nn_90 <- readRDS(paste(model,'_pred_2090_list_nn.rds', sep=''))[, index]
  preds_all_woa <- append(preds_all_woa, pred_nn_woa)
  preds_all_06 <- append(preds_all_06, pred_nn_06)
  preds_all_90 <- append(preds_all_90, pred_nn_90)
  pred_gam_woa<- readRDS(paste(model,'_pred_woa_list_gam.rds', sep=''))[, index]
  pred_gam_06 <- readRDS(paste(model,'_pred_2006_list_gam.rds', sep=''))[, index]
  pred_gam_90 <- readRDS(paste(model,'_pred_2090_list_gam.rds', sep=''))[, index]
  preds_all_woa <- append(preds_all_woa, pred_gam_woa)
  preds_all_06 <- append(preds_all_06, pred_gam_06)
  preds_all_90 <- append(preds_all_90, pred_gam_90)
  
  
  bc<- readRDS('model-mean_bray_curtis.rds')
  
  ndw <- data.frame(Long=bc$Long, Lat=bc$Lat, Pr=pred_woa_list[,index]*100)
  nd <- data.frame(Long=bc$Long, Lat=bc$Lat, Pr=pred_2006_list[,index]*100)
  nd1 <- data.frame(Long=bc$Long, Lat=bc$Lat, Pr=pred_2090_list[,index]*100)
  nd_nc <- data.frame(Long=stats$Long, Lat=stats$Lat, Pr=pred_2006_list_nc[,index]*100)
  nd1_nc <- data.frame(Long=stats$Long, Lat=stats$Lat, Pr=pred_2090_list_nc[,index]*100)
  nd0 <- data.frame(Long=bc$Long, Lat=bc$Lat,Pr=(nd1$Pr-nd$Pr)+100)
  nd3 <- data.frame(Long=bc$Long, Lat=bc$Lat,
                    d_woa=100*pred_woa_list_delta_sdm[,index],
                    d_2006=100*pred_2006_list_delta_sdm[,index],
                    d_2090=100*pred_2090_list_delta_sdm[,index],
                    sd_woa=100*pred_woa_list_sd_sdm[,index],
                    sd_2006=100*pred_2006_list_sd_sdm[,index],
                    sd_2090=100*pred_2090_list_sd_sdm[,index])
  
  
  
  # map of the delta probability of presence of the metacommunity 'k' of size fraction 'frac' between 2006-15 and 2090-99
  data_contour <- vec2grid(nd0,3)
  data_contour1 <- round(data_contour, -1)/10+1
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_delta_2090-2006.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_delta_2090-2006.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set2, levels = c(1:21))
  hide_arctic()
  axis_map0()
  dev.off()
  # map of the probabilty of presence of metacommunity 'k' of size fraction 'frac' using woa2013 data
  data_contour <- vec2grid(ndw,3)
  data_contour1 <- round(data_contour, -1)/10+1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_woa.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_woa.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  #     pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))    
  #     points(x=ndw$Long, y=ndw$Lat, col=set[round(ndw$Pr*100)+1],bg= set[round(ndw$Pr*100)+1],pch=15,cex=2.83)
  #     SDMTools::legend.gradient(pnt,set,limits=c(0,1),title='Probability of presence',cex=0.4)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set0, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  
  # map of the probabilty of presence of metacommunity 'k' of size fraction 'frac' in 2006-15 then 2090-99
  data_contour <- vec2grid(nd,3)
  data_contour1 <- round(data_contour, -1)/10 +1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_2006-15.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_2006-15.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  #     pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))    
  #     points(x=nd$Long, y=nd$Lat, col=set[round(nd$Pr*100)+1],bg= set[round(nd$Pr*100)+1],pch=15,cex=2.83)
  #     SDMTools::legend.gradient(pnt,set,limits=c(0,1),title='Probability of presence',cex=0.4)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set0, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- vec2grid(nd1,3)
  data_contour1 <- round(data_contour,-1)/10+1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_2090-99.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_2090-99.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  #     pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
  #     points(x=nd1$Long, y=nd1$Lat, col=set[round(nd1$Pr*100)+1],bg= set[round(nd1$Pr*100)+1],pch=15,cex=2.83)
  #     SDMTools::legend.gradient(pnt,set,limits=c(0,1),title='Probability of presence',cex=0.4)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set0, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- vec2grid(nd_nc,3)
  data_contour1 <- round(data_contour, -1)/10 +1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_2006-15_nc.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_2006-15_nc.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set0, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- vec2grid(nd1_nc,3)
  data_contour1 <- round(data_contour,-1)/10+1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_2090-99_nc.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_2090-99_nc.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set0, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  preds_all_mat_woa <- matrix(preds_all_woa, ncol=4)
  preds_all_mat_06 <- matrix(preds_all_06, ncol=4)
  preds_all_mat_90 <- matrix(preds_all_90, ncol=4)
  for (i in 1:3){
    for (j in (i+1):4){
      cors_all_woa <-append(cors_all_woa, cor(preds_all_mat_woa[,i], preds_all_mat_woa[,j]))
      cors_all_2006 <-append(cors_all_2006, cor(preds_all_mat_06[,i], preds_all_mat_06[,j]))
      cors_all_2090 <-append(cors_all_2090, cor(preds_all_mat_90[,i], preds_all_mat_90[,j]))
    }
  }
  
  colors_d <-  colorRampPalette(c('blue','green', 'black'))(100)
  colors_sd <-  colorRampPalette(c('blue','green', 'black'))(101)
  
  data_contour <- vec2grid(nd3, 3)
  data_contour1 <- round(data_contour, -1)/10+1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_delta_sdm_model_woa.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_delta_sdm_model_woa.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  #     points(x=nd3$Long, y=nd3$Lat, col=colors_d[nd3$d_woa+1],bg=colors_d[nd3$d_woa+1],pch=15,cex=2.83)
  #     SDMTools::legend.gradient(pnt,colors_d,limits=c(0,1),title='Maximum delta probability\n of presence',cex=0.4)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=colors_d, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- vec2grid(nd3, 4)
  data_contour1 <- round(data_contour, -1)/10+1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_delta_sdm_model_2006.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_delta_sdm_model_2006.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  #     points(x=nd3$Long, y=nd3$Lat, col=colors_d[nd3$d_2006+1],bg=colors_d[nd3$d_2006+1],pch=15,cex=2.83)
  #     SDMTools::legend.gradient(pnt,colors_d,limits=c(0,1),title='Maximum delta probability\n of presence',cex=0.4)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=colors_d, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- vec2grid(nd3, 5)
  data_contour1 <- round(data_contour, -1)/10+1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_delta_sdm_model_2090.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_delta_sdm_model_2090.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  #     points(x=nd3$Long, y=nd3$Lat, col=colors_d[nd3$d_2090+1],bg=colors_d[nd3$d_2090+1],pch=15,cex=2.83)
  #     SDMTools::legend.gradient(pnt,colors_d,limits=c(0,1),title='Maximum delta probability\n of presence',cex=0.4)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=colors_d, levels = c(1:NewMax))
  hide_arctic()
  axis_map0()
  dev.off()
  
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_sd_sdm_model_woa.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_sd_sdm_model_woa.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  indexes <- 100*(nd3$sd_woa-min(nd3$sd_woa))/(max(nd3$sd_woa)-min(nd3$sd_woa))+1
  points(x=nd1$Long, y=nd1$Lat, col=colors_sd[indexes],bg=colors_sd[indexes],pch=15,cex=2.83)
  SDMTools::legend.gradient(pnt,colors_d,limits=c(round(min(nd3$sd_woa), digits=1),round(max(nd3$sd_woa), digits=1)),
                            title='Standard deviation',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_sd_sdm_model_2006.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_sd_sdm_model_2006.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  indexes <- 100*(nd3$sd_2006-min(nd3$sd_2006))/(max(nd3$sd_2006)-min(nd3$sd_2006))+1
  points(x=nd1$Long, y=nd1$Lat, col=colors_sd[indexes],bg=colors_sd[indexes],pch=15,cex=2.83)
  SDMTools::legend.gradient(pnt,colors_d,limits=c(round(min(nd3$sd_2006), digits=1),round(max(nd3$sd_2006), digits=1)),
                            title='Standard deviation',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_sd_sdm_model_2090.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_sd_sdm_model_2090.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  indexes <- 100*(nd3$sd_2090-min(nd3$sd_2090))/(max(nd3$sd_2090)-min(nd3$sd_2090))+1
  points(x=nd1$Long, y=nd1$Lat, col=colors_sd[indexes],bg=colors_sd[indexes],pch=15,cex=2.83)
  SDMTools::legend.gradient(pnt,colors_d,limits=c(round(min(nd3$sd_2090), digits=1),round(max(nd3$sd_2090), digits=1)),
                            title='Standard deviation',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  date<-c('woa','2006','2090')
  for (i in 1:3){
    if (i ==1){
      diffs<-nd3$d_woa
      pr<-ndw$Pr
    } else if (i==2){
      diffs<-nd3$d_2006
      pr<-nd$Pr
    } else if (i==3){
      diffs<-nd3$d_2090
      pr<-nd1$Pr
    }
    uncertainties <- list()
    if (length( diffs[pr>66])>0 ){
      uncertainties[[1]] <- diffs[pr>66]
    } else {
      uncertainties[[1]] <- 1
    }
    if (length(diffs[pr<66 & pr>33]) > 0){
      uncertainties[[2]] <- diffs[pr<66 & pr>33]
    } else{
      uncertainties[[2]] <- 1
    }
    if (length( diffs[pr<33])>0 ){
      uncertainties[[3]] <- diffs[pr<33]
    } else{
      uncertainties[[3]] <- 1
    }
    for (cu in 1:3){
      if ( length(uncertainties[[cu]]) >1 ){
        if (i==1){
          uncertainties_all_woa[[cu]] <- append(uncertainties_all_woa[[cu]], uncertainties[[cu]])
        } else if (i==2){
          uncertainties_all_2006[[cu]] <- append(uncertainties_all_2006[[cu]], uncertainties[[cu]])
        } else if (i==3){
          uncertainties_all_2090[[cu]] <- append(uncertainties_all_2090[[cu]], uncertainties[[cu]])
        }
      }
    }
    
    
    if (frac != '43952'){
      pdf(family="Helvetica",file=paste(model,'_',frac,'_metacommunity_',as.character(k),'_delta_sdm_model_distribution_',date[i],'.pdf', sep=''),,width=20,height=10)
    } else{
      pdf(family="Helvetica",file=paste(model,'_5-20_metacommunity_',as.character(k),'_delta_sdm_model_distribution_',date[i],'.pdf', sep=''),width=20,height=10)
    }
    colors_ryb <- c('red', 'yellow', 'blue')
    if ( length(uncertainties[[1]]) == 1 & length(uncertainties[[2]]) > 1 & length(uncertainties[[3]]) > 1  ){
      plot(c(5,5), main='', , xlim=c(0,100), ylim=c(0,max(density(uncertainties[[2]])$y, density(uncertainties[[3]])$y)),
           xlab='Maximum delta probability of presence',ylab='Density', col='white',xaxs="i" )
    } else if (length(uncertainties[[1]]) > 1 & length(uncertainties[[2]]) == 1 & length(uncertainties[[3]]) > 1 ){
      plot(c(5,5), main='', , xlim=c(0,100), ylim=c(0,max(density(uncertainties[[1]])$y, density(uncertainties[[3]])$y)),
           xlab='Maximum delta probability of presence',ylab='Density', col='white',xaxs="i" )
    } else if (length(uncertainties[[1]]) > 1 & length(uncertainties[[2]]) > 1 & length(uncertainties[[3]]) == 1 ){
      plot(c(5,5), main='', , xlim=c(0,100), ylim=c(0,max(density(uncertainties[[1]])$y, density(uncertainties[[2]])$y)),
           xlab='Maximum delta probability of presence',ylab='Density', col='white',xaxs="i" )
    } else {
      plot(c(5,5), main='', xlim=c(0,100), ylim=c(0,max(density(uncertainties[[1]] , cut=7)$y,density(uncertainties[[2]], cut=7)$y, density(uncertainties[[3]], cut=7)$y )),
           xlab='Maximum delta probability of presence',ylab='Density', col='white',xaxs="i" )
    }
    for (count in 1:3){
      if (length(uncertainties[[count]]) >1){
        d <- density(uncertainties[[count]], cut=7)
        d$y[d$x < 0] <- 0
        d$y[d$x >100] <- 0
        polygon(d, col=scales::alpha(colors_ryb[count], 0.5), border=scales::alpha(colors_ryb[count], 0.5))
        legend('topright', legend = c('P>0.66', '0.33<P<0.66', 'P<0.33'), col=colors_ryb,pch=15,cex=2, bty='n')
      }
    }
    dev.off()
  }
  
  latitude_w <- nd$Lat*cos(nd$Lat*2*pi/360)
  latitude_w1 <- nd1$Lat*cos(nd1$Lat*2*pi/360)
  
  lat_w <- cos(nd$Lat*2*pi/360)
  lat_w1 <- cos(nd1$Lat*2*pi/360)
  
  long0 <- nd$Long
  long1 <- nd$Long
  long0[nd$Long< -80]<- 360-abs(nd$Long[nd$Long< -80])
  long1[nd$Long< -70]<- 360-abs(nd$Long[nd$Long< -70])
  
  # Calculating the centroid of the metacommunity 'k' of size fraction 'frac' in 5 major oceanic basins in 2006-15 and 2090-99
  cond_na <- nd$Long<30 & nd$Long > -80 & nd$Lat >0 & nd$Pr>50
  cond_na1 <- nd1$Long<30 & nd1$Long > -80 & nd1$Lat >0 & nd1$Pr>50
  
  mean_lat_na = sum(latitude_w[cond_na]*nd$Pr[cond_na]/100)/sum(nd$Pr[cond_na]*lat_w[cond_na]/100)
  mean_long_na = sum(nd$Long[cond_na]*nd$Pr[cond_na]/100)/sum(nd$Pr[cond_na]/100)
  mean_lat_na1 = sum(latitude_w1[cond_na1]*nd1$Pr[cond_na1]/100)/sum(nd1$Pr[cond_na1]*lat_w1[cond_na1]/100)
  mean_long_na1 = sum(nd1$Long[cond_na1]*nd1$Pr[cond_na1]/100)/sum(nd1$Pr[cond_na1]/100)
  
  cond_sa <- nd$Long<20 & nd$Long > -70 & nd$Lat < 0 & nd$Pr>50
  cond_sa1 <- nd1$Long<20 & nd1$Long > -70 & nd1$Lat <0 & nd1$Pr>50
  
  mean_lat_sa = sum(latitude_w[cond_sa]*nd$Pr[cond_sa]/100)/sum(nd$Pr[cond_sa]*lat_w[cond_sa]/100)
  mean_long_sa = sum(nd$Long[cond_sa]*nd$Pr[cond_sa]/100)/sum(nd$Pr[cond_sa]/100)
  mean_lat_sa1 = sum(latitude_w1[cond_sa1]*nd1$Pr[cond_sa1]/100)/sum(nd1$Pr[cond_sa1]*lat_w1[cond_sa1]/100)
  mean_long_sa1 = sum(nd1$Long[cond_sa1]*nd1$Pr[cond_sa1]/100)/sum(nd1$Pr[cond_sa1]/100)
  
  cond_np <- (nd$Long>120 | nd$Long < -90) & nd$Lat >0 & nd$Pr>50
  cond_np1 <- (nd1$Long>120 | nd1$Long < -90) & nd1$Lat >0 & nd1$Pr>50
  
  mean_lat_np = sum(latitude_w[cond_np]*nd$Pr[cond_np]/100)/sum(nd$Pr[cond_np]*lat_w[cond_np]/100)
  mean_long_np = sum(long0[cond_np]*nd$Pr[cond_np]/100)/sum(nd$Pr[cond_np]/100)
  if (mean_long_np>180 & !(is.nan(mean_long_np))){
    mean_long_np = -180+(mean_long_np-180)
  }
  mean_lat_np1 = sum(latitude_w1[cond_np1]*nd1$Pr[cond_np1]/100)/sum(nd1$Pr[cond_np1]*lat_w1[cond_np1]/100)
  mean_long_np1 = sum(long0[cond_np1]*nd1$Pr[cond_np1]/100)/sum(nd1$Pr[cond_np1]/100)
  if (mean_long_np1>180 & !(is.nan(mean_long_np1))){
    mean_long_np1 = -180+(mean_long_np1-180)
  }
  
  cond_sp <- (nd$Long>120 | nd$Long < -70) & nd$Lat <0 & nd$Pr>50
  cond_sp1 <- (nd1$Long>120 | nd1$Long < -70) & nd1$Lat <0 & nd1$Pr>50
  
  mean_lat_sp = sum(latitude_w[cond_sp]*nd$Pr[cond_sp]/100)/sum(nd$Pr[cond_sp]*lat_w[cond_sp]/100)
  mean_long_sp = sum(long1[cond_sp]*nd$Pr[cond_sp]/100)/sum(nd$Pr[cond_sp]/100)
  if (mean_long_sp>180 & !(is.nan(mean_long_sp))){
    mean_long_sp = -180+(mean_long_sp-180)
  }
  mean_lat_sp1 = sum(latitude_w1[cond_sp1]*nd1$Pr[cond_sp1]/100)/sum(nd1$Pr[cond_sp1]*lat_w1[cond_sp1]/100)
  mean_long_sp1 = sum(long1[cond_sp1]*nd1$Pr[cond_sp1]/100)/sum(nd1$Pr[cond_sp1]/100)
  if (mean_long_sp1>180 & !(is.nan(mean_long_sp1))){
    mean_long_sp1 = -180+(mean_long_sp1-180)
  }
  
  cond_i <- nd$Long>20 & nd$Long < 120 & nd$Lat <30 & nd$Pr>50
  cond_i1 <- nd1$Long>20 & nd1$Long < 120 & nd1$Lat <30 & nd1$Pr>50
  
  mean_lat_i = sum(latitude_w[cond_i]*nd$Pr[cond_i]/100)/sum(nd$Pr[cond_i]*lat_w[cond_i]/100)
  mean_long_i = sum(nd$Long[cond_i]*nd$Pr[cond_i]/100)/sum(nd$Pr[cond_i]/100)
  mean_lat_i1 = sum(latitude_w1[cond_i1]*nd1$Pr[cond_i1]/100)/sum(nd1$Pr[cond_i1]*lat_w1[cond_i1]/100)
  mean_long_i1 = sum(nd1$Long[cond_i1]*nd1$Pr[cond_i1]/100)/sum(nd1$Pr[cond_i1]/100)
  
  longs <- c(mean_long_na, mean_long_sa, mean_long_np, mean_long_sp, mean_long_i)
  longs1 <- c(mean_long_na1, mean_long_sa1, mean_long_np1, mean_long_sp1, mean_long_i1)
  lats <- c(mean_lat_na, mean_lat_sa, mean_lat_np, mean_lat_sp, mean_lat_i)
  lats1 <- c(mean_lat_na1, mean_lat_sa1, mean_lat_np1, mean_lat_sp1, mean_lat_i1)
  colors <- c('red', 'cyan', 'orange', 'deeppink', 'blue')
  
  locs <- c('North_A', 'South_A', 'North_P', 'South_P', 'Indian') # The 5 major basins
  distances0 <- NULL
  for (i in 1:5){
    if ( !(is.na(lats[i])) & !(is.na(lats1[i])) & !(is.na(longs[i])) & !(is.na(longs1[i]))){
      dis <- getDistance(lats[i], longs[i], lats1[i], longs1[i])
      dis_lat <- getDistance(lats[i], longs[i], lats1[i], longs[i])
      dis_long <- getDistance(lats[i], longs[i], lats[i], longs1[i])
      if (lats1[i] < lats[i]){
        dis_lat <- -dis_lat # southtward or north movement
      }
      if (longs[i]>0 & longs1[i]>0 & longs[i]>longs1[i]){
        dis_long <- -dis_long ## westward or eastward movement
      } else if (longs[i]>0 & longs1[i]<0){
        dis_long <- -dis_long ## westward or eastward movement
      } else if (longs1[i]<0 & longs[i]<0 & longs1[i]<longs[i]){
        dis_long <- -dis_long ## westward or eastward movement
      } 
      distances0 <- rbind(distances0, c(frac, k,model, locs[i],lats[i],longs[i],lats1[i],longs1[i], dis, dis_lat, dis_long))
    }
  }
  distances1 <- rbind(distances1, distances0)
  
  # Delta map with the centroids
  data_contour=vec2grid(nd0,3)
  data_contour1 = round(data_contour, -1)/10+1
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_shifts_delta_2090-2006.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_shifts_delta_2090-2006.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set2, levels = c(1:21))
  # points(x=nd0$Long, y=nd0$Lat, col=set1[nd0$Pr],bg= set1[nd0$Pr],pch=15,cex=2.83)
  points(longs, lats, col = colors, pch =19, cex = 1.2)
  #text(longs, lats, labels = rep('2006', 5), cex = 3)
  points(longs1, lats1, col = colors, pch =17, cex = 1.2)
  #text(longs1, lats1, labels = rep('2090', 5), cex = 3)
  #SDMTools::legend.gradient(pnt,set1,limits=c(-1,1),title='Delta of \nprobability of presence',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  if (frac != '43952'){
    pdf(family="Helvetica",file=paste('projection_',model,'_',frac,'_metacommunity_',as.character(k),'_shifts_delta_2090-2006_arrows.pdf', sep=''),width=10,height=4.065)
  } else{
    pdf(family="Helvetica",file=paste('projection_',model,'_5-20_metacommunity_',as.character(k),'_shifts_delta_2090-2006_arrows.pdf', sep=''),width=10,height=4.065)
  }
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(x=longi,y=lati, z=data_contour1,  col=set0, levels = c(1:21))
  arrows(x0 = longs, y0 = lats, x1 = longs1, y1=lats1, col=colors, lwd=1.2, lenght=0.025)
  #SDMTools::legend.gradient(pnt,set1,limits=c(-1,1),title='Delta of \nprobability of presence',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
}
pdf(family="Helvetica",file=paste(model,'_delta_sdm_model_distribution_woa.pdf', sep=''),width=20,height=10)
colors_ryb <- c('red', 'yellow', 'blue')
plot(c(10,10), main='',  xlim=c(0,100), 
     ylim=c(0,max(density(uncertainties_all_woa[[1]], cut=7)$y, density(uncertainties_all_woa[[2]], cut=7)$y, density(uncertainties_all_woa[[3]], cut=7)$y)), 
     xlab='Maximum delta probability of presence',ylab='Density', col='white',xaxs="i")
for (count in 1:3){
  d <- density(uncertainties_all_woa[[count]], cut=7)
  d$y[d$x < 0] <- 0
  d$y[d$x >100] <- 0
  polygon(d, col=scales::alpha(colors_ryb[count], 0.5), border=scales::alpha(colors_ryb[count], 0.5))
  legend('topright', legend = c('P>0.66', '0.33<P<0.66', 'P<0.33'), col=colors_ryb,pch=15,cex=2, bty='n')
}
dev.off()

pdf(family="Helvetica",file=paste(model,'_delta_sdm_model_distribution_2006.pdf', sep=''),width=20,height=10)
colors_ryb <- c('red', 'yellow', 'blue')
plot(c(10,10), main='',  xlim=c(0,100), 
     ylim=c(0,max(density(uncertainties_all_2006[[1]], cut=7)$y, density(uncertainties_all_2006[[2]], cut=7)$y, density(uncertainties_all_2006[[3]], cut=7)$y)), 
     xlab='Maximum delta probability of presence',ylab='Density', col='white',xaxs="i")
for (count in 1:3){
  d <- density(uncertainties_all_2006[[count]], cut=7)
  d$y[d$x < 0] <- 0
  d$y[d$x >100] <- 0
  polygon(d, col=scales::alpha(colors_ryb[count], 0.5), border=scales::alpha(colors_ryb[count], 0.5))
  legend('topright', legend = c('P>0.66', '0.33<P<0.66', 'P<0.33'), col=colors_ryb,pch=15,cex=2, bty='n')
}
dev.off()

pdf(family="Helvetica",file=paste(model,'_delta_sdm_model_distribution_2090.pdf', sep=''),width=20,height=10)
colors_ryb <- c('red', 'yellow', 'blue')
plot(c(10,10), main='',  xlim=c(0,100), 
     ylim=c(0,max(density(uncertainties_all_2090[[1]], cut=7)$y, density(uncertainties_all_2090[[2]], cut=7)$y, density(uncertainties_all_2090[[3]], cut=7)$y)), 
     xlab='Maximum delta probability of presence',ylab='Density', col='white',xaxs="i")
for (count in 1:3){
  d <- density(uncertainties_all_2090[[count]], cut=7)
  d$y[d$x < 0] <- 0
  d$y[d$x >100] <- 0
  polygon(d, col=scales::alpha(colors_ryb[count], 0.5), border=scales::alpha(colors_ryb[count], 0.5))
  legend('topright', legend = c('P>0.66', '0.33<P<0.66', 'P<0.33'), col=colors_ryb,pch=15,cex=2, bty='n')
}
dev.off()

#colnames(distances1) <- c("Fraction" ,"Metacommunity" , "Model","Location","Lat_2006", 'Long_2006',"Lat_2090", 'Long_2090',"Shift" ,"Lat_shift" ,"Long_shift")
write.table(distances1, file=paste('Distances_shifts_',model,'.txt', sep=''), row.names = F, col.names = F)

cors_mat_woa <- t(matrix(cors_all_woa, ncol=N))
write.table(cors_mat_woa, 'correlation_sdm_models_esm_woa.txt')
cors_mat_2006 <- t(matrix(cors_all_2006, ncol=N))
write.table(cors_mat_2006, 'correlation_sdm_models_esm_2006.txt')
cors_mat_2090 <- t(matrix(cors_all_2090, ncol=N))
write.table(cors_mat_2090, 'correlation_sdm_models_esm_2090.txt')

labs <- NULL
labs0 <- NULL
code_fraction <-c('F', 'E','D','C','B','A')
count=1
for (u in best_models$Fraction){
  gen <- best_models$Gen[count]
  if (u != "43952"){
    labs <-append(labs, paste(u, gen, sep='_'))
    letter <- code_fraction[match(u, fractions)]
    labs0 <- append(labs0, paste(letter,gen, sep=''))
  } else{
    labs <- append(labs, paste('5-20', gen, sep='_'))
    letter <- code_fraction[match(u, fractions)]
    labs0 <- append(labs0, paste(letter,gen, sep=''))
  }
  count=count+1
}

dates <- c('woa', '2006','2090')
for (date in dates){
  br <- c(seq(-0.1, 1, length.out = 16))
  cor_table <- read.table(paste('correlation_sdm_models_esm_',date,'.txt', sep=''))
  data <- as.matrix(cor_table)
  pdf(family="Helvetica",paste('correlation_sdm_models_',model,'_',date,'.pdf', sep=''), height= 7, width=10)
  heatmap.2(data, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
            dendrogram = "none", keysize=1,margins=c(10,22), 
            labCol=c('rf vs gbm', 'gam vs gbm','nn vs gbm','gam vs rf',  'gam vs nn', 'nn vs rf'),
            labRow=labs, col= greenred(15), breaks=br, symkey = F, density.info = 'none')
  dev.off()
  pdf(family="Helvetica",paste('correlation_sdm_models_',model,'_',date,'_1.pdf', sep=''), height= 7, width=10)
  heatmap.2(data, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
            dendrogram = "none", keysize=1,margins=c(10,22), 
            labCol=c('rf vs gbm', 'gam vs gbm','nn vs gbm','gam vs rf',  'gam vs nn', 'nn vs rf'),
            labRow=labs0, col= greenred(15), breaks=br, symkey = F, density.info = 'none', cexRow=1.2)
  dev.off()
}

