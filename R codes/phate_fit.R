#!/usr/bin/env Rscript
library('FactoMineR')
library("stringr")
library('mapproj')
library('maptools')
library('RColorBrewer')
library('parallel')
library('phateR')
library('scatterplot3d')

source('axis_map.R')
cambria <- readRDS('cambria.rds')
type = commandArgs(trailingOnly = T)[1]
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
if (type=='prob'){
  raw_data <-readRDS('predictions_new_stations_model-mean_1deg.rds')
  new_stations_1 <- readRDS('new_stations_1deg.rds')
} else if (type=='dom'){
  raw_data06 <-readRDS('pred_dom06.rds')
  raw_data90 <-readRDS('pred_dom90.rds')
  raw_data <- rbind(raw_data06, raw_data90)
  rm(raw_data06,raw_data90)
  mapping <- readRDS('mapping_lo_lt.rds')
  new_stations_1 <- readRDS('new_stations_1deg.rds')
  new_stations_1 <- new_stations_1[mapping,]
} else if (type=='all'){
  raw_data0690 <-readRDS('predictions_new_stations_model-mean_1deg.rds')
  seq1 <- seq(1, dim(raw_data0690)[1], 2)
  seq2 <- seq(2, dim(raw_data0690)[1], 2)
  raw_data06 <- raw_data0690[seq1,]
  raw_data90 <- raw_data0690[seq2,]
  raw_data_woa <- readRDS('model-mean_pred_woa_list.rds')
  raw_data <- rbind(raw_data_woa,raw_data06, raw_data90)
  new_stations_1 <- readRDS('new_stations_1deg.rds')
} else if (type=='06_prob'){
  raw_data0690 <-readRDS('predictions_new_stations_model-mean_1deg.rds')
  seq1 <- seq(1, dim(raw_data0690)[1], 2)
  raw_data06 <- raw_data0690[seq1,]
  raw_data <- raw_data06
  new_stations_1 <- readRDS('new_stations_1deg.rds')
} else if (type=='90_prob'){
  raw_data0690 <-readRDS('predictions_new_stations_model-mean_1deg.rds')
  seq1 <- seq(2, dim(raw_data0690)[1], 2)
  raw_data90 <- raw_data0690[seq1,]
  raw_data <- raw_data90
  new_stations_1 <- readRDS('new_stations_1deg.rds')
}



if (type %in% c('prob', 'dom', 'discrete')){
  phate_fit <- phate(raw_data, ndim = 3, knn = 2000, decay = 20)
  saveRDS(phate_fit, paste('phate_fit_0690_',type,'.rds', sep=''))
} else if (type=='all'){
  phate_fit <- phate(raw_data, ndim = 3, knn = 3000, decay = 20)
  saveRDS(phate_fit, paste('phate_fit_0690woa_',type,'.rds', sep=''))
} else if (type %in% c('90_prob', '06_prob')){
  phate_fit <- phate(raw_data, ndim = 3, knn = 1000, decay = 20)
  saveRDS(phate_fit, paste('phate_fit_',type,'.rds', sep=''))
}

x <- phate_fit$embedding[,1]
y <- phate_fit$embedding[,2]
z <- phate_fit$embedding[,3]

x_n <- (x - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255
y_n <- (y - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255
z_n <- (z - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255

pdf(family="Helvetica",paste('rgb_phate_1deg_',type,'.pdf', sep=''), width=10,height=10)
plot(x_n, y_n, col=rgb(x_n, y_n, z_n, maxColorValue = 255), xlab = '1st axis', ylab='2nd axis')
dev.off()

colors <-rgb(x_n, y_n, z_n, maxColorValue = 255)

if (type=='prob'){
  seq1=seq(1,dim(raw_data)[1], 2)
  seq2=seq(2,dim(raw_data)[1], 2)
} else if (type=='dom'){
  seq1=1:(dim(raw_data)[1]/2)
  seq2=((dim(raw_data)[1]/2)+1):dim(raw_data)[1]
} else if (type=='all'){
  N=dim(raw_data)[1]/3
  seq1=1:N
  seq2=(N+1):(2*N)
  seq3=(2*N+1):(3*N)
}
if (type %in% c('prob', 'dom', 'discrete')){
  pdf(family="Helvetica",file=paste('rgb_2006_phate_1deg_',type,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(25,10,25,10))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in seq1){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  #axis_map()
  dev.off()
  
  pdf(family="Helvetica",file=paste('rgb_2090_phate_1deg_',type,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(25,10,25,10))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in seq2){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  #axis_map()
  dev.off()
} else if (type=='all'){
  pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(25,10,25,10))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in seq1){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  #axis_map()
  dev.off()
  
  pdf(family="Helvetica",file=paste('rgb_2006_phate_1deg_',type,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(25,10,25,10))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in seq2){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  #axis_map()
  dev.off()
  
  pdf(family="Helvetica",file=paste('rgb_2090_phate_1deg_',type,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(25,10,25,10))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in seq3){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  #axis_map()
  dev.off()
} else if (type %in% c('06_prob', '90_prob')){
  pdf(family="Helvetica",file=paste('rgb_phate_1deg_',type,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(25,10,25,10))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in 1:dim(raw_data)[1]){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  #axis_map()
  dev.off()
}

pdf(family="Helvetica",paste('3D_phate_',type,'.pdf', sep=''), width = 10, height = 10)
scatterplot3d(x = x, y = y, z= z, color = rgb(x_n, y_n, z_n, maxColorValue = 255), 
              xlab = "1st dim", ylab = "2nd dim", zlab = "3rd dim" )
dev.off()

