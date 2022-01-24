#!/usr/bin/env Rscript
library('FactoMineR')
library("stringr")
library('mapproj')
library('maptools')
library('RColorBrewer')
library('parallel')
library('phateR')
library('scatterplot3d')
library('ncdf4')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
source('axis_map0.R')
source('hide_arctic.R')
cambria <- readRDS('cambria.rds')
type = commandArgs(trailingOnly = T)[1]
mod= commandArgs(trailingOnly = T)[2]
n_opt = commandArgs(trailingOnly = T)[3]

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
}

test <- ncdf4::nc_open('Time_Varying_Biomes.nc')
biome <- ncvar_get(test, 'MeanBiomes')
biome[is.nan(biome)]<-0
lt <- ncvar_get(test, 'lat')
lo <- ncvar_get(test, 'lon')

dimred_fit <- readRDS(paste(mod,'_fit_0690_',type,'.rds', sep=''))
clust_medoids <-readRDS(paste(mod,'_cluster_medoids_',type,'_',n_opt,'_0690.rds', sep=''))

clusts_medoids_mat_06 <- matrix(0, ncol = 360, nrow = 180)
clusts_medoids_mat_90 <- matrix(0, ncol = 360, nrow = 180)
clusts_medoids_mat_0690 <- matrix(NA, ncol = 360, nrow = 180)
N <- dim(new_stations_1)[1]
for (j in 1:N){
  lati <- new_stations_1[j,1]
  loni <- new_stations_1[j,2]
  ind_la <- which(lt==lati)
  ind_lo <- which(lo==loni)
  clusts_medoids_mat_06[ind_la, ind_lo]<-clust_medoids[j*2-1]
  clusts_medoids_mat_90[ind_la, ind_lo]<-clust_medoids[j*2]
  clusts_medoids_mat_0690[ind_la, ind_lo]<-as.numeric(!clust_medoids[j*2-1]==clust_medoids[j*2])
}

if (mod=='tsne'){
  data <- dimred_fit$Y
} else if (mod=='phate'){
  data <- dimred_fit$embedding
}

x <- data[,1]
y <- data[,2]
z <- data[,3]

x_n <- (x - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255
y_n <- (y - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255
z_n <- (z - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255

# pdf(family="Helvetica",paste('rgb_phate_1deg_',type,'.pdf', sep=''), width=10,height=10)
# plot(x_n, y_n, col=rgb(x_n, y_n, z_n, maxColorValue = 255), xlab = '1st axis', ylab='2nd axis')
# dev.off()

colors <-rgb(x_n, y_n, z_n, maxColorValue = 255)

if (type=='prob'){
  seq1=seq(1,dim(raw_data)[1], 2)
  seq2=seq(2,dim(raw_data)[1], 2)
} else if (type=='dom'){
  seq1=1:(dim(raw_data)[1]/2)
  seq2=((dim(raw_data)[1]/2)+1):dim(raw_data)[1]
}

pdf(family="Helvetica",file=paste('rgb_2006_',mod,'_1deg_',type,'_', n_opt,'_clusts.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
count = 1
for(j in seq1){
  c = c(new_stations_1[count,1], new_stations_1[count,2])
  points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
  count =count+1
}
for (q in unique(clusts_medoids_mat_06[clusts_medoids_mat_06 !=0])){
  hi<-t(clusts_medoids_mat_06)
  hi[t(clusts_medoids_mat_06)!=q]=0
  contour(lo, lt, hi, col = 'white',add=T, drawlabels = F, cex=0.284)
}
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",file=paste('rgb_2090_',mod,'_1deg_',type,'_', n_opt,'_clusts.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
count = 1
for(j in seq2){
  c = c(new_stations_1[count,1], new_stations_1[count,2])
  points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
  count =count+1
}
for (q in unique(clusts_medoids_mat_90[clusts_medoids_mat_90 !=0])){
  hi<-t(clusts_medoids_mat_90)
  hi[t(clusts_medoids_mat_90)!=q]=0
  contour(lo, lt, hi, col = 'white',add=T, drawlabels = F, cex=0.284)
}
hide_arctic()
axis_map0()
dev.off()


pdf(family="Helvetica",file=paste('rgb_2090-2006_',mod,'_1deg_',type,'_', n_opt,'_clusts.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
.filled.contour(lo, lt,t(clusts_medoids_mat_0690),levels = c(0,1,2),
                col = c('blue','red',scales::alpha('white', 0)) )
hide_arctic()
axis_map0()
dev.off()
# pdf(family="Helvetica",paste('3D_phate_',type,'.pdf', sep=''), width = 10, height = 10)
# scatterplot3d(x = x, y = y, z= z, color = rgb(x_n, y_n, z_n, maxColorValue = 255), 
#               xlab = "1st dim", ylab = "2nd dim", zlab = "3rd dim" )
# dev.off()
