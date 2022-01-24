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
source('hide_arctic0.R')
type = commandArgs(trailingOnly = T)[1]
mod = commandArgs(trailingOnly = T)[2]
opt_n =  commandArgs(trailingOnly = T)[3]

cambria <- readRDS('cambria.rds')
clust_medoids <-readRDS(paste(mod,'_cluster_medoids_',type,'_',opt_n,'_woa.rds', sep=''))
raw_data <- readRDS('pred_domWOA.rds')
raw_data[raw_data!=0]<-1
raw_data1 <- raw_data
raw_data1 <- unique(raw_data1)
new_stations_1 <- readRDS('new_stations_1deg.rds')

if (type=='discrete'){
  func_com <- function(u){
    com <- u[1]
    for (i in u[2:length(u)]){
      com <- paste(com, i, sep='_')
    }
    return(com)
  }
  
  dt =as.data.frame(raw_data)
  fused_raw_data <- apply(dt, 1, func_com)
  fused_raw_data1 <- apply(data.frame(raw_data1),1, func_com)
  clust_medoids_1 <- NULL
  for (i in 1:dim(raw_data)[1]){
    cl <- clust_medoids[which(fused_raw_data[i]==fused_raw_data1)]
    clust_medoids_1 <- append(clust_medoids_1, cl)
  }
  clust_medoids <- clust_medoids_1
  saveRDS(clust_medoids, paste(mod,'_cluster_medoids_',type,'_',opt_n,'_all.rds', sep=''))
}


test <- ncdf4::nc_open('Time_Varying_Biomes.nc')
biome <- ncvar_get(test, 'MeanBiomes')
biome[is.nan(biome)]<-0
lt <- ncvar_get(test, 'lat')
lo <- ncvar_get(test, 'lon')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
antarctic_fronts <- readShapeSpatial('ant_fronts/antarctic_circumpolar_current_fronts.shp')
BGCP_shp <- readShapeSpatial('longhurst_v4_2010/Longhurst_world_v4_2010.shp')
reyg_bgcp <- as.matrix(read.csv2('Reygondeau/PROVINCE_REYGONDEAU.csv', header = F, sep = ','))
reyg_bio <- as.matrix(read.csv2('Reygondeau/BIOME_REYGONDEAU.csv', header = F, sep = ','))
reyg_bgcp19_1 <-as.matrix(readRDS('BGCP_2019_REYGONDEAU_mat_1.rds'))
#reyg_bgcp[is.nan(reyg_bgcp)]<-0
reyg_bio[is.nan(reyg_bio)]<-0
reyg_bgcp19_1[is.nan(reyg_bgcp19_1)]<-0

clusts_medoids_mat <- matrix(0, ncol = 360, nrow = 180)
for (j in 1:length(clust_medoids)){
  lati <- new_stations_1[j,1]
  loni <- new_stations_1[j,2]
  ind_la <- which(lt==lati)
  ind_lo <- which(lo==loni)
  clusts_medoids_mat[ind_la, ind_lo]<-clust_medoids[j]
}
to_plot <-list(biome,reyg_bgcp19_1, reyg_bgcp, reyg_bio)
names <- c('McKingley','Reygondeau_BGCP19','Reygondeau_BGCP','Reygondeau_BIOME')

raw_data <- readRDS('pred_domWOA.rds')
new_stations_1 <- readRDS('new_stations_1deg.rds')

#phate_fit <- phate(raw_data, ndim = 3, knn = 2500, decay = 40)
#saveRDS(phate_fit, 'phate_fit_woa.rds')

phate_fit <- readRDS(paste(mod,'_fit_woa_',type,'.rds', sep=''))
if (mod=='phate'){
  data <- phate_fit$embedding
} else if (mod=='tsne'){
  data <- phate_fit$Y
}

x <- data[,1]
y <- data[,2]
z <- data[,3]

x_n <- (x - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255
y_n <- (y - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255
z_n <- (z - min(x,y,z))/(max(x,y,z) - min(x,y,z))*255

pdf(family="Helvetica",paste('rgb_phate_1deg_',type,'_',mod,'_',opt_n,'_woa.pdf', sep=''), width=10,height=10)
plot(x_n, y_n, col=rgb(x_n, y_n, z_n, maxColorValue = 255), xlab = '1st axis', ylab='2nd axis')
dev.off()

colors <-rgb(x_n, y_n, z_n, maxColorValue = 255)

if (type=='discrete'){
  pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'_',mod,'_',opt_n,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in 1:dim(raw_data)[1]){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    co = colors[which(fused_raw_data[j]==fused_raw_data1)]
    points(x=c[2], y=c[1],  col=co,bg= co,pch=15,cex=0.235)
    count =count+1
  }
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'_',mod,'_',opt_n,'_clusts.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in 1:dim(raw_data)[1]){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    co = colors[which(fused_raw_data[j]==fused_raw_data1)]
    points(x=c[2], y=c[1],  col=co,bg= co,pch=15,cex=0.235)
    count =count+1
  }
  for (q in unique(clusts_medoids_mat[clusts_medoids_mat !=0])){
    hi<-t(clusts_medoids_mat)
    hi[t(clusts_medoids_mat)!=q]=0
    contour(lo, lt, hi, col = 'white',add=T, drawlabels = F, cex=0.284)
  }
  hide_arctic()
  axis_map0()
  dev.off()
  
  for (cou in 1:4){
    pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'_',mod,'_',opt_n,'_',names[cou],'.pdf', sep=''),width=10,height=4.065)
    par(mar=c(0,0,0,0))
    maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
    plot(coastline,lwd=0.0475, col='black', add=T)
    count = 1
    for(j in 1:dim(raw_data)[1]){
      c = c(new_stations_1[count,1], new_stations_1[count,2])
      co = colors[which(fused_raw_data[j]==fused_raw_data1)]
      points(x=c[2], y=c[1],  col=co,bg= co,pch=15,cex=0.235)
      count =count+1
    }
    for (q in unique(clusts_medoids_mat[clusts_medoids_mat !=0])){
      hi<-t(clusts_medoids_mat)
      hi[t(clusts_medoids_mat)!=q]=0
      contour(lo, lt, hi, col = 'white',add=T, drawlabels = F, cex=0.284)
    }
    biogeog <- t(to_plot[[cou]])
    for (q in unique(biogeog[biogeog !=0])){
      hi<-biogeog
      hi[biogeog!=q]=0
      contour(lo, lt, hi, col = 'black',add=T, drawlabels = F, cex=0.284)
    }
    hide_arctic()
    axis_map0()
    dev.off()
  }
} else{
  pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'_',mod,'_',opt_n,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in 1:dim(raw_data)[1]){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'_',mod,'_',opt_n,'_clusts.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in 1:dim(raw_data)[1]){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  for (q in unique(clusts_medoids_mat[clusts_medoids_mat !=0])){
    hi<-t(clusts_medoids_mat)
    hi[t(clusts_medoids_mat)!=q]=0
    contour(lo, lt, hi, col = 'white',add=T, drawlabels = F, cex=0.284)
  }
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'_',mod,'_',opt_n,'_clusts_antarctic_fronts.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  count = 1
  for(j in 1:dim(raw_data)[1]){
    c = c(new_stations_1[count,1], new_stations_1[count,2])
    points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
    count =count+1
  }
  for (q in unique(clusts_medoids_mat[clusts_medoids_mat !=0])){
    hi<-t(clusts_medoids_mat)
    hi[t(clusts_medoids_mat)!=q]=0
    contour(lo, lt, hi, col = 'white',add=T, drawlabels = F, cex=0.284)
  }
  plot(antarctic_fronts, col='black',lwd=2.5, add=T)
  hide_arctic()
  axis_map0()
  dev.off()
  
  for (cou in 1:4){
    pdf(family="Helvetica",file=paste('rgb_woa_phate_1deg_',type,'_',mod,'_',opt_n,'_',names[cou],'.pdf', sep=''),width=10,height=4.065)
    par(mar=c(0,0,0,0))
    maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
    plot(coastline,lwd=0.0475, col='black', add=T)
    count = 1
    for(j in 1:dim(raw_data)[1]){
      c = c(new_stations_1[count,1], new_stations_1[count,2])
      points(x=c[2], y=c[1],  col=colors[j],bg= colors[j],pch=15,cex=0.235)
      count =count+1
    }
    for (q in unique(clusts_medoids_mat[clusts_medoids_mat !=0])){
      hi<-t(clusts_medoids_mat)
      hi[t(clusts_medoids_mat)!=q]=0
      contour(lo, lt, hi, col = 'white',add=T, drawlabels = F, cex=0.284)
    }
    biogeog <- t(to_plot[[cou]])
    for (q in unique(biogeog[biogeog !=0])){
      hi<-biogeog
      hi[biogeog!=q]=0
      contour(lo, lt, hi, col = 'black',add=T, drawlabels = F, cex=0.284)
    }
    hide_arctic()
    axis_map0()
    dev.off()
  }
} 
  
pdf(family="Helvetica",paste('3D_phate_',type,'_',mod,'_',opt_n,'_woa.pdf', sep=''), width = 10, height = 10)
scatterplot3d(x = x, y = y, z= z, color = rgb(x_n, y_n, z_n, maxColorValue = 255), 
              xlab = "1st dim", ylab = "2nd dim", zlab = "3rd dim" )
dev.off()
