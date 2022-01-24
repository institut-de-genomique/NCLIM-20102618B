#!bin/usr/bin/env Rscript
library("gbm")
library('randomForest')
library('mgcv')
library('nnet')
library("dismo")
library("FactoMineR")
# library("factoextra")
# library("readxl")
# library("ggplot2")
# library("reshape2")
library("gplots")
# library("plotly")
library("stringr")
# library("caret")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library("CDFt")
library('plotrix')
library('png')
library('grid')

source('CDFt_BIS.R')
model0<-'model-mean'
nvar <- 7
mean_sds <- readRDS('means_sds_woa_tara.rds')
black_caspien <- expand.grid(seq(28.5, 60.5, 1),seq(37.5, 48.5, 1))
med_sea <- expand.grid(seq(-2.5, 35.5, 1),seq(30.5, 45.5, 1))
closed_sea <- rbind(black_caspien, med_sea)
f2006 = ncdf4::nc_open(paste('all_7_mon_',model0,'_rcp85_200601_201512_clim.nc',sep=''))
f2090 = ncdf4::nc_open(paste('all_7_mon_',model0,'_rcp85_209001_209912_clim.nc',sep=''))
i =1
lati<-seq(-89.5, 89.5, 1)
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
lati_sorted<-seq(-89.5, 60, 1)
longi_sorted =sort(longi)
variables1 <- c('thetao','so', 'si','no3','po4', 'dfe')
data <- lapply(1:nvar, matrix, data= NA, nrow=360, ncol=180)
data1 <- lapply(1:nvar, matrix, data= NA, nrow=360, ncol=180)
local(
  for (v in variables1){
    f2 <- ncdf4::ncvar_get(f2006,v)
    f3 <- ncdf4::ncvar_get(f2090,v)
    
    ## si-no3 matrix
    if (v =='no3'){
      mat = matrix(data=NA, nrow=360, ncol=180)
      mat1 = matrix(data=NA, nrow=360, ncol=180)
      for (k in 1:360){
        for (l in 1:150){
          if ( !( is.na(sum(f2[k,l,], na.rm = T ) ) ) & sum(f2[k,l,], na.rm = T) != 0 ){
            Si_no3 = ( max(f2[k,l,], na.rm=TRUE) - min(f2[k,l,], na.rm=TRUE) )/20.3490 ## Highest SI-no3 in tara
            Si_no3_1 = ( max(f3[k,l,], na.rm=TRUE) - min(f3[k,l,], na.rm=TRUE) )/20.3490
            mat[k,l]=Si_no3
            mat1[k,l]=Si_no3_1
          }
        }
      }
      data[[7]]<<-mat
      data1[[7]]<<-mat1
    }
    
    a <- apply(f2, c(1,2), mean, na.rm = TRUE)
    a1 <- apply(f3, c(1,2), mean, na.rm =TRUE)
    
    data[[i]] <<- a
    data1[[i]] <<- a1
    i =i+1
  }
)
data_scaled <- data
data1_scaled <- data1
for (i in 1:nvar){
  data_scaled[[i]] <- (data_scaled[[i]]-mean_sds[[i+5]][1])/mean_sds[[i+5]][2]
  data1_scaled[[i]] <- (data1_scaled[[i]]-mean_sds[[i+5]][1])/mean_sds[[i+5]][2]
}
saveRDS(data, 'model-mean_2006-15.rds')
saveRDS(data1, 'model-mean_2090-99.rds')
saveRDS(data_scaled, 'model-mean_2006-15_scaled.rds')
saveRDS(data1_scaled, 'model-mean_2090-99_scaled.rds')

fwoa <- ncdf4::nc_open('t_so_si_no3_po4_woa13_annual_mean_2005_2012.nc')
no3_woa_m <- ncdf4::nc_open('no3_mon_woa-13-an_2005_2012_clim.nc')
dfe <- ncdf4::nc_open('FED_PISCES2_1y_FEMIP.nc')
variables_ <- c('t_an','s_an', 'i_an','n_an', 'p_an')
i =1
data_w <- lapply(1:nvar, matrix, data= NA, nrow=360, ncol=180)
local(
  for (v in variables_){
    f2 <- ncdf4::ncvar_get(fwoa,v)
    data_w[[i]] <<- f2
    i =i+1
  })
dfe <- ncdf4::ncvar_get(dfe, 'IRON')
fe <- dfe[,,1]*10^6
data_w[[6]]=fe
## si-no3
no3_m <- ncdf4::ncvar_get(no3_woa_m, 'no3')
mat2 <- matrix(data=NA, nrow=360, ncol=180)
cos_lat <- matrix(data=NA, nrow=360, ncol=180)
local(
  for (k in 1:360){
    for (l in 1:150){
      if (  !(is.na(sum(no3_m[k,l,], na.rm = T))) & sum(no3_m[k,l,], na.rm = T) != 0 ){
        Si_no3 = ( max(no3_m[k,l,], na.rm=TRUE) - min(no3_m[k,l,], na.rm=TRUE) )/20.3490 ## Highest SI-no3 in tara
        mat2[k,l]<<-Si_no3
        cos_lat[k,l]<<-cos(lati[l]*2*pi/360)
      }
    }
  })
data_w[[7]]=mat2


longitude <- matrix(NA,  nrow=360, ncol=180)
latitude <- matrix(NA,  nrow=360, ncol=180)
k=1
for (i in 1:360){
 for (j in 1:180){
  longitude[i,j]=longi[i]
  latitude[i,j]=lati[j]
  k=k+1
 }
}
saveRDS(longitude, 'longitude.rds')
saveRDS(latitude, 'latitude.rds')

## Correction
data_p <- lapply(1:nvar, matrix, data= NA, nrow=360, ncol=180)
data_f <- lapply(1:nvar, matrix, data= NA, nrow=360, ncol=180)
data_p_n <- lapply(1:nvar, matrix, data= NA, nrow=360, ncol=180)
data_f_n <- lapply(1:nvar, matrix, data= NA, nrow=360, ncol=180)
local(
  for (l in 1:nvar){
    f_h = data[[l]]
    f_f = data1[[l]]
    f_w =data_w[[l]]
    
    f_h[is.nan(f_h)]=NA
    f_f[is.nan(f_f)]=NA
    f_w[is.nan(f_w)]=NA
    
    obs <- NULL
    pres <- NULL
    futur <- NULL
    lats <- NULL
    longs <- NULL
    
    for (i in 1:360){
      for (j in 1:150){
        if ( !(is.na(f_h[i,j])) & !(is.na(f_f[i,j])) & !(is.na(f_w[i,j])) & !(longi[i] %in% closed_sea[[1]] & lati[j] %in% closed_sea[[2]])){
          if (l !=4){
            obs <- append(obs, as.numeric(f_w[i,j]))
            pres <- append(pres, as.numeric(f_h[i,j]))
            futur <- append(futur, as.numeric(f_f[i,j]))
            longs <- append(longs, i)
            lats <- append(lats, j)
          } else if (l==4 & f_h[i,j]<35){
            obs <- append(obs, as.numeric(f_w[i,j]))
            pres <- append(pres, as.numeric(f_h[i,j]))
            futur <- append(futur, as.numeric(f_f[i,j]))
            longs <- append(longs, i)
            lats <- append(lats, j)
          }
        }
      }
    }
    
    ## bias correction
    if (l==4){
      test0 <- CDFt_BIS(obs, pres, pres, npas=1000)
      test <- CDFt_BIS(obs, pres, futur, npas=1000)
    } else{
      test0 <- CDFt(obs, pres, pres, npas=1000)
      test <- CDFt(obs, pres, futur, npas=1000)
    }
    
    
    for (k in 1:length(longs)){
      data_f_n[[l]][longs[k], lats[k]] <<- test$DS[k]
      data_p_n[[l]][longs[k], lats[k]] <<- test0$DS[k]
      if (test$DS[k]<=0 & l != 1){
        data_f[[l]][longs[k], lats[k]] <<- min(test$DS[test$DS>=0], na.rm =T)
      } else{
        data_f[[l]][longs[k], lats[k]] <<- test$DS[k]
      }
      if (test0$DS[k]<=0 & l !=1){
        data_p[[l]][longs[k], lats[k]] <<- min(test0$DS[test0$DS>=0], na.rm =T)
      } else{
        data_p[[l]][longs[k], lats[k]] <<- test0$DS[k]
      }
    } 
})
rm(f2006, f2090, fwoa)
data_p_scaled <- data_p
data_f_scaled <- data_f
data_w_scaled <- data_w
for (i in 1:nvar){
  data_p_scaled[[i]] <- (data_p_scaled[[i]]-mean_sds[[i+5]][1])/mean_sds[[i+5]][2]
  data_f_scaled[[i]] <- (data_f_scaled[[i]]-mean_sds[[i+5]][1])/mean_sds[[i+5]][2]
  data_w_scaled[[i]] <- (data_w_scaled[[i]]-mean_sds[[i+5]][1])/mean_sds[[i+5]][2]
  print(mean(data_p_scaled[[i]], na.rm = T))
  print(mean(data_f_scaled[[i]], na.rm = T))
  print(mean(data_w_scaled[[i]], na.rm = T))
  print(sd(data_p_scaled[[i]], na.rm = T))
  print(sd(data_f_scaled[[i]], na.rm = T))
  print(sd(data_w_scaled[[i]], na.rm = T))
}
saveRDS(data_p, 'data_p.rds')
saveRDS(data_f, 'data_f.rds')
saveRDS(data_w, 'data_w.rds')
saveRDS(data_p_n, 'data_p_n.rds')
saveRDS(data_f_n, 'data_f_n.rds')
saveRDS(data_p_scaled, 'data_p_scaled.rds')
saveRDS(data_f_scaled, 'data_f_scaled.rds')
saveRDS(data_w_scaled, 'data_w_scaled.rds')

