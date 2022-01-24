library('ncdf4')
library('maps')
library("readxl")
library('mapproj')
library('mapplots')
library('SDMTools')
library('stringr')
library('ggplot2')
library('reshape2')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")

df <- read_xlsx('Genocenoses_env_parameters_all_tara.xlsx')

data06 <- readRDS('data_p.rds')
data90 <- readRDS('data_f.rds')

variables <- c('thetao','so', 'si','no3','po4')


long<-seq(0.5,359.5,1)
long[181:360]=as.numeric(long[181:360])-360
lat<-seq(-89.5, 89.5, 1)


find_param1 <- function(lg, lt, data){
  val <- mean(data[(lg-1):(lg+1),(lt-1):(lt+1) ], na.rm=T)
  if (is.nan(val) | is.na(val)){
    val <- mean(data[(lg-2):(lg+2),(lt-2):(lt+2) ], na.rm=T)
  }
  return(val)
}

df1 <- NULL
df2 <- NULL
stations <- unique(df$Station)
stations <- stations[grepl('SUR', stations)]
for (st in stations){
  d <- df[df$Station == st,]
  n <- nrow(d)
  
  month<- df$month[df$Station==st][1]
  latitude <- df$Lat[df$Station==st][1]
  longitude <- df$Long[df$Station==st][1]
  
  close_lat <- abs(lat-latitude)
  lat_ind <- which.min(close_lat)
  close_long <- abs(long-longitude)
  long_ind <- which.min(close_long)
  
  values <- NULL
  values <- append(values, c(st, as.numeric(latitude), as.numeric(longitude)))
  ## Temperature,so, NO3, Si, PO4, SI_NO3
  for (z in 1:7){
    v <- data06[[z]][long_ind,lat_ind]
    v1 <- data90[[z]][long_ind,lat_ind]
    if (is.na(v)){
      v <- find_param1(long_ind, lat_ind, data06[[z]])
    }
    if (is.na(v1)){
      v1 <- find_param1(long_ind, lat_ind, data90[[z]])
    }
    v <- as.numeric(v)
    values <-append(values, v)
    values <-append(values, v1)
  }
  
  df1 <- rbind(df1, values)
  
  d<- as.data.frame(d)
  for (u in 1:n){
    Fraction <- d[u,2][[1]]
    Genocenose <- d[u,3][[1]]
    values1 <- c(values[1], Fraction, Genocenose, values[2:length(values)])
    df2 <- rbind(df2, values1)
  }
}



df1 <- as.data.frame(df1)
colnames(df1) <- c('Station', 'Lat', 'Long','T_06','T_90','Sal_06','Sal_90',
                   'Si_06','Si_90','NO3_06','NO3_90', 'Phos_06','Phos_90',
                   'Fe_06','Fe_90' , 'SI_NO3_06','SI_NO3_90')#, 'SI_NO3')
df2 <- as.data.frame(df2)
colnames(df2) <- c('Station','Fraction','Genocenose', 'Lat', 'Long',
                   'T_06','T_90','Sal_06','Sal_90','Si_06','Si_90','NO3_06','NO3_90', 
                   'Phos_06','Phos_90','Fe_06','Fe_90' , 'SI_NO3_06','SI_NO3_90')#, 'SI_NO3')


write.table(df1, 'env_parameters_all_2006-15_2090-99.txt', row.names = F)
write.table(df2, 'Genocenoses_env_parameters_all_2006-15_2090-99.txt', row.names = F)

