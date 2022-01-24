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

type = commandArgs(trailingOnly = T)[1]
if (type=='dom'){
  raw_data_woa <- readRDS('pred_domWOA.rds')
} else if (type=='discrete'){
  raw_data_woa <- readRDS('pred_domWOA.rds')
  raw_data_woa[raw_data_woa!=0]<-1
  raw_data_woa <- unique(raw_data_woa)
} else if (type=='prob'){
  raw_data_woa <- readRDS('model-mean_pred_woa_list.rds')
}
new_stations_1 <- readRDS('new_stations_1deg.rds')

if (type=='discrete'){
  knn =5
  dc=10
} else{
  knn=1000
  dc=10
}
set.seed(2)
phate_fit <- phate(raw_data_woa, ndim = 3, knn = knn, decay = dc)
saveRDS(phate_fit, paste('phate_fit_woa_',type,'.rds', sep=''))


