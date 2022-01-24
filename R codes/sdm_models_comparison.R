#!bin/usr/bin/env Rscript
# sdm_models_comparison.R
# Code to compare the outputs of 4 machine learning algorithms (randomforest, General additive models, neural networks and gradient boosting machines) 
# used to define environmental niches of the 'metacommunities'
# They are compared by calculating the pearson correlation coefficient of the predicted values of the training dataset
library("gbm")
library("dismo")
library("FactoMineR")
library("factoextra")
# library("readxl")
library("ggplot2")
library("reshape2")
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library('plotrix')
library('randomForest')
library('mgcv')
library('nnet')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
df <- readRDS('Genocenoses_env_parameters_woa_scaled.rds')
cambria <- readRDS('cambria.rds')
fractions <- c('180-2000', '20-180', '43952','0.8-5','0.22-3', '0-0.2')
variables <- c(6:11,13)
best_models <- read.table('best_selected_models.txt', header = T)
best_models<-as.data.frame(best_models)
N<-dim(best_models)[1] 
row.names(best_models)<-c(1:N)
models_bt <- readRDS('models_bt.rds')
models_gam <- readRDS('models_gam.rds')
models_nn <- readRDS('models_nn.rds')
models_rf <- readRDS('models_rf.rds')
good <- readRDS('good_models.rds')
cors <- NULL
indexes <- NULL
count <-1
for (fraction in fractions){
  df1 <- df[df$Fraction== fraction,]
  df1 <- df1[!is.na(df1$Genocenose),]
  for (k in best_models$Gen[best_models$Fraction==fraction]){
    df2 <- df1
    df2$Genocenose <- as.integer(df2$Genocenose == k)
    df2 <- df2[sample(1:nrow(df2)),]
    df2 <- as.data.frame(df2)
    for (i in 6:14){
      df2[,i] <- na.roughfix(df2[,i])
    }
    gam_model <- models_gam[[count]]
    bt_model <- models_bt[[count]]
    rf_model <- models_rf[[count]]
    nn_model <- models_nn[[count]]
    # Clculating the fitted values of the different models and calculating their pairwise correlation
    if (!is.na(gam_model)){
      test <- mgcv::predict.gam(gam_model, df2[,variables], type='response')
    }
    test1 <- stats::predict(rf_model, df2[,variables], type='prob')
    test1 <- test1[,2]
    if (!is.na(bt_model)){
      test2 <- gbm::predict.gbm(bt_model, df2[,variables], type='response', n.trees = bt_model$gbm.call$best.trees)
    }
    if (!is.na(nn_model)){
      test3 <-stats::predict(nn_model, df2[,variables], type='raw')[,2]
    }
    cors <- append(cors, fraction)
    cors <- append(cors, k)
    if (!is.na(bt_model)){
      cors <- append(cors, cor(test1, test2))
      if (!is.na(gam_model)){
        cors <- append(cors, cor(test, test2))
      } else{
        cors <- append(cors, NA)
      }
      if (!is.na(nn_model)){
        cors <- append(cors, cor(test3, test2))
      } else{
        cors <- append(cors, NA)
      }
    } else{
      cors <- append(cors, NA)
      cors <- append(cors, NA)
      cors <- append(cors, NA)
    }
    if (!is.na(gam_model)){
      cors <- append(cors, cor(test1, test))
      if (!is.na(nn_model)){
        cors <- append(cors, cor(test3, test))
      } else{
        cors <- append(cors, NA)
      }
    } else{
      cors <- append(cors, NA)
      cors <- append(cors, NA)
    }
    if (!is.na(nn_model)){
      cors <- append(cors, cor(test3, test1))
    } else{
      cors <- append(cors, NA)
    }
    count=count+1
  }
}
cors <- matrix(cors, ncol=8, byrow=T)
cors1<-NULL
for (i in c(3:8)){
  cors1 <- cbind(cors1, as.numeric(cors[,i]))
}

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
# Plotting the correaltion matrix
pdf(family="Helvetica",'correlation_sdm_models.pdf', height= 7, width=10)
br <- c(seq(-0.1, 1, length.out = 16))
heatmap.2(cors1, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1,margins=c(10,22), labCol=c('rf vs gbm', 'gam vs gbm','nn vs gbm','gam vs rf',  'gam vs nn', 'nn vs rf'),
          labRow=labs, col= greenred(15), breaks=br, symkey = F)
dev.off()
pdf(family="Helvetica",'correlation_sdm_models_1.pdf', height= 7, width=10)
br <- c(seq(-0.1, 1, length.out = 16))
heatmap.2(cors1, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1,margins=c(10,22), labCol=c('rf vs gbm', 'gam vs gbm','nn vs gbm','gam vs rf',  'gam vs nn', 'nn vs rf'),
          labRow=labs0, col= greenred(15), breaks=br, symkey = F, cexRow=1.2)
dev.off()

