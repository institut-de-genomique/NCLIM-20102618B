#!bin/usr/bin/env Rscript
library("gbm")
library('randomForest')
library('mgcv')
library('nnet')
library("dismo")
library("FactoMineR")
# library("readxl")
library('ggplot2')
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library("CDFt")
library('parallel')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
df <- readRDS("Genocenoses_env_parameters_woa_scaled.rds")

best_models <- read.table('best_selected_models.txt', header = T)
best_models <- as.data.frame(best_models)
N <-dim(best_models[1])
row.names(best_models)<-c(1:N)
variables <- c(6:11,13)
variables1 <- c('thetao','so', 'si','no3','po4', 'dfe')
fractions <- c('180-2000', '20-180', '43952', '0.8-5', '0.22-3', '0-0.2')
models_bt <- list()
models_nn <- list()
models_rf <- list()
models_gam <- list()
good <- NULL
count <- 1
for (frac in fractions){
  fraction <- frac
  df1 <- df[df$Fraction== fraction,]
  df1 <- df1[!is.na(df1$Genocenose),]
  
  goods <- best_models$Gen[best_models$Fraction==fraction]
  for (k in goods)
    local({
      df2 <- df1
      df2$Genocenose <- as.integer(df2$Genocenose == k)
      df2 <- as.data.frame(df2)
      set.seed(42)
      df2 <- df2[sample(1:nrow(df2)),]
      for (i in 6:14){
        df2[,i] <- randomForest::na.roughfix(df2[,i])
      }
      set.seed(42)
      if (!is.na(best_models$comp_bt[best_models$Fraction==fraction & best_models$Gen==k])){
        flags <- NULL
        set.seed(42)
        bt_model <- NULL
        possibleError <- tryCatch(
          bt_model <- dismo::gbm.step(data=df2, gbm.x = variables,gbm.y = 3,tree.complexity = best_models$comp_bt[best_models$Fraction==fraction & best_models$Gen==k],
                                      learning.rate = best_models$lr_bt[best_models$Fraction==fraction & best_models$Gen==k],
                                      n.minobsinnode = best_models$mo_bt[best_models$Fraction==fraction & best_models$Gen==k], 
                                      bag.fraction = 0.6, max.trees = 15000, verbose = F, plot.main = F),
          error=function(e) flag <- FALSE
        )
      } else{
        bt_model<-NA
      }
      models_bt[[count]]<<- bt_model
      if (best_models$ok_gam[count]==1){
          if (fraction!='0-0.2' & fraction != '43952'){
            set.seed(42)
            gam_model <- mgcv::gam(Genocenose ~ s(T)   + s(NO3) + s(Fe) + s(Phos) + s(Si) + s(SI_NO3) + s(Sal), data = df2, family = 'binomial',
                                   mehtod='ML', optimizer = c('outer', 'newton'), control = mgcv::gam.control(epsilon = 10^-6))
          } else{
            set.seed(42)
            gam_model <- mgcv::gam(Genocenose ~ s(T, k=2) + s(NO3, k=2) + s(Fe, k=2) + s(Phos,k=2) + s(Si,k=2) + s(SI_NO3,k=2) + s(Sal,k=2), 
                                   data = df2, family = 'binomial',mehtod='ML', optimizer = c('outer', 'newton'), control = mgcv::gam.control(epsilon = 10^-6))
          }
      } else{
          gam_model <- NA
      }
      models_gam[[count]] <<- gam_model
      set.seed(42)
      nn_model <- nnet::nnet(nnet::class.ind(df2$Genocenose) ~ T  + Sal + Si + NO3 + Phos  + Fe  + SI_NO3, data = df2, softmax=T,
                               size=best_models$size_nn[best_models$Fraction==fraction & best_models$Gen==k], 
                               decay=best_models$decay_nn[best_models$Fraction==fraction & best_models$Gen==k],
                               maxit=best_models$mxit_nn[best_models$Fraction==fraction & best_models$Gen==k])
      models_nn[[count]] <<- nn_model
      df2$Genocenose <- as.factor(df2$Genocenose)
      set.seed(42)
      rf_model <- randomForest::randomForest(data=df2[,variables], x = df2[,variables], y=df2$Genocenose,importance=T, proximity =T, 
                                             mtry=best_models$mtry_rf[best_models$Fraction==fraction & best_models$Gen==k], 
                                             ntree=best_models$ntree_rf[best_models$Fraction==fraction & best_models$Gen==k])
      models_rf[[count]] <<- rf_model
      g <- paste(as.character(k), '_',fraction, ' ', sep='')
      good <<- append(good, g)
      count <<- count+1
    })
}
saveRDS(good, 'good_models.rds')
saveRDS(models_bt, 'models_bt.rds')
saveRDS(models_nn, 'models_nn.rds')
saveRDS(models_rf, 'models_rf.rds')
saveRDS(models_gam, 'models_gam.rds')
