#!bin/usr/bin/env Rscript
# models_optimisation_bt.R
# Code to perform cross-validation of the General additive model for each metacommunities using 30 random cross-validations and calculte the model performance (AUC)
library("gbm")
library("dismo")
#library("readxl")
library("ggplot2")
library("reshape2")
library("gplots")
# library("plotly")
library("stringr")
# library("caret")
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
df <- readRDS("Genocenoses_env_parameters_woa_scaled.rds")
fractions <- c('180-2000', '20-180', '43952','0.8-5','0.22-3', '0-0.2')

best_models_gam <- function(id){
  optimisation_gam <- function(i, variables){
    gam_model <- NULL 
    flag <- TRUE
    df2$Genocenose <- as.factor(df2$Genocenose)
    # For size fraction '0-0.2' and '5-20' ('43952'),
    # the number of splined is fixed to 2 per parameters otherwise the model would have more parameters than observations
    set.seed(42)
    if (fraction != '0-0.2' & fraction!='43952'){
      possibleError <- tryCatch(
        gam_model <- mgcv::gam(Genocenose ~ s(T)   + s(NO3) + s(Fe) + s(Phos) + s(Si) + s(SI_NO3) + s(Sal), data = df2, family = 'binomial',mehtod='ML', optimizer = c('outer', 'newton'), control = mgcv::gam.control(epsilon = 10^-6)),
        error=function(e) flag <- FALSE
      )
    } else{
      possibleError <- tryCatch(
        gam_model <- mgcv::gam(Genocenose ~ s(T, k=2)   + s(NO3, k=2) + s(Fe, k=2) + s(Phos,k=2) + s(Si,k=2) + s(SI_NO3,k=2) + s(Sal,k=2), data = df2, family = 'binomial',mehtod='ML', optimizer = c('outer', 'newton'), control = mgcv::gam.control(epsilon = 10^-6)),
        error=function(e) flag <- FALSE
      )
    }
    if (flag != FALSE & !is.null(gam_model)){
      test3 <- mgcv::predict.gam(gam_model, df2[,variables], type='response')
      TSSs <- 0
      cv_gam <- function(sample){
        df3 <- df2[sample,]
        set.seed(42)
        df3$Genocenose <- as.factor(df3$Genocenose)
        print(sum(as.numeric(df3$Genocenose)-1))
        flag <- TRUE
        gam_model_cv <- NULL
        if (fraction != '0-0.2' & fraction != '43952'){
          possibleError <- tryCatch(
            gam_model_cv <- mgcv::gam(Genocenose ~ s(T)   + s(NO3) + s(Fe) + s(Phos) + s(Si) + s(SI_NO3) + s(Sal), data = df3, family = 'binomial',mehtod='ML', optimizer = c('outer', 'newton'), control = mgcv::gam.control(epsilon = 10^-6)),
            error=function(e) flag <- FALSE
          )
        } else{
          possibleError <- tryCatch(
            gam_model_cv <- mgcv::gam(Genocenose ~ s(T, k=2)   + s(NO3, k=2) + s(Fe, k=2) + s(Phos,k=2) + s(Si,k=2) + s(SI_NO3,k=2) + s(Sal,k=2), data = df3, family = 'binomial',mehtod='ML', optimizer = c('outer', 'newton'), control = mgcv::gam.control(epsilon = 10^-6)),
            error=function(e) flag <- FALSE
          )
        }
        if (flag==TRUE & !is.null(gam_model_cv)){
          preds <- mgcv::predict.gam(gam_model_cv, df2[!(c(1:nrow(df2)) %in% sample),variables], type='link')
          preds[preds>0]=1
          preds[preds<0]=0
          d <- cbind(as.numeric(df2$Genocenose[!(c(1:nrow(df2)) %in% sample)]), as.numeric(preds)-1)
          pres <- d[d[,1]==1, 2]
          abs <- d[d[,1]==0, 2]
          sens <- sum(pres)/(sum(pres)+(length(pres)-sum(pres)))
          spec <- (length(abs)-sum(abs))/( (length(abs)-sum(abs)) +  sum(abs) )
          TSS <- sens+spec-1
          e <- dismo::evaluate(pres, abs)
          auc <- e@auc
          cor <- e@cor
        } else{
          TSS<-NA
          auc<-NA
          cor<-NA
        }
        cv_gam <- c(TSS, auc, cor)
      }
      
      set.seed(42)
      samples_list <- rep(list(), 30)
      df2$Genocenose <- as.numeric(levels(df2$Genocenose))[df2$Genocenose]
      for (u in 1:30){
        df3<-NULL
        while (sum(df3[,3], na.rm = T) < 2 | sum(df3[,3], na.rm = T) == sum(df2[,3], na.rm = T)){
          samples <- sample(nrow(df2), 0.85*nrow(df2))
          df3 <- df2[samples,]
        }
        samples_list[[u]] <- samples
      }
      score_list <- lapply(samples_list, FUN = cv_gam)
      
      tss_list <- NULL
      auc_list <- NULL
      cor_list <- NULL
      for (vector in score_list){
        tss_list <- append(tss_list,vector[1])
        auc_list <- append(auc_list,vector[2])
        cor_list <- append(cor_list,vector[3])
      }
      TSSs <- mean(tss_list, na.rm=T)
      AUCs <- mean(auc_list, na.rm=T)
      CORs <- mean(cor_list, na.rm=T)
    } else{
      TSSs <- NA
      AUCs <- NA
      CORs <- NA
      rmse1 <- NA
    }
    if (!is.null(gam_model)){
      rmse1 <- Metrics::rmse(test3[!is.na(test3)], as.numeric(df2$Genocenose[!is.na(test3)]))
    } else{
      rmse1<- NA
    }
    optimisation_gam <- c(TSSs,AUCs, CORs, rmse1)
  }
  fraction <- strsplit(id, '_')[[1]][1]
  k <- as.integer(strsplit(id, '_')[[1]][2])
  df1 <- df[df$Fraction== fraction,]
  df1 <- df1[!is.na(df1$Genocenose),]
  
  df2 <- df1
  df2$Genocenose <- as.integer(df2$Genocenose == k)
  
  set.seed(42)
  df2 <- df2[sample(1:nrow(df2)),]
  df2 <- as.data.frame(df2)
  
  for (i in 6:14){
    df2[,i] <- randomForest::na.roughfix(df2[,i])
  }     
  i=1
  vec <- optimisation_gam(i, variables = variables)
  tss <- vec[1]
  auce <- vec[2]
  core <- vec[3]
  rmse <- vec[4]
  g <- which.max(tss)
  if (length(g)>0){
    best_model<- c(fraction,k, tss, rmse, auce, core)
  } else{
    if (!is.na(rmse)){
      best_model <- c(fraction,k, NA,rmse,NA, NA)
    } else{
      best_model <-  c(fraction,k, NA,NA, NA, NA)
    }
  }
  write(c(fraction, k), 'follow_gam.txt', append=T)
  return(best_model)
}

variables <- c(6:11,13)
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("df",  'fractions', 'variables'))
clusters <- sort(unique(df$id))
score_list <- parSapply(cl = cl,clusters, FUN = best_models_gam)
stopCluster(cl)
# Saving the results
best_models_gam1 <-t(score_list)
colnames(best_models_gam1) <- c('Fraction','Gen', 'tss', 'rmse', 'auc', 'cor')
best_models_gam1 <- as.data.frame(best_models_gam1)
write.table(best_models_gam1, 'best_models_gam.txt')
