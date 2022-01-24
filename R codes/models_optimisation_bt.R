#!bin/usr/bin/env Rscript
# models_optimisation_bt.R
# Code to optimize hyperparameters of gradient boosting machines for each metacommunities using 30 random cross-validations
library("gbm")
library("dismo")
library("FactoMineR")
library("factoextra")
#library("readxl")
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
df <- readRDS("Genocenoses_env_parameters_woa_scaled.rds")
fractions <- c('180-2000', '20-180', '43952','0.8-5','0.22-3', '0-0.2')


# bt optimisation grid: space of parameters tested
btGrid <- expand.grid(interaction.depth=c(1,3,5),
                         shrinkage=c(0.001,0.01),
                         n.minobsinnode=seq(1,10,2))
best_models_bt <- function(id){
  optimisation_bt <- function(i, btGrid, variables){
    flags <- NULL
    bt_models <- list()
    count <- 1
    # 3 tests of fitting the model are tried because gbm.step algorithm depends on initial conditions
    set.seed(42)
    for (u in 1:3){
      flag <- TRUE
      bt_model <- NULL
      possibleError <- tryCatch(
        bt_model <- dismo::gbm.step(data=df2, gbm.x = variables, gbm.y = 3, tree.complexity = btGrid[i,1],learning.rate = btGrid[i,2],  n.minobsinnode = btGrid[i,3], bag.fraction = 0.6, max.trees = 15000, verbose = F, plot.main = F),
        error=function(e) flag <- FALSE
      )
      if (is.null(bt_model) & flag==T){
        flag <- F
      }
      bt_models[[count]] <- bt_model
      flags <- append(flags, flag)
      count <- count+1
    }
    index <- which(flags==T)[1]
    if (!is.na(index)){
      bt_model<-bt_models[[index]]
      write(index, 'indexes.txt', append=T)
    } else{
      bt_model <- NULL
    }
    
    # Performing 30 random cross-validations to calculate the mean AUC (area under ROC curve) which is the parameter we chose to optimize
    if (!is.null(bt_model)){
      test3 <-  gbm::predict.gbm(bt_model, df2[,variables], type='response', n.trees = bt_model$gbm.call$best.trees)
      TSSs <- 0
      cv_bt <- function(sample, i){
        df3 <- df2[sample,]
        flag0 <- TRUE
        bt_model_cv <- NULL
        set.seed(42)
        possibleError <- tryCatch(
          bt_model_cv <- dismo::gbm.step(data=df3, gbm.x = variables, gbm.y = 3, tree.complexity = btGrid[i,1],learning.rate = btGrid[i,2],  n.minobsinnode = btGrid[i,3], bag.fraction = 0.6, max.trees = 15000, verbose = F, plot.main = F),
          error=function(e) flag0 <- FALSE
        )
        if (flag0 != FALSE & !is.null(bt_model_cv)){
          preds <- stats::predict(bt_model_cv, df2[!(c(1:nrow(df2)) %in% sample),variables], type='link', n.trees = bt_model_cv$gbm.call$best.trees)
          preds[preds>0]=1
          preds[preds<0]=0
          d <- cbind(df2$Genocenose[!(c(1:nrow(df2)) %in% sample)], preds)
          pres <- d[d[,1]==1, 2]
          abs <- d[d[,1]==0, 2]
          sens <- sum(pres)/(sum(pres)+(length(pres)-sum(pres)))
          spec <- (length(abs)-sum(abs))/( (length(abs)-sum(abs)) +  sum(abs) )
          TSS <- sens+spec-1
          e <- dismo::evaluate(pres, abs)
          auc <- e@auc
          cor <- e@cor
          cv_bt <- c(TSS, auc, cor)
        } else{
          cv_bt <- NA
        }
      }
      set.seed(42)
      samples_list <- rep(list(), 30)
      for (u in 1:30){
        df3<-NULL
        while (sum(df3[,3], na.rm = T) < 2 | sum(df3[,3], na.rm = T) == sum(df2[,3], na.rm = T)){
          samples <- sample(nrow(df2), 0.85*nrow(df2))
          df3 <- df2[samples,]
        }
        samples_list[[u]] <- samples
      }
      score_list <- lapply(samples_list, FUN = cv_bt, i=i)
      
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
    if (!is.null(bt_model)){
      rmse1 <- Metrics::rmse(test3[!is.na(test3)], df2$Genocenose[!is.na(test3)])
    } else{
      rmse1<- NA
    }
    optimisation_bt <- c(TSSs,AUCs, CORs, rmse1)
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
  tss <- NULL
  core <- NULL
  auce <- NULL
  rmse <- NULL
  for (i in 1:dim(btGrid)[1]){
    set.seed(42)
    vec <- optimisation_bt(i, btGrid = btGrid, variables=variables)
    tss <- append(tss,vec[1])
    auce <- append(auce, vec[2])
    core <- append(core, vec[3])
    rmse <- append(rmse, vec[4])
  }
  # Finding the parameters combination for which the mean auc of the cross-validation is maximized
  g <- which.max(auce)
  if (length(g)==0){
    if (sum(rmse, na.rm = T) != 0){
      g <- which.min(rmse)
      best_model <-  c(fraction,k,  btGrid[g,1], btGrid[g,2],btGrid[g,3],NA, rmse[g], NA, NA)
    } else{
      best_model <-  c(fraction,k,  NA, NA,NA,NA, NA, NA, NA)
    }
  } else{
    best_model <- c(fraction,k, btGrid[g,1], btGrid[g,2],btGrid[g,3],  tss[g], rmse[g], auce[g], core[g])
  }
  write(c(fraction, k), 'follow_bt.txt', append=T)
  return(best_model)
}
variables <- c(6:11,13)
# Initiate cluster
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("df","btGrid",  'fractions', 'variables'))
clusters <- sort(unique(df$id))
score_list <- parSapply(cl = cl, clusters, FUN = best_models_bt)
stopCluster(cl)
# Saving the results
best_models_bt1 <-t(score_list)
colnames(best_models_bt1) <- c('Fraction','Gen', 'comp', 'lr','mo', 'tss', 'rmse', 'auc', 'cor')
best_models_bt1 <- as.data.frame(best_models_bt1)
write.table(best_models_bt1, 'best_models_bt.txt')

