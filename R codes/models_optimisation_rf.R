#!bin/usr/bin/env Rscript
# models_optimisation_bt.R
# Code to optimize hyperparameters of random forest for each metacommunities using 30 random cross-validations
library("gbm")
library("dismo")
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

df <- readRDS('Genocenoses_env_parameters_woa_scaled.rds')
fractions <- c('180-2000', '20-180', '43952','0.8-5','0.22-3', '0-0.2')

# rf optimisation grid: space of parameters tested
rfGrid <-  expand.grid(mtry = seq(from = 1, to = 8, by = 1),
                       ntree = c(seq(1, 10, by= 2) %o% 10^(2:3)))

best_models_rf <- function(id){
  optimisation_rf <- function(i, rfGrid, variables){
    flag <- TRUE
    rf_model <- NULL  
    df2$Genocenose <- as.factor(df2$Genocenose)
    # Fitting the model on the whole dataset
    set.seed(42)
    possibleError <- tryCatch(
      rf_model <- randomForest::randomForest(data=df2[,variables], x = df2[,variables], y=df2$Genocenose,importance=T, proximity =T, mtry=rfGrid[i,1], ntree = rfGrid[i,2]),
      error=function(e) flag <- FALSE
    )
    # Performing 30 cross-validations to calculate the mean AUC (area under ROC curve) which is the parameter we chose to optimize
    if (flag != FALSE){
      test3 <- stats::predict(rf_model, df2[,variables], type='prob')
      test3 <- test3[,2]
      TSSs <- 0
      cv_rf <- function(sample, i){
        df3 <- df2[sample,]
        set.seed(42)
        df3$Genocenose <- as.factor(df3$Genocenose)
        rf_model_cv <- randomForest::randomForest(data=df3[,variables], x = df3[,variables], y=df3$Genocenose,importance=T, proximity =T, mtry=rfGrid[i,1], ntree = rfGrid[i,2])
        preds <- stats::predict(rf_model_cv, df2[!(c(1:nrow(df2)) %in% sample),variables], type='class')
        d <- cbind(as.numeric(df2$Genocenose[!(c(1:nrow(df2)) %in% sample)]), as.numeric(preds)-1)
        pres <- d[d[,1]==1, 2]
        abs <- d[d[,1]==0, 2]
        sens <- sum(pres)/(sum(pres)+(length(pres)-sum(pres)))
        spec <- (length(abs)-sum(abs))/( (length(abs)-sum(abs)) +  sum(abs) )
        TSS <- sens+spec-1
        e <- dismo::evaluate(pres, abs)
        auc <- e@auc
        cor <- e@cor
        return(c(TSS, auc, cor))
      }
      
      set.seed(42)
      samples_list <- rep(list(), 30)
      df2$Genocenose <- as.numeric(levels(df2$Genocenose))[df2$Genocenose]
      for (u in 1:30){
        df3<-NULL
        while (sum(df3[,3], na.rm = T) < 2 | sum(df3[,3], na.rm = T) == sum(df2[,3], na.rm = T)){
          samples <- sample(nrow(df2), 0.75*nrow(df2))
          df3 <- df2[samples,]
        }
        samples_list[[u]] <- samples
      }
      score_list <- lapply(samples_list, FUN = cv_rf, i=i)
      
      
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
      rmse1 <- Metrics::rmse(test3[!is.na(test3)], as.numeric(df2$Genocenose[!is.na(test3)]))
    } else{
      TSSs <- NA
      AUCs <- NA
      CORs <- NA
      rmse1 <- NA
    }
    optimisation_rf <- c(TSSs,AUCs, CORs, rmse1)
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
  for (i in 1:dim(rfGrid)[1]){
    vec <- optimisation_rf(i, rfGrid = rfGrid, variables = variables)
    tss <- append(tss,vec[1])
    auce <- append(auce, vec[2])
    core <- append(core, vec[3])
    rmse <- append(rmse, vec[4])
  }
  # Finding the parameters combination for which the mean auc of the cross-validation is maximized
  g <- which.max(auce)
  if (length(g)>0){
    best_model<- c(fraction, k, rfGrid[g,1], rfGrid[g,2], tss[g], rmse[g], auce[g], core[g])
  } else{
    best_model<- c(fraction,k, NA,NA,NA,NA, NA, NA)
  }
  write(c(fraction, k), 'follow_rf.txt', append=T)
  return(best_model)
}
variables <- c(6:11,13)
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("df","rfGrid",  'fractions', 'variables'))
clusters <- sort(unique(df$id))
score_list <- parSapply(cl = cl, clusters, FUN = best_models_rf)
stopCluster(cl)
# Saving the results
best_models_rf1 <-t(score_list)
colnames(best_models_rf1) <- c('Fraction','Gen', 'mtry', 'ntree', 'tss', 'rmse', 'auc', 'cor')
best_models_rf1 <- as.data.frame(best_models_rf1)
write.table(best_models_rf1, 'best_models_rf.txt')
