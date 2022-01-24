library('DALEX')
library('readxl')
library('gbm')
library('mgcv')
library('nnet')
library('RColorBrewer')
library('ingredients')
#setwd("~/MobaXterm/home/final_model_2")
source('vioplot.R')
cambria <- readRDS('cambria.rds')

models_bt <- readRDS('models_bt.rds')
models_gam <- readRDS('models_gam.rds')
models_rf <- readRDS('models_rf.rds')
models_nn <- readRDS('models_nn.rds')

pred_bt <- function(model, newdata)  {
  results <- gbm::predict.gbm(model, newdata, type='response', n.trees = bt_model$gbm.call$best.trees)
  return(results)
}
pred_gam <- function(model, newdata)  {
  results <- as.vector(mgcv::predict.gam(model, newdata, type='response'))
  return(results)
}
pred_nn <- function(model, newdata)  {
  results <- stats::predict(model, newdata, type='raw')[,2]
  return(results)
}
pred_rf <- function(model, newdata)  {
  results <- stats::predict(model, newdata, type='prob')[,2]
  return(results)
}

df <- readRDS('Genocenoses_env_parameters_woa_scaled.rds')
best_models <- read.table('best_selected_models.txt', header = T)
best_models <- as.data.frame(best_models)
n_mods <- dim(best_models)[1]
row.names(best_models)<-c(1:n_mods)
variables <- c(6:11, 13)
N <- length(variables)
data_rl = data.frame()
for (j in 1:n_mods){
  fraction <- best_models$Fraction[j]
  k=best_models$Gen[j]
  df1 <- df[df$Fraction== fraction,]
  df1 <- df1[!is.na(df1$Genocenose),]
  df2 <- df1
  df2$Genocenose <- as.integer(df2$Genocenose == k)
  df2 <- as.data.frame(df2)
  df2 <- df2[sample(1:nrow(df2)),]
  for (i in 6:14){
    df2[,i] <- randomForest::na.roughfix(df2[,i])
  }
  bt_model = models_bt[[j]]
  gam_model = models_gam[[j]]
  rf_model = models_rf[[j]]
  nn_model = models_nn[[j]]
  
  if (!is.na(bt_model[1])){
    exp_bt <- explain(bt_model, data = df2[, variables], y=df2$Genocenose,
                       predict_function = pred_bt)
    imp_bt <- feature_importance(exp_bt, loss_function = loss_root_mean_square)
    rel_inf_bt<- imp_bt$dropout_loss[2:(N+1)]*100/sum(imp_bt$dropout_loss[2:(N+1)])
    vars_bt <- as.character(imp_bt$variable[2:(N+1)])
  } else{
    rel_inf_bt <- rep(NA, N)
  }
  if (!is.na(gam_model[1])){
    exp_gam <- explain(gam_model, data = df2[, variables], y=df2$Genocenose,
                      predict_function = pred_gam)
    imp_gam <- feature_importance(exp_gam, loss_function = loss_root_mean_square)
    rel_inf_gam <- imp_gam$dropout_loss[2:(N+1)]*100/sum(imp_gam$dropout_loss[2:(N+1)])
    vars_gam <- as.character(imp_gam$variable[2:(N+1)])
  } else{
    rel_inf_gam <- rep(NA, N)
  }
  exp_rf <- explain(rf_model, data = df2[, variables], y=df2$Genocenose,
                    predict_function = pred_rf)
  exp_nn <- explain(nn_model, data = df2[, variables], y=df2$Genocenose,
                    predict_function = pred_nn)
  
  imp_rf <- feature_importance(exp_rf, loss_function = loss_root_mean_square)
  imp_nn <- feature_importance(exp_nn, loss_function = loss_root_mean_square)
  rel_inf_rf <- imp_rf$dropout_loss[2:(N+1)]*100/sum(imp_rf$dropout_loss[2:(N+1)])
  vars_rf <- as.character(imp_rf$variable[2:(N+1)])
  rel_inf_nn <- imp_nn$dropout_loss[2:(N+1)]*100/sum(imp_nn$dropout_loss[2:(N+1)])
  vars_nn <- as.character(imp_nn$variable[2:(N+1)])
  rel_inf_bt <- rel_inf_bt[match(vars_nn, vars_rf)]
  rel_inf_gam <- rel_inf_gam[match(vars_nn, vars_gam)]
  rel_inf_rf<- rel_inf_rf[match(vars_nn, vars_rf)]
  vars <- vars_nn
  if (j==1){
    data_rl  <- data.frame(rep(fraction, N), rep(k,N),vars, rel_inf_gam, rel_inf_nn, rel_inf_rf, rel_inf_bt)
  } else{
    data_rl <- rbind(data_rl,
                     data.frame(rep(fraction, N), rep(k,N),vars, 
                                rel_inf_gam, rel_inf_nn, rel_inf_rf, rel_inf_bt))
  }
  print(j)
}
data_rl$Inf.<-apply(data_rl[,4:7], 1, mean, na.rm=T)
data_rl$sd <-apply(data_rl[,4:7], 1, sd, na.rm=T)


x1 <- data_rl$Inf.[data_rl$vars=='T']
x2 <- data_rl$Inf.[data_rl$vars=='Sal']
x3 <- data_rl$Inf.[data_rl$vars=='Si']
x4 <- data_rl$Inf.[data_rl$vars=='NO3']
x5 <- data_rl$Inf.[data_rl$vars=='Phos']
x6 <- data_rl$Inf.[data_rl$vars=='Fe']
# x7 <- data_rl$Inf.[data_rl$vars=='Chla']
x8 <- data_rl$Inf.[data_rl$vars=='SI_NO3']
colores = brewer.pal(N, "Dark2")
pdf(family="Helvetica",'relative_influence_violin_dalex.pdf')
vioplot(x1, x2, x3,x4,x5,x6,x8, names=c('T', 'Sal', 'Si', 'NO3', 'PO4', 'Fe','SI_NO3'), col =colores, ylim=c(0,100))
dev.off()
for (frac in unique(data_rl$rep.fraction..N.)){
  x1 <- data_rl$Inf.[data_rl$vars=='T' & data_rl$rep.fraction..N.==frac]
  x2 <- data_rl$Inf.[data_rl$vars=='Sal' & data_rl$rep.fraction..N.==frac]
  x3 <- data_rl$Inf.[data_rl$vars=='Si' & data_rl$rep.fraction..N.==frac]
  x4 <- data_rl$Inf.[data_rl$vars=='NO3' & data_rl$rep.fraction..N.==frac]
  x5 <- data_rl$Inf.[data_rl$vars=='Phos' & data_rl$rep.fraction..N.==frac]
  x6 <- data_rl$Inf.[data_rl$vars=='Fe' & data_rl$rep.fraction..N.==frac]
  # x7 <- data_rl$Inf.[data_rl$vars=='Chla' & data_rl$rep.fraction..N.==frac]
  x8 <- data_rl$Inf.[data_rl$vars=='SI_NO3' & data_rl$rep.fraction..N.==frac]
  colores = brewer.pal(N, "Dark2")
  pdf(family="Helvetica",paste('Relative_influences_violin_dalex_',frac,'.pdf', sep=''))
  vioplot(x1, x2, x3,x4,x5,x6,x8, names=c('T', 'Sal', 'Si', 'NO3', 'PO4', 'Fe','SI_NO3'), col =colores)
  dev.off()
}
data_glob_relinf <- NULL
for (frac in unique(data_rl$rep.fraction..N.)){
  dt <- data_rl[data_rl$rep.fraction..N.==frac,]
  print(pairwise.wilcox.test(dt$Inf., g = dt$vars))
  for (p in unique(data_rl$vars)){
    data_glob_relinf <- append(data_glob_relinf, c(frac, p, mean(dt$Inf.[dt$vars==p]),
                                                   median(dt$Inf.[dt$vars==p])))
  }
}
for (p in unique(data_rl$vars)){
  data_glob_relinf <- append(data_glob_relinf, c('all', p, mean(data_rl$Inf.[data_rl$vars==p]),
                                                 median(data_rl$Inf.[data_rl$vars==p])))
}
data_glob_relinf <- matrix(data_glob_relinf, ncol=4, byrow = T)
colnames(data_glob_relinf)<- c('frac', 'param', 'mean', 'median')
write.table(data_glob_relinf, 'data_glob_relinf_dalex.txt')
write.table(data_rl, 'data_relinf_dalex.txt')
pairwise.wilcox.test(data_rl$Inf., g = data_rl$vars)

for (frac in unique(data_rl$rep.fraction..N.)){
  dt <- data_rl[data_rl$rep.fraction..N.==frac,]
  print(pairwise.wilcox.test(dt$Inf., g = dt$vars))
}

x1 <- data_rl$sd[data_rl$vars=='T']
x2 <- data_rl$sd[data_rl$vars=='Sal']
x3 <- data_rl$sd[data_rl$vars=='Si']
x4 <- data_rl$sd[data_rl$vars=='NO3']
x5 <- data_rl$sd[data_rl$vars=='Phos']
x6 <- data_rl$sd[data_rl$vars=='Fe']
# x7 <- data_rl$sd[data_rl$vars=='Chla']
x8 <- data_rl$sd[data_rl$vars=='SI_NO3']
colores = brewer.pal(N, "Dark2")
pdf(family="Helvetica",'relative_influence_violin_sd_dalex.pdf')
vioplot(x1, x2, x3,x4,x5,x6,x8, names=c('T', 'Sal', 'Si', 'NO3', 'Phos', 'Fe', 'SI_NO3'), col =colores)
dev.off()

data_frc <- rep(list(NA), N)
count=1
for (fr in unique(data_rl$rep.fraction..N.)){
  data_frc[[count]]=data_rl$sd[data_rl$rep.fraction..7.==fr]
  count=count+1
}

#pdf(family="Helvetica",'relative_influence_violin_sdfrac_dalex.pdf')
#vioplot(data_frc[[1]], data_frc[[2]],data_frc[[3]],data_frc[[4]],data_frc[[5]],data_frc[[6]],
#        names=unique(data_rl$rep.fraction..7.), col =colores)
#dev.off()

