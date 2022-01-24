#!bin/usr/bin/env Rscript
library("gbm")
library('randomForest')
library('mgcv')
library('nnet')
library("dismo")
library("FactoMineR")
# library("factoextra")
#library("readxl")
# library("ggplot2")
# library("reshape2")
library("gplots")
# library("plotly")
library("stringr")
# library("caret")
library('mapproj')
library('mapplots')
library('maptools')
library('SDMTools')
library('RColorBrewer')
library('ncdf4')
library("CDFt")
library('plotrix')
library('png')
library('extrafont')
library('grid')
library('matlab')

source('axis_map0.R')
source('hide_arctic.R')
cambria <- readRDS('cambria.rds')
model0<-'model-mean'
nvar <- 7
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
best_models <- read.table('best_selected_models.txt', header = T)
best_models <- as.data.frame(best_models)
N <- dim(best_models)[1]
black_caspien <- expand.grid(seq(28.5, 60.5, 1),seq(37.5, 48.5, 1))
med_sea <- expand.grid(seq(-2.5, 35.5, 1),seq(30.5, 45.5, 1))
closed_sea <- rbind(black_caspien, med_sea)
df <- readRDS('Genocenoses_env_parameters_woa.rds')
f2006 = ncdf4::nc_open(paste('all_7_mon_',model0,'_rcp85_200601_201512_clim.nc',sep=''))
f2090 = ncdf4::nc_open(paste('all_7_mon_',model0,'_rcp85_209001_209912_clim.nc',sep=''))
i =1
lati<-seq(-89.5, 89.5, 1)
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
lati_sorted<-seq(-89.5, 60, 1)
longi_sorted =sort(longi)
variables <- c(6:11, 13)

data_p <- readRDS('data_p.rds')
data_f <- readRDS('data_f.rds')
data_w <- readRDS('data_w.rds')
data <- readRDS('model-mean_2006-15.rds')
data1 <- readRDS('model-mean_2090-99.rds')
data_p_scaled <- readRDS('data_p_scaled.rds')
data_f_scaled <- readRDS('data_f_scaled.rds')
# One parameter varies 
drivers_al_val <- list()
for (i in 1:nvar){
  data_d <- NULL
  for (j in 1:nvar){
    if (j==i){
      data_d <- append(data_d, as.vector(data_f_scaled[[j]]))
    } else{
      data_d <- append(data_d, as.vector(data_p_scaled[[j]]))
    }
  }
  data_d <- matrix(data_d, ncol = nvar, byrow = F)
  drivers_al_val[[i]]<-data_d[!is.na(apply(data_d, 1, sum)),]
}

# delta of the parameters
deltas<- NULL
for (i in 1:nvar){
  pres <- as.vector(data_p[[i]])
  fut <- as.vector(data_f[[i]])
  delta <- fut-pres
  deltas <- cbind(deltas, delta)
}
deltas0 <- deltas[!is.na(apply(deltas, 1, sum)),]
saveRDS(deltas, 'deltas.rds')
longitude <- readRDS('longitude.rds')
latitude <- readRDS('latitude.rds')
longitude0 <- as.vector(longitude)[!is.na(apply(deltas, 1, sum))]
latitude0 <- as.vector(latitude)[!is.na(apply(deltas, 1, sum))]
cells <- paste(latitude0, longitude0, sep='_')
bray_curtis <- readRDS('model-mean_bray_curtis.rds')
bray_curtis_f_eez<- readRDS('model-mean_bray_curtis_cond_fish_eez.rds')
pred_2006_list<-readRDS(paste(model0, '_pred_2006_list.rds', sep=''))
pred_2090_list<-readRDS(paste(model0, '_pred_2090_list.rds', sep=''))
pred_2090_list_noT<-readRDS(paste(model0, '_pred_2090_list_noT.rds', sep=''))
mapping <- match(cells, bray_curtis$cell)
rev_mapping <- match(bray_curtis$cell,cells )
saveRDS(mapping,'mapping_lo_lt.rds')
saveRDS(rev_mapping,'rev_mapping_lo_lt.rds')
pred_2006_list <- pred_2006_list[mapping,]
pred_2090_list <- pred_2090_list[mapping,]
pred_2090_list_noT <- pred_2090_list_noT[mapping,]
bray_curtis <- bray_curtis[mapping,]
bray_curtis_f_eez <- bray_curtis_f_eez[mapping,]
models_bt <- readRDS('models_bt.rds')
models_gam <- readRDS('models_gam.rds')
models_nn <- readRDS('models_nn.rds')
models_rf <- readRDS('models_rf.rds')

dom_com_06 <- readRDS('dominant_communities_2006.rds')
dom_com_90 <- readRDS('dominant_communities_2090.rds')
dom_com_90_noT <- readRDS('dominant_communities_2090_noT.rds')
weights_cos <- readRDS('model-mean_weights_cos.rds')
number_of_changes <- readRDS('number_of_changes.rds')
number_of_changes <- number_of_changes[mapping]
dom_com_06 <- dom_com_06[mapping,]
dom_com_90 <- dom_com_90[mapping,]
dom_com_90_noT <- dom_com_90_noT[mapping,]
weights_cos <- weights_cos[mapping]
fractions <- c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
color_drivers <- brewer.pal(length(variables), "Dark2")
pred_new06 <- NULL
pred_new90 <- NULL
pred_new90_noT <- NULL
for (u in 1:N){
  com <- strsplit(colnames(pred_2006_list)[u], '_')[[1]][2]
  frac <- which(fractions==strsplit(colnames(pred_2006_list)[u], '_')[[1]][1])
  sel_06 <- dom_com_06[,frac*2-1]==com
  sel_90 <- dom_com_90[,frac*2-1]==com
  sel_90_noT <- dom_com_90_noT[,frac*2-1]==com
  vec_2006 <- rep(0, dim(pred_2006_list)[1])
  vec_2006[sel_06] <- pred_2006_list[sel_06, u]
  vec_2090 <- rep(0, dim(pred_2090_list)[1])
  vec_2090[sel_90] <- pred_2090_list[sel_90, u]
  vec_2090_noT <- rep(0, dim(pred_2090_list_noT)[1])
  vec_2090_noT[sel_90_noT] <- pred_2090_list_noT[sel_90_noT, u]
  pred_new06 <-append(pred_new06, vec_2006)
  pred_new90 <-append(pred_new90, vec_2090)
  pred_new90_noT <-append(pred_new90_noT, vec_2090_noT)
}
pred_new06 <- matrix(pred_new06, ncol=N)
pred_new90 <- matrix(pred_new90, ncol=N)
pred_new90_noT <- matrix(pred_new90_noT, ncol=N)
saveRDS(pred_new06, 'pred_dom06.rds')
saveRDS(pred_new90, 'pred_dom90.rds')
saveRDS(pred_new90_noT, 'pred_dom90_noT.rds')
pred_new06 <- readRDS('pred_dom06.rds')
pred_new90 <- readRDS('pred_dom90.rds')
pred_new90_noT <- readRDS('pred_dom90_noT.rds')

bc_new <- apply(abs(pred_new06-pred_new90), 1, sum)/apply(pred_new06+pred_new90, 1, sum)
saveRDS(bc_new, 'model-mean_bray-curtis_dom.rds')

bc_new_noT <- apply(abs(pred_new06-pred_new90_noT), 1, sum)/apply(pred_new06+pred_new90_noT, 1, sum)
saveRDS(bc_new_noT, 'model-mean_bray-curtis_dom_noT.rds')


values <- 100*(bc_new-min(bc_new))/((max(bc_new)-min(bc_new)))+1
set = jet.colors(101)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_high_fish_dom_2006-2090.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
selection <- bray_curtis_f_eez$cond_fish
selection[is.na(selection)]<-FALSE
points(x=bray_curtis$Long[selection], y =bray_curtis$Lat[selection],
       col=set[values[selection]],bg= set[values[selection]],
       pch=15,cex=0.229)
#SDMTools::legend.gradient(pnt,set,limits=c(0,100-round(min(bray_curtis$b_c), digits = 1)),title='Bray-curtis \ndissimilarity index',cex=0.4)
hide_arctic()
axis_map0()
dev.off()

values <- 100*(bc_new-min(bc_new))/((max(bc_new)-min(bc_new)))+1
set = jet.colors(101)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_eez_dom_2006-2090.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
cond_eez <- bray_curtis_f_eez$cond_eez
selection <- cond_eez
points(x=bray_curtis$Long[selection], y =bray_curtis$Lat[selection],
       col=set[values[selection]],bg= set[values[selection]],
       pch=15,cex=0.229)
#SDMTools::legend.gradient(pnt,set,limits=c(0,100-round(min(bray_curtis$b_c), digits = 1)),title='Bray-curtis \ndissimilarity index',cex=0.4)
hide_arctic()
axis_map0()
dev.off()

drivers_solo_niche <- vector("list", N)
drivers_solo_niche_glob <- vector("list", N)
drivers_glob <-matrix(0, ncol=length(variables), nrow = length(cells))
selections <- vector("list", N)
local(
for (u in 1:N){
  com <- strsplit(colnames(pred_2006_list)[u], '_')[[1]][2]
  frac <- which(fractions==strsplit(colnames(pred_2006_list)[u], '_')[[1]][1])
  bt_model <- models_bt[[u]]
  gam_model <- models_gam[[u]]
  nn_model <- models_nn[[u]]
  rf_model <- models_rf[[u]]
  selection <- abs(pred_2006_list[,u]-pred_2090_list[,u])>0 & (dom_com_06[,frac*2-1]==com | dom_com_90[,frac*2-1]==com) & bc_new > 1/6
  selections[[u]]<<- selection
  if (sum(selection, na.rm=T)>0){
    lat_selection <- latitude0[selection]
    long_selection <- longitude0[selection]
    predictions_drivers <- NULL
    for (i in 1:length(variables)){
      newdata <- as.data.frame(drivers_al_val[[i]][selection,])
      colnames(newdata)<-colnames(df)[variables]
      if (!is.na(bt_model)){
        pred_bt <- gbm::predict.gbm(bt_model,newdata , n.trees=bt_model$gbm.call$best.trees, type="response")
      }
      pred_rf <- stats::predict(rf_model, newdata, type='prob')[,2]
      pred_nn <- stats::predict(nn_model, newdata, type='raw')[,2]
      if (!is.na(gam_model)){
        pred_gam <- mgcv::predict.gam(gam_model, newdata, type='response')
      }
      if (!is.na(bt_model) & !is.na(gam_model)){
        pred <- (pred_gam+pred_rf+pred_bt+pred_nn)/4
      } else if (is.na(bt_model) & !is.na(gam_model)){
        pred <- (pred_gam+pred_rf+pred_nn)/3
      } else if (!is.na(bt_model) & is.na(gam_model)){
        pred <- (pred_rf+pred_bt+pred_nn)/3
      }
      predictions_drivers <- append(predictions_drivers, pred)
    }
    predictions_drivers <- matrix(predictions_drivers, ncol=length(variables))
    sum_driv <- rep(0, sum(selection))
    for (e in 1:length(variables)){
      sum_driv <- sum_driv+abs(predictions_drivers[,e]-pred_2006_list[selection,u])
    }
    for (e in 1:length(variables)){
       drivers_solo_niche[[u]] <<- append(drivers_solo_niche[[u]], abs(predictions_drivers[,e]-pred_2006_list[selection,u])/sum_driv )
    }
    drivers_solo_niche[[u]] <<- matrix(drivers_solo_niche[[u]], ncol=length(variables))
    drivers_solo_niche[[u]][is.na(apply(drivers_solo_niche[[u]] ,1,sum)),] <<- 0
    #drivers_solo_niche[[u]] <- drivers_solo_niche[[u]][!is.na(apply(drivers_solo_niche[[u]] ,1,sum)),]
    for (e in 1:length(variables)){
    	drivers_solo_niche_glob[[u]] <<- append(drivers_solo_niche_glob[[u]], 
                                              sum(drivers_solo_niche[[u]][,e]*abs(pred_2006_list[selection,u]-pred_2090_list[selection,u])*cos(lat_selection*2*pi/360)))
        #drivers_solo_niche_glob[[u]] <<- append(drivers_solo_niche_glob[[u]], 
         #                                       sum(drivers_solo_niche[[u]][,e]*cos(lat_selection*2*pi/360)))
    }
    drivers_glob[selection,]<<- drivers_glob[selection,]+drivers_solo_niche[[u]]*sum_driv
  }
})
#rm(models_bt, models_gam, models_nn, models_rf)
saveRDS(drivers_solo_niche, 'drivers_solo_niche.rds')
saveRDS(drivers_glob, 'drivers_glob.rds')
saveRDS(drivers_solo_niche_glob, 'drivers_solo_niche_glob.rds')
saveRDS(selections, 'selections.rds')

drivers_solo_niche <- readRDS('drivers_solo_niche.rds')
selections <- readRDS('selections.rds')
for (u in 1:N){
  col_vec <- as.vector(apply(drivers_solo_niche[[u]], 1, which.max))
  alpha_vec <- apply(drivers_solo_niche[[u]], 1, max)
  selection <- selections[[u]]
  #     pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[u], '_drivers_alpha_dom.pdf', sep='')
  #         ,width=10,height=4.065)
  #     maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #     points(longitude0[selection], latitude0[selection], col=scales::alpha(color_drivers[col_vec], alpha = alpha_vec), pch=15,cex=0.229)
  #     legend(x = 60, y = 60, legend = colnames(df)[6:13], col=color_drivers,pch=15,cex=8, ncol=2)
  #     dev.off()
  data_contour<- rep(list(matrix(NA, ncol=150, nrow=360)),length(variables))
  alpha_contour <- matrix(NA, ncol=150, nrow=360)
  data_contour1 <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:length(col_vec)){
    dd=col_vec[i]
    data_contour[[dd]][which(longi_sorted==longitude0[selection][i]),which(lati==latitude0[selection][i])]= dd
    alpha_contour[which(longi_sorted==longitude0[selection][i]),which(lati==latitude0[selection][i])]= round(alpha_vec[i]*100)+1
    data_contour1[which(longi_sorted==longitude0[selection][i]),which(lati==latitude0[selection][i])]= 1
  }
  for (i in 1:length(variables)){
    data_contour[[i]][is.na(data_contour[[i]])]<-0
  }
  
  pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[u], '_drivers_dom_contour.pdf', sep=''), width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  set_col=color_drivers
  ec=scales::alpha('white',0)
  for (i in 1:length(variables)){
    .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
  }
  contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
  data_contour1[is.na(data_contour1)]<-0
  png("contour1.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  set_col=color_drivers
  ec=scales::alpha('white',0)
  for (i in 1:length(variables)){
    .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
  }
  hide_arctic()
  axis_map0()
  # contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
  dev.off()
  set_al <- NULL
  for (i in seq(0,1,0.01)){
    set_al <-append(set_al, rgb(i,0,0))
  }
  png('alpha1.png', width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
  .filled.contour(longi_sorted, lati_sorted, alpha_contour,
                  col=set_al, levels = 1:101)
  hide_arctic()
  #axis_map0()
  dev.off()
  
  im <- readPNG('contour1.png')
  alph <- readPNG('alpha1.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[u], '_drivers_alpha_dom_contour.pdf', sep=''), width=10,height=4.065)
  grid.raster(w)
  dev.off()
}


drivers_solo_niche_glob <- readRDS('drivers_solo_niche_glob.rds')
fractions0 <- readRDS('fractions0.rds')
contributions_var <- matrix(0, ncol = length(variables), nrow=6)
for (i in 1:6){
  frac <- fractions0[[i]]
  for (k in frac){
    if (!is.na(drivers_solo_niche_glob[[k]])){
      contributions_var[i,] <- contributions_var[i,]+drivers_solo_niche_glob[[k]]
    }
  }
}
sum_row <- apply(contributions_var, 1, sum)
for (i in 1:6){
  contributions_var[i,]<-contributions_var[i,]/sum_row[i] 
}
row.names(contributions_var)<-fractions
pdf(family="Helvetica",'drivers_fraction.pdf')
barplot(t(contributions_var), col = brewer.pal(length(variables), 'Dark2'), 
        legend = colnames(df)[variables], main='',
        xlim=c(0, 10), cex.axis = 1.3, cex.lab =1.3, cex.names =0.79)
dev.off()
contrib_drivers <- t(contributions_var)
contrib_drivers <- cbind(contrib_drivers, apply(t(contributions_var),1,mean))
write.table(contrib_drivers, 'contributions_drivers_fraction.txt')

drivers_solo_niche_noT <- vector("list", N)
drivers_solo_niche_glob_noT <- vector("list", N)
selections_noT <- vector("list", N)
drivers_glob_noT <-matrix(0, ncol=length(variables)-1, nrow = length(cells))
local(
  for (u in 1:N){
    com <- strsplit(colnames(pred_2006_list)[u], '_')[[1]][2]
    frac <- which(fractions==strsplit(colnames(pred_2006_list)[u], '_')[[1]][1])
    bt_model <- models_bt[[u]]
    gam_model <- models_gam[[u]]
    nn_model <- models_nn[[u]]
    rf_model <- models_rf[[u]]
    selection <- abs(pred_2006_list[,u]-pred_2090_list_noT[,u])>0 & (dom_com_06[,frac*2-1]==com | dom_com_90_noT[,frac*2-1]==com) & bc_new_noT > 1/6
    selections_noT[[u]] <<- selection
    if (sum(selection, na.rm=T)>0){
      lat_selection <- latitude0[selection]
      long_selection <- longitude0[selection]
      predictions_drivers <- NULL
      for (i in 2:length(variables)){
        newdata <- as.data.frame(drivers_al_val[[i]][selection,])
        colnames(newdata)<-colnames(df)[variables]
        if (!is.na(bt_model)){
          pred_bt <- gbm::predict.gbm(bt_model,newdata , n.trees=bt_model$gbm.call$best.trees, type="response")
        }
        pred_rf <- stats::predict(rf_model, newdata, type='prob')[,2]
        pred_nn <- stats::predict(nn_model, newdata, type='raw')[,2]
        if (!is.na(gam_model)){
          pred_gam <- mgcv::predict.gam(gam_model, newdata, type='response')
        }
        if (!is.na(bt_model) & !is.na(gam_model)){
          pred <- (pred_gam+pred_rf+pred_bt+pred_nn)/4
        } else if (is.na(bt_model) & !is.na(gam_model)){
          pred <- (pred_gam+pred_rf+pred_nn)/3
        } else if (!is.na(bt_model) & is.na(gam_model)){
          pred <- (pred_rf+pred_bt+pred_nn)/3
        }
        predictions_drivers <- append(predictions_drivers, pred)
      }
      predictions_drivers <- matrix(predictions_drivers, ncol=length(variables)-1)
      sum_driv <- rep(0, sum(selection))
      for (e in 2:length(variables)){
        sum_driv <- sum_driv+abs(predictions_drivers[,e-1]-pred_2006_list[selection,u])
      }
      for (e in 2:length(variables)){
        drivers_solo_niche_noT[[u]] <<- append(drivers_solo_niche_noT[[u]], abs(predictions_drivers[,e-1]-pred_2006_list[selection,u])/sum_driv )
      }
      drivers_solo_niche_noT[[u]] <<- matrix(drivers_solo_niche_noT[[u]], ncol=length(variables)-1)
      drivers_solo_niche_noT[[u]][is.na(apply(drivers_solo_niche_noT[[u]] ,1,sum)),] <<- 0
      #drivers_solo_niche_noT[[u]] <<- drivers_solo_niche_noT[[u]][!is.na(apply(drivers_solo_niche_noT[[u]] ,1,sum)),]
      for (e in 2:length(variables)){
         drivers_solo_niche_glob_noT[[u]] <<- append(drivers_solo_niche_glob_noT[[u]], 
                                                     sum(drivers_solo_niche_noT[[u]][,e-1]*abs(pred_2006_list[selection,u]-pred_2090_list_noT[selection,u])*cos(lat_selection*2*pi/360)))
        #drivers_solo_niche_glob_noT[[u]] <<- append(drivers_solo_niche_glob_noT[[u]], 
         #                                           sum(drivers_solo_niche_noT[[u]][,e-1]*cos(lat_selection*2*pi/360)))
      }
      drivers_glob_noT[selection,]<<- drivers_glob_noT[selection,]+drivers_solo_niche_noT[[u]]*sum_driv
    }
})
rm(models_bt, models_gam, models_nn, models_rf)
saveRDS(drivers_solo_niche_noT, 'drivers_solo_niche_noT.rds')
saveRDS(drivers_glob_noT, 'drivers_glob_noT.rds')
saveRDS(drivers_solo_niche_glob_noT, 'drivers_solo_niche_glob_noT.rds')
saveRDS(selections_noT, 'selections_noT')

drivers_solo_niche_noT <- readRDS('drivers_solo_niche_noT.rds')
selections_noT <- readRDS('selections_noT')
for (u in 1:N){
  col_vec <- as.vector(apply(drivers_solo_niche_noT[[u]], 1, which.max))+1
  alpha_vec <- apply(drivers_solo_niche_noT[[u]], 1, max)
  selection <- selections_noT[[u]]
  #     pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[u], '_drivers_alpha_dom.pdf', sep='')
  #         ,width=10,height=4.065)
  #     maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  #     points(longitude0[selection], latitude0[selection], col=scales::alpha(color_drivers[col_vec], alpha = alpha_vec), pch=15,cex=0.229)
  #     legend(x = 60, y = 60, legend = colnames(df)[6:13], col=color_drivers,pch=15,cex=7, ncol=2)
  #     dev.off()
  data_contour<- rep(list(matrix(NA, ncol=150, nrow=360)),length(variables))
  alpha_contour <- matrix(NA, ncol=150, nrow=360)
  data_contour1 <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:length(col_vec)){
    dd=col_vec[i]
    data_contour[[dd]][which(longi_sorted==longitude0[selection][i]),which(lati==latitude0[selection][i])]= dd
    alpha_contour[which(longi_sorted==longitude0[selection][i]),which(lati==latitude0[selection][i])]= round(alpha_vec[i]*100)+1
    data_contour1[which(longi_sorted==longitude0[selection][i]),which(lati==latitude0[selection][i])]= 1
  }
  for (i in 1:length(variables)){
    data_contour[[i]][is.na(data_contour[[i]])]<-0
  }
  
  pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[u], '_drivers_dom_contour_noT.pdf', sep=''), width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  set_col=color_drivers
  ec=scales::alpha('white',0)
  for (i in 1:length(variables)){
    .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
  }
  contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
#   data_contour1[is.na(data_contour1)]<-0
#   png("contour1.png", width=10,height=4.065, res=200, units = 'in')
#   maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
#   set_col=color_drivers
#   ec=scales::alpha('white',0)
#   for (i in 1:8){
#     .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
#   }
#   # contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
#   dev.off()
#   set_al <- NULL
#   for (i in seq(0,1,0.01)){
#     set_al <-append(set_al, rgb(i,0,0))
#   }
#   png('alpha1.png', width=10,height=4.065, res=200, units = 'in')
#   maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
#   .filled.contour(longi_sorted, lati_sorted, alpha_contour,
#                   col=set_al, levels = 1:101)
#   dev.off()
#   
#   im <- readPNG('contour1.png')
#   alph <- readPNG('alpha1.png')
#   w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
#   pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[u], '_drivers_alpha_dom_contour_noT.pdf', sep=''), width=10,height=4.065)
#   grid.raster(w)
#   dev.off()
}


drivers_solo_niche_glob_noT <- readRDS('drivers_solo_niche_glob_noT.rds')
fractions0 <- readRDS('fractions0.rds')
contributions_var_noT <- matrix(0, ncol = length(variables)-1, nrow=6)
for (i in 1:6){
  frac <- fractions0[[i]]
  for (k in frac){
    if (!is.na(drivers_solo_niche_glob_noT[[k]])){
      contributions_var_noT[i,] <- contributions_var_noT[i,]+drivers_solo_niche_glob_noT[[k]]
    }
  }
}
sum_row_noT <- apply(contributions_var_noT, 1, sum)
for (i in 1:6){
  contributions_var_noT[i,]<-contributions_var_noT[i,]/sum_row_noT[i] 
}
row.names(contributions_var_noT)<-fractions
pdf(family="Helvetica",'drivers_fraction_noT.pdf')
barplot(t(contributions_var_noT), col = brewer.pal(length(variables), 'Dark2')[2:length(variables)], 
        legend = colnames(df)[variables[2:length(variables)]], main='',
        xlim=c(0, 10), cex.axis = 1.3, cex.lab =1.3, cex.names =0.79)
dev.off()
contrib_drivers_noT <- t(contributions_var_noT)
contrib_drivers_noT <- cbind(contrib_drivers_noT, apply(t(contributions_var_noT),1,mean))
write.table(contrib_drivers_noT, 'contributions_drivers_fraction_noT.txt')


bc_new <-readRDS('model-mean_bray-curtis_dom.rds')
values <- 100*(bc_new-min(bc_new))/((max(bc_new)-min(bc_new)))+1
#values <- 100*(bc_new-min(bc_new))/((max(bc_new)-min(bc_new)))+1
data_contour <- matrix(NA, ncol=150, nrow=360)
data_contour1 <- matrix(NA, ncol=150, nrow=360)
for (i in 1:dim(bray_curtis)[1]){
  data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= round(values[i])
  #data_contour1[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= round(values[i], -1)/10+1
}
data_contour1 <- round(data_contour, -1)/10+1
set = jet.colors(101)
set1 = jet.colors(10)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_dom_2006-2090.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
# points(x=bray_curtis$Long, y =bray_curtis$Lat, col=set[values],bg= set[values],pch=15,cex=0.229)
SDMTools::legend.gradient(pnt,set1,limits=c(0,1),title='Bray-curtis \ndissimilarity index',cex=0.4)
.filled.contour(x=longi_sorted,y=lati_sorted, z=data_contour1,  col=set1, levels = c(1:11))
hide_arctic()
axis_map0()
dev.off()


bc_new_noT <-readRDS('model-mean_bray-curtis_dom_noT.rds')
values_noT <- 100*(bc_new_noT-min(bc_new_noT))/((max(bc_new_noT)-min(bc_new_noT)))+1
data_contour <- matrix(NA, ncol=150, nrow=360)
for (i in 1:dim(bray_curtis)[1]){
  data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati_sorted==bray_curtis$Lat[i])]= round(values_noT[i])
}
set = jet.colors(101)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_dom_2006-2090_noT.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
# points(x=bray_curtis$Long, y =bray_curtis$Lat, col=set[values_noT],bg= set[values_noT],pch=15,cex=0.229)
SDMTools::legend.gradient(pnt,set1,limits=c(0,1-round(min(bc_new_noT), digits = 1)),title='Bray-curtis \ndissimilarity index',cex=0.4)
.filled.contour(x=longi_sorted,y=lati_sorted, z=data_contour,  col=set, levels = c(1:101))
hide_arctic()
axis_map0()
dev.off()

cond_trop <- bray_curtis$Lat < 25 & bray_curtis$Lat > -25
cond_n_temp <- bray_curtis$Lat >25
cond_s_temp <- bray_curtis$Lat > -60 & bray_curtis$Lat < -25
cond_polar <- bray_curtis$Lat < -60
cond_all <- bray_curtis$Lat > -100
cond_climate <- list(cond_trop, cond_n_temp, cond_s_temp, cond_polar, cond_all)
climates <- c('Tropical', 'North_temperate', 'South_temperate', 'Polar', 'All')

cond_na <- bray_curtis$Long<30 & bray_curtis$Long > -80 & bray_curtis$Lat >0 
cond_sa <- bray_curtis$Long<20 & bray_curtis$Long > -70 & bray_curtis$Lat < 0
cond_np <- (bray_curtis$Long>120 | bray_curtis$Long < -90) & bray_curtis$Lat >0
cond_sp <- (bray_curtis$Long>120 | bray_curtis$Long < -70) & bray_curtis$Lat <0
cond_i <- bray_curtis$Long>20 & bray_curtis$Long < 120 & bray_curtis$Lat <30
basins_conds <- list(cond_na, cond_sa, cond_np, cond_sp, cond_i)
basins <- c('North_A', 'South_A', 'North_P', 'South_P', 'Indian')

glob_stats_fish_eez_basins <- NULL
cond_eez<- bray_curtis_f_eez$cond_eez
cond_fish<- bray_curtis_f_eez$cond_fish
c<-1
for (cond in basins_conds){
  mean_bc <- sum((bc_new[cond]*100)*weights_cos[cond])/sum(weights_cos[cond])
  mean_bc_noT <- sum((bc_new_noT[cond]*100)*weights_cos[cond])/sum(weights_cos[cond])
  mean_bc_eez <- sum((bc_new[cond & cond_eez])*weights_cos[cond & cond_eez])/sum(weights_cos[cond & cond_eez])
  mean_bc_fish <- sum((bc_new[cond & cond_fish])*weights_cos[cond & cond_fish], na.rm = T)/sum(weights_cos[cond & cond_fish], na.rm=T)
  glob_stats_fish_eez_basins <- rbind(glob_stats_fish_eez_basins, c(basins[c], mean_bc,mean_bc_noT,mean_bc_eez, mean_bc_fish))
  c<-c+1
}
colnames(glob_stats_fish_eez_basins)<- c('basins', 'mean_bc','mean_bc_noT','eez','fish' )
write.table(glob_stats_fish_eez_basins, 'glob_stats_fish_eez_basins_bc_dom.txt', col.names = T)

glob_stats_fish_eez_climate <- NULL
c<-1
for (cond in cond_climate){
  mean_bc <- sum((bc_new[cond]*100)*weights_cos[cond])/sum(weights_cos[cond])
  mean_bc_noT <- sum((bc_new_noT[cond]*100)*weights_cos[cond])/sum(weights_cos[cond])
  mean_bc_eez <- sum((bc_new[cond & cond_eez])*weights_cos[cond & cond_eez])/sum(weights_cos[cond & cond_eez])
  mean_bc_fish <- sum((bc_new[cond & cond_fish])*weights_cos[cond & cond_fish], na.rm = T)/sum(weights_cos[cond & cond_fish], na.rm=T)
  glob_stats_fish_eez_climate <- rbind(glob_stats_fish_eez_climate, c(climates[c], mean_bc,mean_bc_noT, mean_bc_eez, mean_bc_fish))
  c<-c+1
}
colnames(glob_stats_fish_eez_climate)<- c('climate', 'mean_bc','mean_bc_noT','eez','fish' )
write.table(glob_stats_fish_eez_climate, 'glob_stats_fish_eez_climate_bc_dom.txt', col.names = T)

drivers_glob <- readRDS('drivers_glob.rds')
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
drivers_glob0 <- t(apply(drivers_glob, 1, function(x){x/sum(x)}))*100
drivers_glob0 <- cbind(drivers_glob0, latitude0)
drivers_glob0 <- cbind(drivers_glob0, longitude0)
drivers_glob1 <- drivers_glob0[!is.na(apply(drivers_glob0,1,sum)),]
colors2 <- colorRampPalette(c('lightblue','lightgreen', 'orange'))(101)
alpha_col <- bc_new[!is.na(apply(drivers_glob0,1,sum))]
count <- 1
for (i in variables){
  pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_',colnames(df)[i],'_alone_dom.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  points(drivers_glob1[,9], drivers_glob1[,8], 
         col=scales::alpha(colors2[round(drivers_glob1[,count])+1], alpha_col), pch=15,cex=0.229)
  SDMTools::legend.gradient(pnt,colors2,limits=c(0,100),title=paste(colnames(df)[i],'\nDriver impact (%)', sep=''),cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  data_contour<- matrix(NA, ncol=150, nrow=360)
  alpha_contour <- matrix(NA, ncol=150, nrow=360)
  for (j in 1:dim(drivers_glob1)[1]){
    data_contour[which(longi_sorted==drivers_glob1[j,9]),which(lati_sorted==drivers_glob1[j,8])]= round(drivers_glob1[j,count])+1
    alpha_contour[which(longi_sorted==drivers_glob1[j,9]),which(lati_sorted==drivers_glob1[j,8])]= round(alpha_col[j]*100)+1
  }
  png("contour1.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=colors2, levels=c(1:101))
  hide_arctic()
  axis_map0()
  dev.off()
  set_al <- NULL
  for (u in seq(0,1,0.01)){
    set_al <-append(set_al, rgb(u,0,0))
  }
#   png('alpha1.png', width=10,height=4.065, res=200, units = 'in')
#   par(mar=c(0,0,0,0))
#   maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
#   .filled.contour(longi_sorted, lati_sorted, alpha_contour,
#                   col=set_al, levels = 1:101)
#   axis_map0()
#   dev.off()
  
  im <- readPNG('contour1.png')
  #alph <- readPNG('alpha1.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_',colnames(df)[i],'_alone_dom_contour.pdf', sep=''),width=10,height=4.065)
  grid.raster(w)
  dev.off()
  count <- count+1
}

max_driv <- apply(drivers_glob1[,1:length(variables)], 1, which.max)
alpha_driv <- apply(drivers_glob1[,1:length(variables)], 1, max)
color_drivers <- brewer.pal(length(variables), "Dark2")
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_dom_alpha.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
points(drivers_glob1[,9], drivers_glob1[,8], col=scales::alpha(color_drivers[max_driv], alpha_driv/100), pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(df)[variables], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_dom.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
points(drivers_glob1[,9], drivers_glob1[,8], col=color_drivers[max_driv], pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(df)[variables], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

data_contour<- rep(list(matrix(NA, ncol=150, nrow=360)),length(variables))
alpha_contour <- matrix(NA, ncol=150, nrow=360)
data_contour1 <- matrix(NA, ncol=150, nrow=360)
for (i in 1:length(max_driv)){
  dd=max_driv[i]
  data_contour[[dd]][which(longi_sorted==drivers_glob1[i,9]),which(lati==drivers_glob1[i,8])]= dd
  alpha_contour[which(longi_sorted==drivers_glob1[i,9]),which(lati==drivers_glob1[i,8])]= round(alpha_driv[i])+1
  data_contour1[which(longi_sorted==drivers_glob1[i,9]),which(lati==drivers_glob1[i,8])]=1
}
data_contour1[is.na(data_contour1)]<-0
for (i in 1:length(variables)){
  data_contour[[i]][is.na(data_contour[[i]])]<-0
}
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_dom_contour.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
set_col=color_drivers
ec=scales::alpha('white',0)
for (i in 1:length(variables)){
  .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
}
# .filled.contour(x=longi_sorted,y=lati_sorted, z=data_contour,  col=color_drivers, levels = c(1:8))
contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
hide_arctic()
axis_map0()
dev.off()

png("contour1.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
set_col=color_drivers
ec=scales::alpha('white',0)
for (i in 1:length(variables)){
  .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
}
# .filled.contour(longi_sorted, lati_sorted,data_contour,
#                 col=color_drivers, levels=c(1:8))
# contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
hide_arctic()
axis_map0()
dev.off()
set_al <- NULL
for (i in seq(0,1,0.01)){
  set_al <-append(set_al, rgb(i,0,0))
}
png('alpha1.png', width=10,height=4.065, res=200, units = 'in', family="Helvetica")
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
.filled.contour(longi_sorted, lati_sorted, alpha_contour,
                col=set_al, levels = 1:101)
hide_arctic()
# contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='black',add=T, lwd=0.45,drawlabels=F)
#axis_map()
dev.off()

im <- readPNG('contour1.png')
alph <- readPNG('alpha1.png')
w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_dom_contour_alpha.pdf', sep=''),width=10,height=4.065)
grid.raster(w)
dev.off()

contrib_vars_dom <- NULL
for (i in 1:length(variables)){
  varia<-colnames(df)[variables[i]]
  contrib<-sum(cos(drivers_glob1[,8]*2*pi/360)*drivers_glob1[,i]*bc_new[!is.na(apply(drivers_glob0,1,sum))])/sum(cos(drivers_glob1[,8]*2*pi/360)*bc_new[!is.na(apply(drivers_glob0,1,sum))])
  contrib_vars_dom <- rbind(contrib_vars_dom, c(varia, contrib))
}
write.table(contrib_vars_dom, 'contributions_drivers_dom.txt')

drivers_glob2 <- drivers_glob1
max_driv <- apply(drivers_glob2[,1:length(variables)], 1, which.max)
data_contour<- rep(list(matrix(NA, ncol=150, nrow=360)),length(variables))
alpha_contour <- matrix(NA, ncol=150, nrow=360)
data_contour1 <- matrix(NA, ncol=150, nrow=360)
for (i in 1:length(max_driv)){
  dd=max_driv[i]
  data_contour[[dd]][which(longi_sorted==drivers_glob2[i,9]),which(lati==drivers_glob2[i,8])]= dd
  ##alpha_contour[which(longi_sorted==drivers_glob2[i,10]),which(lati==drivers_glob2[i,9])]= round(alpha_driv[i])+1
  data_contour1[which(longi_sorted==drivers_glob2[i,9]),which(lati==drivers_glob2[i,8])]=1
}
data_contour1[is.na(data_contour1)]<-0
for (i in 1:length(variables)){
  data_contour[[i]][is.na(data_contour[[i]])]<-0
}
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_dom_contour_bcfilter.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
set_col=color_drivers
ec=scales::alpha('white',0)
for (i in 1:length(variables)){
  .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
}
hide_arctic()
axis_map0()
#.filled.contour(x=longi_sorted,y=lati_sorted, z=data_contour,  col=color_drivers, levels = c(1:8))
#contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
dev.off()

drivers_glob_noT <- readRDS('drivers_glob_noT.rds')
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
drivers_glob0_noT <- t(apply(drivers_glob_noT, 1, function(x){x/sum(x)}))*100
drivers_glob0_noT <- cbind(drivers_glob0_noT, latitude0)
drivers_glob0_noT <- cbind(drivers_glob0_noT, longitude0)
drivers_glob1_noT <- drivers_glob0_noT[!is.na(apply(drivers_glob0_noT,1,sum)),]
colors2 <- colorRampPalette(c('lightblue','lightgreen', 'orange'))(101)

drivers_glob2_noT <- drivers_glob1_noT
max_driv <- apply(drivers_glob2_noT[,1:length(variables)-1], 1, which.max)+1
data_contour<- rep(list(matrix(NA, ncol=150, nrow=360)),length(variables))
alpha_contour <- matrix(NA, ncol=150, nrow=360)
data_contour1 <- matrix(NA, ncol=150, nrow=360)
for (i in 1:length(max_driv)){
  dd=max_driv[i]
  data_contour[[dd]][which(longi_sorted==drivers_glob2_noT[i,8]),which(lati==drivers_glob2_noT[i,7])]= dd
  ##alpha_contour[which(longi_sorted==drivers_glob2[i,10]),which(lati==drivers_glob2[i,9])]= round(alpha_driv[i])+1
  data_contour1[which(longi_sorted==drivers_glob2_noT[i,8]),which(lati==drivers_glob2_noT[i,7])]=1
}
data_contour1[is.na(data_contour1)]<-0
for (i in 2:length(variables)){
  data_contour[[i]][is.na(data_contour[[i]])]<-0
}
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_dom_contour_bcfilter_noT.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
set_col=color_drivers
ec=scales::alpha('white',0)
for (i in 2:length(variables)){
  .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
}
# .filled.contour(x=longi_sorted,y=lati_sorted, z=data_contour,  col=color_drivers, levels = c(1:length(variables)))
#contour(x=longi_sorted,y=lati_sorted, z=data_contour1,col='white',add=T, lwd=0.45,drawlabels=F)
hide_arctic()
axis_map0()
dev.off()
# correlation_maps

# mapping0 <- match(bray_curtis$cell, cells)
# com <- which(colnames(pred_2006_list)=='180-2000_5')
# for (u in 1:8){
#   lati_square <- NULL
#   longi_square <- NULL
#   cor_loc<- NULL
#   for (i in seq(min(bray_curtis$Lat)+3, max(bray_curtis$Lat)-3, 1)){
#     for (j in seq(min(bray_curtis$Long)+3, max(bray_curtis$Long)-3, 1)){
#       square0 <- which(bray_curtis$Lat>=i-3 & bray_curtis$Lat<=i+3 & bray_curtis$Long>=j-3 & bray_curtis$Long<=j+3)
#       #square <- mapping[square0]
#       if (length(square0)>20){
#         selec <- which((pred_2090_list[square0,com] - pred_2006_list[square0,com])>0.2 | (pred_2090_list[square0,com] - pred_2006_list[square0,com])< -0.2)
#         if (length(selec)>20){
#           cor_local <- cor(pred_2090_list[square0,com]-pred_2006_list[square0,com], deltas0[square0,u], method = 'spearman')*100+100
#           lati_square <- append(lati_square, i)
#           longi_square <- append(longi_square, j)
#           cor_loc <- append(cor_loc, cor_local)
#         }
#       }
#     }
#   }
#   pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
#   col_scale <- redgreen(201)
#   pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[com], '_', colnames(df)[u+5], '.pdf', sep='')
#       ,width=10,height=4.065)
#   maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
#   points(longi_square, lati_square, col=col_scale[round(cor_loc)+1], pch=15,cex=0.235)
#   SDMTools::legend.gradient(pnt,col_scale,limits=c(-1, 1),title='Local correlation',cex=0.4)
#   dev.off()
# 
#   
#   cell_selec <- paste(lati_square,longi_square, sep='_')
#   set1 = colorRampPalette(c("darkgreen","light blue", "dark violet"))(201)
#   deltas_p <- (pred_2090_list[match(cell_selec, cells),com]-pred_2006_list[match(cell_selec, cells),com])*100+100
#   pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[com], '_', colnames(df)[u+5],'_delta_p.pdf', sep='')
#       ,width=10,height=4.065)
#   maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
#   points(longi_square, lati_square, col=set1[round(deltas_p)+1], pch=15,cex=0.235)
#   SDMTools::legend.gradient(pnt,set1,limits=c(-1, 1),title='Delta P',cex=0.4)
#   dev.off()
#   
#   deltas_param <- deltas0[,u][match(cell_selec, cells)]
#   max_ok <- max(abs(max(deltas_param, na.rm=T)), abs(min(deltas_param, na.rm=T)))
#   max_ok <- 7.224227
#   values_d <- round(seq(-max_ok,max_ok,2*max_ok/99), digits = 8)
#   col_index <- NULL
#   for (j in 1:length(deltas_param)){
#     if (!is.na(deltas_param[j])){
#       col_index <- append(col_index, which.min(abs(values_d-deltas_param[j])))
#     } else{
#       col_index <- append(col_index, NA)
#     }
#   }
#   colfunc <- colorRampPalette(c('blue', 'lightblue','white', 'orange', 'red'))(100)
#   pdf(family="Helvetica",file=paste(colnames(pred_2006_list)[com], '_', colnames(df)[u+5],'_delta_T.pdf', sep='')
#       ,width=10,height=4.065)
#   maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
#   points(longitude0[match(cell_selec, cells)], latitude0[match(cell_selec, cells)], col=colfunc[col_index],bg=colfunc[col_index], pch=15,cex=0.235)
#   SDMTools::legend.gradient(pnt,colfunc,limits=c(-max_ok, max_ok),title='Delta T',cex=0.4)
#   dev.off()
#   
# }
# 
# cors <- NULL
# for (u in 1:30){
#   for (j in 1:7){
#     selec <- pred_2006_list[,u]>0.5 | pred_2090_list[,u]>0.5
#     delta_pred <- pred_2090_list[selec,u]-pred_2006_list[selec,u]
#     selec0 <- mapping0[which(selec)]
#     cors <- append(cors, cor(delta_pred, deltas0[selec0,j], method = 'pearson'))
#   }
# }
# cors <- matrix(cors,ncol=7, byrow = T)
# pdf(family="Helvetica",'correlation_variables_delta_p.pdf')
# heatmap.2(cors, trace="none",symm=TRUE, Rowv = NA, Colv = NA,
#           dendrogram = "none", keysize=1,margins=c(10,10), 
#           labCol=colnames(df)[6:13],
#           labRow=colnames(pred_2090_list), col= redgreen(15), 
#           symkey = F, density.info = 'none')
# dev.off()
# pdf(family="Helvetica",'correlation_variables_delta_p_0.pdf')
# heatmap.2(cors, trace="none",
#           dendrogram = "none", 
#           labCol=colnames(df)[6:13],margins=c(10,10),
#           labRow=colnames(pred_2090_list), col= redgreen(15))
# dev.off()

delta_t_n <- as.vector(data_f[[1]][1:360,115:150]-data_p[[1]][1:360,115:150])
delta_t_s <- as.vector(data_f[[1]][1:360,30:65]-data_p[[1]][1:360,30:65])

cos_lat <- matrix(data=NA, nrow=360, ncol=180)
local(
for (k in 1:360){
  for (l in 1:150){
    cos_lat[k,l]<<-cos(lati[l]*2*pi/360)
  }
}
)

w_n <- as.vector(cos_lat[1:360,115:150])
w_s <- as.vector(cos_lat[1:360,30:65])
mean_n <- sum(delta_t_n*w_n, na.rm=T)/sum(w_n[!is.na(delta_t_n)], na.rm=T)
mean_s <- sum(delta_t_s*w_s, na.rm=T)/sum(w_s[!is.na(delta_t_s)], na.rm=T)
write(c(mean_n, mean_s), 'mean_delta_T.txt')

pchs=c(15:18, 0:3)
for (i in 1:length(variables)){
  obs <- as.vector(data_w[[i]])
  mod <- as.vector(data[[i]])
  mod_c <- as.vector(data_p[[i]])
  obs0 <- obs[!is.na(obs) & !is.na(mod) &!is.na(mod_c)]
  mod0 <- mod[!is.na(obs) & !is.na(mod) &!is.na(mod_c)]
  mod_c0 <- mod_c[!is.na(obs) & !is.na(mod) &!is.na(mod_c)]
  pdf(family="Helvetica",paste('taylor_diagram',  colnames(df)[variables][i],'.pdf',sep=''))
  taylor.diagram(obs0,mod0, ref.sd=T,col='red', pch=pchs[i], pcex=2)
  taylor.diagram(obs0,mod_c0, pch=pchs[i], pcex=2, col='blue', add=T)
  dev.off()
}

color_var <- brewer.pal(length(variables), 'Dark2')
color_var1 <- brewer.pal(length(variables), 'Set2')
source('taylor.diagram_1.R')
pdf(family="Helvetica",'taylor_diagram_all.pdf')
for (i in c(6, 1:5, 7)){
  obs <- as.vector(data_w[[i]])
  mod <- as.vector(data[[i]])
  mod_c <- as.vector(data_p[[i]])
  obs0 <- obs[!is.na(obs) & !is.na(mod) &!is.na(mod_c)]
  mod0 <- mod[!is.na(obs) & !is.na(mod) &!is.na(mod_c)]
  mod_c0 <- mod_c[!is.na(obs) & !is.na(mod) &!is.na(mod_c)]
  obs1 <- (obs0-mean(obs0))/sd(obs0)
  mod1 <- (mod0-mean(obs0))/sd(obs0)
  mod_c1 <- (mod_c0-mean(obs0))/sd(obs0)
  if (i==6){
    taylor.diagram_1(obs1,mod1, ref.sd=T,col=color_var[i], pch=pchs[i], pcex=2)
    taylor.diagram_1(obs1,mod_c1, pch=pchs[i], pcex=2, col=color_var1[i], add=T)
  } else{
    taylor.diagram_1(obs1,mod1, ref.sd=F,add=T,col=color_var[i], pch=pchs[i], pcex=2)
    taylor.diagram_1(obs1,mod_c1, pch=pchs[i], pcex=2, col=color_var1[i], add=T)
  }
}
dev.off()
pdf(family="Helvetica",'the_taylor_legend.pdf')
plot(0,0, col='white', xaxt='n', yaxt='n', xlim=c(0, 5), ylim=c(0,5))
legend('topleft', legend=colnames(df)[variables], pch=pchs, ncol=2, col = color_var, bty='n')
dev.off()
pdf(family="Helvetica",'the_taylor_legend1.pdf')
plot(0,0, col='white', xaxt='n', yaxt='n', xlim=c(0, 5), ylim=c(0,5))
legend('topleft', legend=colnames(df)[variables], pch=pchs, ncol=2, col = color_var1, bty='n')
dev.off()

long<-seq(0.5,359.5,1)
long[181:360] =as.numeric(long[181:360])-360
lat<-seq(-89.5,89.5,1)
for (v in 1:length(variables)){
  maxi <- max(data_w[[v]],data_f[[v]],data_p[[v]],data[[v]],data1[[v]],   na.rm = T)
  mini<- min(data_w[[v]],data_f[[v]],data_p[[v]],data[[v]],data1[[v]],  na.rm = T)
  values <- round(seq(mini,maxi,(maxi-mini)/99), digits = 8)
  #colors3 <- rev(heat.colors(100))
  colors3 <- colorRampPalette(c('darkorchid4', 'blue', 'green', 'orange', 'red'))(100)
  data_plot <- NULL
  data_plot_f <- NULL
  data_plot_p <- NULL
  data_plot_f_n <- NULL
  data_plot_p_n <- NULL
  for (lg in 1:360){
    for (lt in 1:150){
      vec <- NULL
      vec_f <- NULL
      vec_p <- NULL
      vec_fn <- NULL
      vec_pn <- NULL
      if (!is.na(data_w[[v]][lg, lt]) & !(is.na(data_f[[v]][lg, lt])) & !(is.na(data_p[[v]][lg, lt]))){
        vec <- append(vec, c(long[lg], lat[lt], which.min(abs(values-data_w[[v]][lg, lt])),data_w[[v]][lg, lt]))
        vec_f <- append(vec_f, c(long[lg], lat[lt], which.min(abs(values-data_f[[v]][lg, lt])), data_f[[v]][lg, lt]))
        vec_p <- append(vec_p, c(long[lg], lat[lt], which.min(abs(values-data_p[[v]][lg, lt])), data_p[[v]][lg, lt]))
        vec_fn <- append(vec_fn, c(long[lg], lat[lt], which.min(abs(values-data1[[v]][lg, lt])), data1[[v]][lg, lt]))
        vec_pn <- append(vec_pn, c(long[lg], lat[lt], which.min(abs(values-data[[v]][lg, lt])), data[[v]][lg, lt]))
        data_plot <- rbind(data_plot, vec)
        data_plot_f <- rbind(data_plot_f, vec_f)
        data_plot_p <- rbind(data_plot_p, vec_p)
        data_plot_f_n <- rbind(data_plot_f_n, vec_fn)
        data_plot_p_n <- rbind(data_plot_p_n, vec_pn)
      }
    }
  }
  pnt=cbind(x =c(100,105,105,100), y =c(30,30,50,50))
  data_plot <- as.data.frame(data_plot)
  colnames(data_plot)<- c('lg', 'lt', 'pred', 'val')
  data_plot_f <- as.data.frame(data_plot_f)
  colnames(data_plot_f)<- c('lg', 'lt', 'pred', 'val')
  data_plot_p <- as.data.frame(data_plot_p)
  colnames(data_plot_p)<- c('lg', 'lt', 'pred', 'val')
  data_plot_f_n <- as.data.frame(data_plot_f_n)
  colnames(data_plot_f_n)<- c('lg', 'lt', 'pred', 'val')
  data_plot_p_n <- as.data.frame(data_plot_p_n)
  colnames(data_plot_p_n)<- c('lg', 'lt', 'pred', 'val')
  mini <- round(mini, digits = attr(regexpr("(?<=\\.)0+", mini, perl = TRUE), "match.length")+2)
  maxi <- round(maxi, digits = attr(regexpr("(?<=\\.)0+", maxi, perl = TRUE), "match.length")+2)

  pdf(family="Helvetica",file=paste(model0,'_',colnames(df[0, variables])[v],'_2006_woa.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  points(data_plot$lg, data_plot$lt, col =colors3[data_plot$pred],bg=colors3[data_plot$pred], pch=15, cex=0.4)
  SDMTools::legend.gradient(pnt,colors3,limits=c(mini,maxi),title=colnames(df[0, variables])[v],cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste(model0,'_',colnames(df[0, variables])[v],'_2090_cor.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  points(data_plot_f$lg, data_plot_f$lt, col =colors3[data_plot_f$pred],bg=colors3[data_plot_f$pred], pch=15, cex=0.4)
  SDMTools::legend.gradient(pnt,colors3,limits=c(mini,maxi),title=colnames(df[0, variables])[v],cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste(model0,'_',colnames(df[0, variables])[v],'_2006_cor.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  points(data_plot_p$lg, data_plot_p$lt, col =colors3[data_plot_p$pred],bg=colors3[data_plot_p$pred], pch=15, cex=0.4)
  SDMTools::legend.gradient(pnt,colors3,limits=c(mini,maxi),title=colnames(df[0, variables])[v],cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste(model0,'_',colnames(df[0, variables])[v],'_2090.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  points(data_plot_f_n$lg, data_plot_f_n$lt, col =colors3[data_plot_f_n$pred],bg=colors3[data_plot_f_n$pred], pch=15, cex=0.4)
  SDMTools::legend.gradient(pnt,colors3,limits=c(mini,maxi),title=colnames(df[0, variables])[v],cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste(model0,'_',colnames(df[0, variables])[v],'_2006.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  points(data_plot_p_n$lg, data_plot_p_n$lt, col =colors3[data_plot_p_n$pred],bg=colors3[data_plot_p_n$pred], pch=15, cex=0.4)
  SDMTools::legend.gradient(pnt,colors3,limits=c(mini,maxi),title=colnames(df[0, variables])[v],cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  
  if (v == 6){
    maxi <- max((data_plot_f$val-data_plot_p$val)[(data_plot_f$val-data_plot_p$val)<0.01])
    mini <- min(data_plot_f$val-data_plot_p$val)
  } else if (v==3){
    maxi <- max((data_plot_f$val-data_plot_p$val)[(data_plot_f$val-data_plot_p$val)<3])
    mini <- min((data_plot_f$val-data_plot_p$val)[(data_plot_f$val-data_plot_p$val)>-3])
  } else{
    maxi <- max(data_plot_f$val-data_plot_p$val)
    mini <- min(data_plot_f$val-data_plot_p$val)
  }
  max_ok <- max(c(abs(maxi), abs(mini)))
  values_d <- round(seq(-max_ok,max_ok,2*max_ok/99), digits = 8)
  col_index <- NULL
  for (j in 1:dim(data_plot_f)[1]){
    col_index <- append(col_index, which.min(abs(values_d-(data_plot_f$val-data_plot_p$val)[j])))
  }
  colfunc <- colorRampPalette(c('blue', 'lightblue','white', 'orange', 'red'))(100)
  mini <- round(mini, digits = attr(regexpr("(?<=\\.)0+", mini, perl = TRUE), "match.length")+2)
  maxi <- round(maxi, digits = attr(regexpr("(?<=\\.)0+", maxi, perl = TRUE), "match.length")+2)
  max_ok <- round(max_ok, digits = attr(regexpr("(?<=\\.)0+", max_ok, perl = TRUE), "match.length")+1)
  data_contour <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(data_plot_f)[1]){
    data_contour[which(longi_sorted==data_plot_f$lg[i]),which(lati==data_plot_f$lt[i])]= col_index[i]
  }
  pdf(family="Helvetica",file=paste(model0,'_',colnames(df[0, variables])[v],'_2090_2006.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(x=longi_sorted, y=lati_sorted, z=data_contour, levels = c(1:100), col=colfunc)
  
  # points(data_plot_f$lg, data_plot_f$lt, col =colfunc[col_index],bg=colfunc[col_index], pch=15, cex=1)
  SDMTools::legend.gradient(pnt,colfunc,limits=c(-maxi,maxi),title=colnames(df[0, variables])[v],cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  cond_oor <- data_plot_f$val < min(data_plot_p$val, na.rm=T) | data_plot_f$val > max(data_plot_p$val, na.rm=T)
  co <- 1
  for (i in unique(data_plot_p$lg)){
    q <- sum(data_plot_p$lg==i)
    if (co%%4==1 ){
      cond_oor[data_plot_p$lg==i][c(1:q)[!(seq(1,q,1) %in% seq(1, q,4))]] <- F
    } else if(co%%4==3){
      cond_oor[data_plot_p$lg==i][c(1:q)[!(seq(1,q,1) %in% seq(3, q,4))]] <- F
    } else {
      cond_oor[data_plot_p$lg==i] <- F
    } 
    co <- co+1
  }
  
  pdf(family="Helvetica",file=paste(model0,'_',colnames(df[0, variables])[v],'_2090_2006_oor.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(x=longi_sorted, y=lati_sorted, z=data_contour, levels = c(1:100), col=colfunc)
  points(x=data_plot_p$lg[cond_oor], y =data_plot_p$lt[cond_oor], col='black',pch=8,cex=0.25)
  # points(data_plot_f$lg, data_plot_f$lt, col =colfunc[col_index],bg=colfunc[col_index], pch=15, cex=1)
  #SDMTools::legend.gradient(pnt,colfunc,limits=c(-maxi,maxi),title=colnames(df[0, variables])[v],cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
}

set = colorRampPalette(c('blue','yellow', 'red'))(101)
set1 = colorRampPalette(c('magenta', 'magenta4'))(21)
varis <- c('thetao','so', 'no3','po4','si', 'dfe','chl', 'si-no3')
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
lati<-seq(-89.5, 89.5, 1)
data_p_n <- readRDS('data_p_n.rds')
data_f_n <- readRDS('data_f_n.rds')
data_list<-list(data_f_n, data_p_n)
count=1
for (data in data_list){
  for (u in 1:length(variables)){
    if (u != 1 & min(data[[u]], na.rm = T)<0){
      if (count==1){
        pdf(family="Helvetica",file=paste(model0,'_', varis[u], '_2090_neg.pdf', sep=''),width=10,height=4.065)
      } else{
        pdf(family="Helvetica",file=paste(model0,'_', varis[u], '_2006_neg.pdf', sep=''),width=10,height=4.065)
      }
      par(mar=c(0,0,0,0))
      maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
      plot(coastline,lwd=0.0475, col='black', add=T)
      pnt=cbind(x =c(100,105,105,100), y =c(30,30,50,50))
      pnt1=cbind(x =c(80,85,85,80), y =c(30,30,50,50))
      mx = max(data[[u]], na.rm = T)
      if (u== 1){
        mi = min(data[[u]], na.rm = T)
        mi_n = 1
      } else{
        mi = min(data[[u]][data[[u]]>0], na.rm = T)
        mi_n = min(data[[u]], na.rm = T)
        if (mi_n<0){
          vals1 <- seq(mi_n,0, abs(mi_n/20) )
        } else{
          vals1 <- NULL
        }
      }
      vals <- seq(mi, mx, (mx-mi)/100)
      for (i in 1:360){
        for (j in 1:150){
          if (!(is.na(data[[u]][i,j]))){
            if (u!=1 & data[[u]][i,j]<0){
              id <- which.min(abs(vals1 - data[[u]][i,j]))
              points(longi[i], lati[j], col=set1[id],pch=15,cex=0.4)
            } else{
              id <- which.min(abs(vals - data[[u]][i,j]))
              points(longi[i], lati[j], col=set[id],pch=15,cex=0.4)
            }
          }
        }
      }
      SDMTools::legend.gradient(pnt,set,limits=c(round(mi, digits = 3),round(mx, digits = 3)),title=varis[u],cex=0.4)
      if (u != 1 & mi_n <0){
        SDMTools::legend.gradient(pnt1,set1,limits=c(round(mi_n, digits = 3),round(0, digits = 3)),title=paste('negative',varis[u]),cex=0.4)
      }
      hide_arctic()
      axis_map0()
      dev.off()
    }
  }
  count=count+1
}


