#!bin/usr/bin/env Rscript
library('Rcpp')
library("gbm")
library('randomForest')
library('mgcv')
library('nnet')
library("dismo")
library("FactoMineR")
library("gplots")
library("stringr")
library('mapproj')
library('mapplots')
library('maptools')
library('SDMTools')
library('RColorBrewer')
library('extrafont')
library('ncdf4')
library("CDFt")
library('isofor')
library('parallel')
library('png')
library('grid')
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
source('vioplot.R')
source('hide_arctic.R')
source('axis_map0.R')
cambria <-readRDS('cambria.rds')
test <- ncdf4::nc_open('Time_Varying_Biomes.nc')
biome <- ncvar_get(test, 'MeanBiomes')
biome[is.nan(biome)]<-0
biome <-t(biome)
lt <- ncvar_get(test, 'lat')
lo <- ncvar_get(test, 'lon')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
BGCP_shp <- readShapeSpatial('longhurst_v4_2010/Longhurst_world_v4_2010.shp')
reyg_bgcp <- as.matrix(read.csv2('Reygondeau/PROVINCE_REYGONDEAU.csv', header = F, sep = ','))
reyg_bio <- as.matrix(read.csv2('Reygondeau/BIOME_REYGONDEAU.csv', header = F, sep = ','))
reyg_bgcp19_1 <-as.matrix(readRDS('BGCP_2019_REYGONDEAU_mat_1.rds'))
#reyg_bgcp[is.nan(reyg_bgcp)]<-0
reyg_bio[is.nan(reyg_bio)]<-0
reyg_bgcp19_1[is.nan(reyg_bgcp19_1)]<-0
to_plot <-list(reyg_bgcp19_1, reyg_bgcp, reyg_bio)
names <- c('Reygondeau_BGCP19','Reygondeau_BGCP','Reygondeau_BIOME' )

lati<-seq(-89.5, 60, 1)
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
lati_sorted<-seq(-89.5, 60, 1)
longi_sorted =sort(longi)

model0<-'model-mean'

# communities0<-readRDS(paste(model0,'_communities_woa.rds', sep=''))
# communities<-readRDS(paste(model0,'_communities_2006.rds', sep=''))
# communities1<-readRDS(paste(model0,'_communities_2090.rds', sep=''))
# mix0 <- readRDS(paste(model0,'_communities_mix_woa.rds', sep=''))
# mix <- readRDS(paste(model0,'_communities_mix_2006.rds', sep=''))
# mix1 <- readRDS(paste(model0,'_communities_mix_2090.rds', sep=''))
# delta_mix<-readRDS(paste(model0,'_delta_mix.rds', sep=''))
bray_curtis<- readRDS(paste(model0,'_bray_curtis.rds', sep=''))
pred_woa_list <-readRDS(paste(model0, '_pred_woa_list.rds', sep=''))
pred_2006_list<-readRDS(paste(model0, '_pred_2006_list.rds', sep=''))
pred_2090_list<-readRDS(paste(model0, '_pred_2090_list.rds', sep=''))
pred_2090_list_noT<-readRDS(paste(model0, '_pred_2090_list_noT.rds', sep=''))
# drivers<-readRDS(paste(model0, '_drivers.rds', sep=''))
drivers_al<-readRDS(paste(model0, '_drivers_al.rds', sep=''))
weights_cos<-readRDS(paste(model0, '_weights_cos.rds', sep=''))
DS <- read.table('Distances_shifts_model-mean.txt', header = F)
colnames(DS)<- c("Fraction" ,"Metacommunity" , "Model","Location","Lat_2006", 'Long_2006',"Lat_2090", 'Long_2090',"Shift" ,"Lat_shift" ,"Long_shift")
frcs <- c('0-0.2','0.22-3','0.8-5', '180-2000', '20-180', '43952')
DS$Fraction <- frcs[DS$Fraction]

data_fisheries<- readRDS('data_fisheries.rds')


fractions <- c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')

best_models <- read.table('best_selected_models.txt', header = T)
best_models <- as.data.frame(best_models)
N <- dim(best_models)[1]
row.names(best_models)<-c(1:N)
fractions0 <- unique(best_models$Fraction)

df <- readRDS("Genocenoses_env_parameters_woa.rds")
variables <- c(6:11,13)

color1 = colorRampPalette(c("cyan", "dark blue"))(13)
color2 = colorRampPalette(c("darkgreen","light blue", "dark violet"))(15)
colors2 <- colorRampPalette(c('lightblue','lightgreen', 'orange'))(100)

cond_na <- bray_curtis$Long<30 & bray_curtis$Long > -80 & bray_curtis$Lat >0 
cond_sa <- bray_curtis$Long<20 & bray_curtis$Long > -70 & bray_curtis$Lat < 0
cond_np <- (bray_curtis$Long>120 | bray_curtis$Long < -90) & bray_curtis$Lat >0
cond_sp <- (bray_curtis$Long>120 | bray_curtis$Long < -70) & bray_curtis$Lat <0
cond_i <- bray_curtis$Long>20 & bray_curtis$Long < 120 & bray_curtis$Lat <30
basins_conds <- list(cond_na, cond_sa, cond_np, cond_sp, cond_i)
basins <- c('North_A', 'South_A', 'North_P', 'South_P', 'Indian')

tot_area_bc <- sum(111*111*weights_cos)
covered_areas_median_change <- 100*sum(111*111*weights_cos*as.numeric((100-bray_curtis$b_c)>median(100-bray_curtis$b_c)))/tot_area_bc

bray_curtis_vs_fish <- list(NULL,NULL,NULL,NULL, NULL)
color_code<- c('blue', 'green', 'yellow', 'orange', 'red')
bray_curtis$fish_decile <- data_fisheries$deciles[match(bray_curtis$cell, data_fisheries$Cell)]
areas_fish <- list(NULL,NULL,NULL,NULL, NULL)
tot_area <- 0
bc_all <- NULL
percent_affected_area <- NULL
cond_fish <- bray_curtis$fish_decile==9 | bray_curtis$fish_decile==10 | bray_curtis$fish_decile==8 | bray_curtis$fish_decile==7
for (u in 1:5){
  bc_vals<-bray_curtis$b_c[bray_curtis$fish_decile==2*u-1 | bray_curtis$fish_decile==2*u]
  areas <- as.numeric(bray_curtis$fish_decile==2*u-1 | bray_curtis$fish_decile==2*u)*111*111*weights_cos
  areas <-areas[bray_curtis$fish_decile==2*u-1 | bray_curtis$fish_decile==2*u]
  bc_vals <- bc_vals[!is.na(bc_vals)]
  ok_areas <- areas[!is.na(areas)]
  areas_fish[[u]] <- ok_areas
  bray_curtis_vs_fish[[u]]<-bc_vals
  bc_all <-append(bc_all, bc_vals)
  percent_affected_area<- append(percent_affected_area,sum(ok_areas[(100-bc_vals)>(100-median(bray_curtis$b_c))])/sum(ok_areas))
  tot_area <- tot_area+sum(ok_areas[100-bc_vals>median(100-bray_curtis$b_c)])
}
bray_curtis$cond_fish <- cond_fish

high_fish <- 100*(sum(areas_fish[[5]][(100-bray_curtis_vs_fish[[5]])>median(100-bray_curtis$b_c)])+
                  sum(areas_fish[[4]][(100-bray_curtis_vs_fish[[4]])>median(100-bray_curtis$b_c)]))/(sum(areas_fish[[5]])+sum(areas_fish[[4]]))

pdf(family="Helvetica",file=paste(model0,'_fish_vs_bc.pdf', sep=''),width=20,height=10)
plot(c(10,10), main='',  xlim=c(0,90), 
     ylim=c(0,max(density(100-bray_curtis_vs_fish[[1]], cut=7)$y, density(100-bray_curtis_vs_fish[[2]], cut=7)$y, 
                  density(100-bray_curtis_vs_fish[[3]], cut=7)$y, density(100-bray_curtis_vs_fish[[4]], cut=7)$y,
                  density(100-bray_curtis_vs_fish[[5]], cut=7)$y)), xlab='Bray-Curtis similarity index',ylab='Density', col='white',xaxs="i")
for (count in 1:5){
  d <- density(100-bray_curtis_vs_fish[[count]], cut=7)
  d$y[d$x < 0] <- 0
  d$y[d$x >100] <- 0
  polygon(d, col=scales::alpha(color_code[count], 0.5), border=scales::alpha(color_code[count], 0.5))
  legend('topright', legend = c('Low', '', '','', 'High'), col=color_code,pch=15,cex=2, bty='n', title='Fisheries')
}
dev.off()

pdf(family="Helvetica",file=paste(model0,'_fish_vs_bc_violin.pdf', sep=''),width=20,height=15)
vioplot(100-bray_curtis_vs_fish[[1]], 100-bray_curtis_vs_fish[[2]],100-bray_curtis_vs_fish[[3]],
        100-bray_curtis_vs_fish[[4]],100-bray_curtis_vs_fish[[5]], col=scales::alpha(color_code, 0.6), 
        names = c('Low', '', 'Medium', '', 'High'), cex = 2.5)
dev.off()

eez_zones <- readRDS('model-mean_eez_points1.rds')
eez_cells <-readRDS('model-mean_eez_points.rds')
cond_eez <- bray_curtis$cell %in% eez_cells
bray_curtis$cond_eez <- cond_eez
saveRDS(bray_curtis, 'model-mean_bray_curtis_cond_fish_eez.rds')
sum((100-bray_curtis$b_c[cond_eez])*weights_cos[cond_eez])/sum(weights_cos[cond_eez])

cond_trop <- bray_curtis$Lat < 25 & bray_curtis$Lat > -25
cond_n_temp <- bray_curtis$Lat >25
cond_s_temp <- bray_curtis$Lat > -60 & bray_curtis$Lat < -25
cond_polar <- bray_curtis$Lat < -60
cond_climate <- list(cond_trop, cond_n_temp, cond_s_temp, cond_polar)
climates <- c('Tropical', 'North_temperate', 'South_temperate', 'Polar')
sum((100-bray_curtis$b_c[cond_trop])*weights_cos[cond_trop])/sum(weights_cos[cond_trop])
sum((100-bray_curtis$b_c[cond_n_temp])*weights_cos[cond_n_temp])/sum(weights_cos[cond_n_temp])
sum((100-bray_curtis$b_c[cond_s_temp])*weights_cos[cond_s_temp])/sum(weights_cos[cond_s_temp])
sum((100-bray_curtis$b_c[cond_polar])*weights_cos[cond_polar])/sum(weights_cos[cond_polar])

glob_stats_fish_eez_basins <- NULL
c<-1
for (cond in basins_conds){
  mean_bc <- sum((100-bray_curtis$b_c[cond])*weights_cos[cond])/sum(weights_cos[cond])
  mean_bc_eez <- sum((100-bray_curtis$b_c[cond & cond_eez])*weights_cos[cond & cond_eez])/sum(weights_cos[cond & cond_eez])
  mean_bc_fish <- sum((100-bray_curtis$b_c[cond & cond_fish])*weights_cos[cond & cond_fish], na.rm = T)/sum(weights_cos[cond & cond_fish], na.rm=T)
  glob_stats_fish_eez_basins <- rbind(glob_stats_fish_eez_basins, c(basins[c], mean_bc,mean_bc_eez, mean_bc_fish))
  c<-c+1
}
write.table(glob_stats_fish_eez_basins, 'glob_stats_fish_eez_basins.txt')

glob_stats_fish_eez_climate <- NULL
c<-1
for (cond in cond_climate){
  mean_bc <- sum((100-bray_curtis$b_c[cond])*weights_cos[cond])/sum(weights_cos[cond])
  mean_bc_eez <- sum((100-bray_curtis$b_c[cond & cond_eez])*weights_cos[cond & cond_eez])/sum(weights_cos[cond & cond_eez])
  mean_bc_fish <- sum((100-bray_curtis$b_c[cond & cond_fish])*weights_cos[cond & cond_fish], na.rm = T)/sum(weights_cos[cond & cond_fish], na.rm=T)
  glob_stats_fish_eez_climate <- rbind(glob_stats_fish_eez_climate, c(climates[c], mean_bc,mean_bc_eez, mean_bc_fish))
  c<-c+1
}
write.table(glob_stats_fish_eez_climate, 'glob_stats_fish_eez_climate.txt')

# covered_areas <- NULL
# for (k in unique(communities$Metacommunity)){
#   cond0 <- communities0$Metacommunity == k
#   cond <- communities$Metacommunity == k
#   cond1 <- communities1$Metacommunity == k
#   area0 <- sum(as.numeric(cond0)*111*111*weights_cos)
#   area <- sum(as.numeric(cond)*111*111*weights_cos)#area in km^2
#   area1 <- sum(as.numeric(cond1)*111*111*weights_cos)#area in km^2
#   delta <- area1 - area
#   covered_areas <- rbind(covered_areas, c(k, area0,area, area1, delta))
# }
# colnames(covered_areas)<- c('Metacommunity', 'Area woa','Area 2006', 'Area 2090', 'Delta area (2090-2006)')
# write.table(covered_areas, paste(model0,'_covered_areas_2006vs2090.txt', sep=''), sep="\t", row.names = F)
# 
# covered_areas <- NULL
# for (k in unique(communities1$Metacommunity)){
#   cond0 <- communities0$Metacommunity == k
#   cond <- communities$Metacommunity == k
#   cond1 <- communities1$Metacommunity == k
#   area0 <- sum(as.numeric(cond0)*111*111*weights_cos)
#   area <- sum(as.numeric(cond)*111*111*weights_cos)#area in km^2
#   area1 <- sum(as.numeric(cond1)*111*111*weights_cos)#area in km^2
#   delta <- area1 - area
#   covered_areas <- rbind(covered_areas, c(k, area0,area, area1, delta))
# }
# colnames(covered_areas)<- c('Metacommunity','Area woa' ,'Area 2006', 'Area 2090', 'Delta area (2090-2006)')
# write.table(covered_areas, paste(model0,'_covered_areas_2090vs2006.txt', sep=''), sep="\t", row.names = F)


# coverage_mix <- NULL
# for (k in 1:max(max(mix$Number, mix1$Number))){
#   cond0 <- mix0$Number == k
#   cond <- mix$Number == k
#   cond1 <- mix1$Number == k
#   area0 <- sum(as.numeric(cond0)*111*111*mix$cosLat)
#   area <- sum(as.numeric(cond)*111*111*mix$cosLat) #area in km^2
#   area1 <- sum(as.numeric(cond1)*111*111*mix1$cosLat)#area in km^2
#   delta <- area1-area
#   coverage_mix <- rbind(coverage_mix, c(k-1, area0,area, area1, delta))
# }
# colnames(coverage_mix)<- c('Number','Area woa' ,'Area 2006', 'Area 2090', 'Delta area (2090-2006)')
# write.table(coverage_mix, paste(model0,'_coverage_mix.txt', sep=''), sep="\t", row.names = F)
# 
# pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_mix_woa.pdf', sep=''),width=10,height=4.065)
# maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
# points(x=mix0$Long, y =mix0$Lat, col=color1[mix0$Number],bg= color1[mix0$Number],pch=15,cex=0.229)
# legend(x = 90, y =60, legend = c(0:max(mix0$Number)), col=color1[1:(max(mix0$Number)+1)],pch=15,cex=8)
# text(x= 100, y = 65, labels = 'Number of \nmetacommunities' ,cex=6)
# dev.off()
# 
# pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_mix_2006-15.pdf', sep=''),width=10,height=4.065)
# maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
# points(x=mix$Long, y =mix$Lat, col=color1[mix$Number],bg= color1[mix$Number],pch=15,cex=0.229)
# legend(x = 90, y =60, legend = c(0:max(mix$Number)), col=color1[1:(max(mix$Number)+1)],pch=15,cex=8)
# text(x= 100, y = 65, labels = 'Number of \nmetacommunities' ,cex=6)
# dev.off()
# 
# pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_mix_2090-99.pdf', sep=''),width=10,height=4.065)
# maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
# points(x=mix1$Long, y =mix1$Lat, col=color1[mix1$Number],bg= color1[mix1$Number],pch=15,cex=0.229)
# legend(x = 90, y =60, legend = c(0:max(mix1$Number)), col=color1[1:(max(mix1$Number)+1)],pch=15,cex=8)
# text(x= 100, y = 65, labels = 'Number of \nmetacommunities' ,cex=6)
# dev.off()

# colnames(delta_mix) = c('Long', 'Lat', 'delta')
# delta_mix <- as.data.frame(delta_mix)
# pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_delta_mix_2090-99.pdf', sep=''),width=10,height=4.065)
# maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
# points(x=delta_mix$Long, y =delta_mix$Lat, col=color2[8+delta_mix$delta],bg= color2[8+delta_mix$delta],pch=15,cex=0.229)
# legend(x = 90, y =60, legend = c(-5:5), col=color2[(min(8+delta_mix$delta)):(max(8+delta_mix$delta))],pch=15,cex=8)
# text(x= 100, y = 65, labels = 'Delta number \nof metacommunities' ,cex=6)
# dev.off()

# communities$Metacommunity <- as.character(communities$Metacommunity)
# communities1$Metacommunity <- as.character(communities1$Metacommunity)
# lats <- NULL
# longs <- NULL
# cols <- NULL
# for (i in 1:dim(communities)[1]){
#   if (communities$Metacommunity[i] %in% unique(communities1$Metacommunity)){
#     if (communities1$Metacommunity[i] == communities$Metacommunity[i]){
#       lats <- append(lats, as.numeric(mix$Lat[i]))
#       longs <- append(longs, as.numeric(mix$Long[i]))
#       cols <- append(cols, 'blue')
#     } else if ( !(communities1$Metacommunity[i] %in% unique(communities$Metacommunity)) ){
#       lats <- append(lats, as.numeric(mix$Lat[i]))
#       longs <- append(longs, as.numeric(mix$Long[i]))
#       cols <- append(cols,'darkorchid1')
#     } else {
#       lats <- append(lats, as.numeric(mix$Lat[i]))
#       longs <- append(longs, as.numeric(mix$Long[i]))
#       cols <- append(cols,'light blue')
#     }
#   } else {
#     lats <- append(lats, as.numeric(mix$Lat[i]))
#     longs <- append(longs, as.numeric(mix$Long[i]))
#     cols <- append(cols, 'red')
#   }
# }
# pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_changes.pdf', sep=''),width=10,height=4.065)
# maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
# points(x=longs, y =lats, col=cols,bg=cols,pch=15,cex=0.229)
# legend(x = 60, y =60, legend = c('Community found only in 2006','Community found only in 2090','Change of community','Same community'), col=c('red', 'darkorchid1','light blue', 'blue'),pch=15,cex=8)
# dev.off()

values <- 100-100*(bray_curtis$b_c-min(bray_curtis$b_c))/((max(bray_curtis$b_c)-min(bray_curtis$b_c)))+1
data_contour <- matrix(NA, ncol=150, nrow=360)
for (i in 1:dim(bray_curtis)[1]){
  data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= round(values[i])+1
}
set = colorRampPalette(c('blue','yellow', 'red'))(101)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_2006-2090.pdf', sep='')
    ,width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
# points(x=bray_curtis$Long, y =bray_curtis$Lat, col=set[values],bg= set[values],pch=15,cex=0.229)
# SDMTools::legend.gradient(pnt,set,limits=c(0,100-round(min(bray_curtis$b_c), digits = 1)),title='Bray-curtis \ndissimilarity index',cex=0.4)
.filled.contour(x=longi_sorted,y=lati, z=data_contour,  col=set, levels = c(1:100))
hide_arctic()
axis_map0()
dev.off()

values <- 100-100*(bray_curtis$b_c-min(bray_curtis$b_c))/((max(bray_curtis$b_c)-min(bray_curtis$b_c)))+1
set = colorRampPalette(c('blue','yellow', 'red'))(101)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_high_fish_2006-2090.pdf',
               sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
selection <- bray_curtis$fish_decile==9 | bray_curtis$fish_decile==10
selection[is.na(selection)]<-FALSE
points(x=bray_curtis$Long[selection], y =bray_curtis$Lat[selection],
       col=set[values[selection]],bg= set[values[selection]],
       pch=15,cex=0.229)
SDMTools::legend.gradient(pnt,set,limits=c(0,100-round(min(bray_curtis$b_c), digits = 1)),title='Bray-curtis \ndissimilarity index',cex=0.4)
hide_arctic()
axis_map0()
dev.off()

values <- 100-100*(bray_curtis$b_c-min(bray_curtis$b_c))/((max(bray_curtis$b_c)-min(bray_curtis$b_c)))+1
set = colorRampPalette(c('blue','yellow', 'red'))(101)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_eez_2006-2090.pdf',
               sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
selection <- cond_eez
points(x=bray_curtis$Long[selection], y =bray_curtis$Lat[selection],
       col=set[values[selection]],bg= set[values[selection]],
       pch=15,cex=0.229)
SDMTools::legend.gradient(pnt,set,limits=c(0,100-round(min(bray_curtis$b_c), digits = 1)),title='Bray-curtis \ndissimilarity index',cex=0.4)
hide_arctic()
axis_map0()
dev.off()

# maxi <- max(drivers[[1]][,3],drivers[[2]][,3],drivers[[3]][,3],drivers[[4]][,3],drivers[[5]][,3],drivers[[6]][,3],drivers[[7]][,3],na.rm=T )
# mini <- min(drivers[[1]][,3],drivers[[2]][,3],drivers[[3]][,3],drivers[[4]][,3],drivers[[5]][,3],drivers[[6]][,3],drivers[[7]][,3],na.rm=T )
# colors2 <- colorRampPalette(c('lightblue','lightgreen', 'orange'))(100)
# values <- seq(mini, maxi, (maxi-mini)/99)
# pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
# for (e in 1:length(variables)){
#   data_d <- drivers[[e]]
#   data_d <- as.data.frame(data_d)
#   colnames(data_d)<- c('long', 'lat', 'pred', 'alph')
#   data_d$pred <- round(data_d$pred*100)+1
#   alpha_col <- data_d$alph+(1-max(data_d$alph, na.rm=T))
#   pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_',colnames(df)[variables][e],'.pdf', sep=''),width=10,height=4.065)
#   maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
#   points(data_d$long, data_d$lat, col=scales::alpha(colors2[data_d$pred], alpha_col), pch=15,cex=0.229)
#   SDMTools::legend.gradient(pnt,colors2,limits=c(round(mini*100),round(maxi*100)),title=paste(colnames(df)[variables][e],'\nDriver impact (%)', sep=''),cex=0.4)
#   dev.off()
# }

maxi <- max(drivers_al[[1]][,3],drivers_al[[2]][,3],drivers_al[[3]][,3],drivers_al[[4]][,3],drivers_al[[5]][,3],drivers_al[[6]][,3],drivers_al[[7]][,3],na.rm=T )
mini <- min(drivers_al[[1]][,3],drivers_al[[2]][,3],drivers_al[[3]][,3],drivers_al[[4]][,3],drivers_al[[5]][,3],drivers_al[[6]][,3],drivers_al[[7]][,3],na.rm=T )
values <- seq(mini, maxi, (maxi-mini)/99)
data_d_al <- NULL
for (e in 1:length(variables)){
  data_d <- drivers_al[[e]]
  data_d_al <-cbind(data_d_al,data_d[,3])
  data_d <- as.data.frame(data_d)
  colnames(data_d)<- c('lat', 'long', 'pred', 'alph')
  data_d$pred <- round(data_d$pred*100)+1
  alpha_col <- data_d$alph+(1-max(data_d$alph, na.rm=T))
  pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_',colnames(df)[variables][e],'_alone.pdf', 
                 sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  points(data_d$long, data_d$lat, col=scales::alpha(colors2[data_d$pred], alpha_col), pch=15,cex=0.229)
  SDMTools::legend.gradient(pnt,colors2,limits=c(round(mini*100),round(maxi*100)),
                            title=paste(colnames(df)[variables][e],'\nDriver impact (%)', sep=''),cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
}
data_d_al <- cbind(data_d_al, data_d[,1:2])
data_d_al <- as.data.frame(data_d_al)
colnames(data_d_al) <- c(colnames(df)[variables], 'lat','long')
data_d_al$maxi <- apply(data_d_al[,1:7], 1, which.max)
color_drivers <- brewer.pal(length(variables), "Dark2")
alphas_drivers <- apply(data_d_al[,1:7], 1, max)
second_max <- function(vec){second_max<-which(order(vec, decreasing=T)==2)}
second_max_value <- function(vec){second_max_value<-vec[order(vec, decreasing=T)==2]}
data_d_al$second_maxi <- apply(data_d_al[,1:7], 1, second_max)
alphas_second <- apply(data_d_al[,1:7], 1, second_max_value)
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
points(data_d_al$long, data_d_al$lat, col=scales::alpha(color_drivers[data_d_al$maxi], alphas_drivers), pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(data_d_al)[1:7], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

data_contour<- matrix(NA, ncol=150, nrow=360)
data_contour_bc <- matrix(NA, ncol=150, nrow=360)
for (i in 1:dim(data_d_al)[1]){
  data_contour[which(longi_sorted==data_d_al$long[i]),which(lati==data_d_al$lat[i])]= data_d_al$maxi[i]
  if ((1-bray_curtis$b_c[i])>1/6){
    data_contour_bc[which(longi_sorted==data_d_al$long[i]),which(lati==data_d_al$lat[i])]= data_d_al$maxi[i]
  }
}
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_contour.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
.filled.contour(x=longi_sorted, y=lati, z=data_contour, col=color_drivers, levels = c(1:8))
contour(x=longi_sorted,y=lati, z=data_contour,col='white',add=T, lwd=0.4,drawlabels=F)
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_contour_coff.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
.filled.contour(x=longi_sorted, y=lati, z=data_contour_bc, col=color_drivers, levels = c(1:8))
contour(x=longi_sorted,y=lati, z=data_contour_bc,col='white',add=T, lwd=0.4,drawlabels=F)
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_1.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
selec <- alphas_drivers>1.5/7
points(data_d_al$long[selec], data_d_al$lat[selec], col=scales::alpha(color_drivers[data_d_al$maxi[selec]], alphas_drivers[selec]), pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(data_d_al)[1:7], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_max_alone_2.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
selec <- alphas_drivers>1.5/7
points(data_d_al$long[selec], data_d_al$lat[selec], col=color_drivers[data_d_al$maxi[selec]], pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(data_d_al)[1:7], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_second_max_alone.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
points(data_d_al$long, data_d_al$lat, col=scales::alpha(color_drivers[data_d_al$second_maxi], 1), pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(data_d_al)[1:7], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

remaining_expl <- 1.5*(1-alphas_drivers)/7
pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_second_max_alone_1.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
points(data_d_al$long[alphas_second>remaining_expl], data_d_al$lat[alphas_second>remaining_expl], col=scales::alpha(color_drivers[data_d_al$second_maxi[alphas_second>remaining_expl]], 1), pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(data_d_al)[1:7], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

pdf(family="Helvetica",file=paste('projection_',model0,'_drivers_second_max_alone_2.pdf', sep=''),
    width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
points(data_d_al$long[alphas_second>1/7], data_d_al$lat[alphas_second>1/7], col=scales::alpha(color_drivers[data_d_al$second_maxi[alphas_second>1/7]], 1), pch=15,cex=0.229)
legend(x = 60, y = 60, legend = colnames(data_d_al)[1:7], col=color_drivers,pch=15,cex=8, ncol=2)
hide_arctic()
axis_map0()
dev.off()

means_d <-NULL
vals_drivers <- NULL
for (i in 1:7){
  mean_d <- sum(drivers_al[[i]][,3]*weights_cos*bray_curtis$b_c)/sum(weights_cos*bray_curtis$b_c)
  means_d <- rbind(means_d, c(colnames(df)[variables[i]], mean_d))
  vals_drivers<-cbind(vals_drivers,drivers_al[[i]][,3])
}
write.table(means_d, 'drivers_impact.txt')
pdf(family="Helvetica",'correlation_drivers_max.pdf')
heatmap.2(cor(vals_drivers), trace="none",symm=TRUE, Rowv = NA, Colv = NA,
          dendrogram = "none", keysize=1,margins=c(10,22), 
          labCol=colnames(df)[variables],
          labRow=colnames(df)[variables], col= greenred(15), symkey = F, density.info = 'none')
dev.off()

community_change <- bray_curtis[,1:2]
fractions <- c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
labs <- NULL
count=1
for (u in best_models$Fraction){
  gen <- best_models$Gen[count]
  if (u != "43952"){
    labs <-append(labs, paste(u, gen, sep='_'))
  } else{
    labs <- append(labs, paste('5-20', gen, sep='_'))
  }
  count=count+1
}



changes_metac<-NULL
changes_metac_basins<-NULL
dominant_communities_woa <- NULL
dominant_communities_2006 <- NULL
dominant_communities_2090 <- NULL
dominant_communities_2090_noT <- NULL
shifts <- NULL
areas <- NULL
#color_sets <-list()
color_sets <- readRDS('colors_provinces.rds')
number_of_changes <- rep(0, dim(bray_curtis)[1])
fractions0 <- unique(best_models$Fraction)
colors1 <- c('saddlebrown', 'red', 'dodgerblue2', 'darkgreen', 'darkviolet', 'darkorange')
for (frac in fractions0){
  goods <- best_models$Gen[best_models$Fraction==frac]
  #color_sets[[which(fractions0==frac)]]<-colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(max(goods))
  frac_ok <- fractions[which(fractions0==frac)]
  preds0 <- pred_woa_list[,grep(frac_ok, labs)]
  values0 <- apply(preds0,1, which.max)
  alphas0 <- apply(preds0,1, max)
  dom_com0 <- goods[values0]
  dominant_communities_woa<-cbind(dominant_communities_woa, dom_com0)
  dominant_communities_woa<-cbind(dominant_communities_woa, alphas0)
  preds <- pred_2006_list[,grep(frac_ok, labs)]
  values <- apply(preds,1, which.max)
  alphas <- apply(preds,1, max)
  dom_com <- goods[values]
  dominant_communities_2006<-cbind(dominant_communities_2006, dom_com)
  dominant_communities_2006<-cbind(dominant_communities_2006, alphas)
  preds1 <- pred_2090_list[,grep(frac_ok, labs)]
  values1 <- apply(preds1,1, which.max)
  alphas1 <- apply(preds1,1, max)
  dom_com1 <- goods[values1]
  dominant_communities_2090<-cbind(dominant_communities_2090, dom_com1)
  dominant_communities_2090<-cbind(dominant_communities_2090, alphas1 )
  preds2 <- pred_2090_list_noT[,grep(frac_ok, labs)]
  values2 <- apply(preds2,1, which.max)
  alphas2 <- apply(preds2,1, max)
  dom_com2 <- goods[values2]
  dominant_communities_2090_noT<-cbind(dominant_communities_2090_noT, dom_com2)
  dominant_communities_2090_noT<-cbind(dominant_communities_2090_noT, alphas2 )
  number_of_changes <- number_of_changes+as.numeric(values!=values1)
  area_no_change <- sum(as.numeric(values==values1)*111*111*weights_cos)/3607000
  area_change <- sum(as.numeric(values!=values1)*111*111*weights_cos)/3607000
  tot <- area_no_change+area_change
  area_no_change1 <- area_no_change*100/tot
  area_change1<-area_change*100/tot
  changes_metac<-append(changes_metac, c(frac_ok,area_change1, area_no_change1))
  count=1
  for (condition in basins_conds){
    area_no_change <- sum(as.numeric(values==values1)*111*111*weights_cos*as.numeric(condition)*pred_2006_list[,grep(frac_ok, labs)])
    area_change <- sum(as.numeric(values!=values1)*111*111*weights_cos*as.numeric(condition)*pred_2090_list[,grep(frac_ok, labs)])
    tot <- area_no_change+area_change
    area_no_change1 <- area_no_change*100/tot
    area_change1<-area_change*100/tot
    changes_metac_basins<-append(changes_metac_basins, c(frac_ok,basins[count],area_change1, area_no_change1))
    for (k in goods){
      area <- sum(as.numeric(dom_com0==k)*111*111*weights_cos*as.numeric(condition))
      area1 <- sum(as.numeric(dom_com1==k)*111*111*weights_cos*as.numeric(condition))
      if (length(area)>0 & area>2*10^6 & area1>2*10^6){
        areas <- rbind(areas, c(frac_ok, k, basins[count],area, area1))
        shifts <- rbind(shifts, DS[DS$Fraction==frac & DS$Metacommunity==k & DS$Location==basins[count],])
      }
    }
    count=count+1
  }
}
saveRDS(number_of_changes, 'number_of_changes.rds')
saveRDS(dominant_communities_2090, 'dominant_communities_2090.rds')
saveRDS(dominant_communities_2090_noT, 'dominant_communities_2090_noT.rds')
saveRDS(dominant_communities_2006, 'dominant_communities_2006.rds')
saveRDS(dominant_communities_woa, 'dominant_communities_woa.rds')
areas<- as.data.frame(areas)
areas$V4<-as.numeric(levels(areas$V4))[areas$V4]
areas$V5<-as.numeric(levels(areas$V5))[areas$V5]

saveRDS(shifts, 'shifts_dominant.rds')

area_numb_changes<-NULL
tot_area = sum(111*111*weights_cos)
for (j in 1:6){
  area_change <- sum(111*111*weights_cos*as.numeric(number_of_changes>=j))
  area_numb_changes<-rbind(area_numb_changes, c(j,area_change/tot_area))
}
write.table(area_numb_changes, 'area_numb_changes.txt')
area_numb_changes<-read.table('area_numb_changes.txt')

data_contour <- matrix(NA, ncol=150, nrow=360)
for (i in 1:dim(bray_curtis)[1]){
  dc=number_of_changes[i]
  data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= dc
}
color_pal = colorRampPalette(c('blue','yellow','red', 'darkred'))(8)
pdf(family="Helvetica",file=paste('number_of_changes.pdf', sep=''),width=10,height=4.065)
par(mar=c(0,0,0,0))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
.filled.contour(longi_sorted, lati_sorted,data_contour,  col=color_pal,levels=c(0:7))
legend(x = 75, y = 60, legend = c(0:6), col=color_pal,pch=15,cex=8, ncol=2)
# contour(x=longi_sorted,y=lati_sorted, z=data_contour,col='white',add=T, lwd=0.4,drawlabels=F)
hide_arctic()
axis_map0()
dev.off()

delta_dom <- function(x){
  return(x[order(x, decreasing = T)][1]-x[order(x, decreasing = T)][2])
}
mean_dom <- function(x){
  return(mean(x[order(x, decreasing = T)][1]-x[order(x, decreasing = T)][2:length(x)]))
}
fractions1 <- readRDS('fractions0.rds')
delta_vec_all <- NULL
delta_vec_all_w <- NULL
delta_vec_all_f <- NULL
mean_delta_vec_all <- NULL
alpha_vec_all <- NULL
for (frac in fractions1){
  delta_vec<-apply(pred_2006_list[,frac], 1, delta_dom)
  delta_vec_w<-apply(pred_woa_list[,frac], 1, delta_dom)
  delta_vec_f<-apply(pred_2090_list[,frac], 1, delta_dom)
  mean_vec <- apply(pred_2006_list[,frac], 1, mean_dom)
  max_vec <- apply(pred_2006_list[,frac], 1, max)
  delta_vec_all <- append(delta_vec_all, delta_vec)
  delta_vec_all_w <- append(delta_vec_all_w, delta_vec_w)
  delta_vec_all_f <- append(delta_vec_all_f, delta_vec_f)
  alpha_vec_all <- append(alpha_vec_all,max_vec)
  mean_delta_vec_all <- append(mean_delta_vec_all, mean_vec)
}
delta_vec_all<-matrix(delta_vec_all, ncol=6)
delta_vec_all <- delta_vec_all*100
delta_vec_all_w<-matrix(delta_vec_all_w, ncol=6)
delta_vec_all_w <- delta_vec_all_w*100
delta_vec_all_f<-matrix(delta_vec_all_f, ncol=6)
delta_vec_all_f <- delta_vec_all_f*100
mean_delta_vec_all<-matrix(mean_delta_vec_all, ncol=6)
mean_delta_vec_all <- mean_delta_vec_all*100
alpha_vec_all<-matrix(alpha_vec_all, ncol=6)
#col_vec <- colorRampPalette(c('red','orange','yellow','green', 'blue'))(101)
col_vec <- gray.colors(101)
col_vec1 <- col_vec[c(0:10)*10+1]
for (i in 1:6){
  data_contour_m <- matrix(NA, ncol=150, nrow=360)
  data_contour_w <- matrix(NA, ncol=150, nrow=360)
  data_contour_06 <- matrix(NA, ncol=150, nrow=360)
  data_contour_90 <- matrix(NA, ncol=150, nrow=360)
  for (j in 1:dim(bray_curtis)[1]){
    dc_m=mean_delta_vec_all[j,i]
    dc_w=delta_vec_all_w[j,i]
    dc_06=delta_vec_all[j,i]
    dc_90=delta_vec_all_f[j,i]
    data_contour_m[which(longi_sorted==bray_curtis$Long[j]),which(lati==bray_curtis$Lat[j])]= round(dc_m,-1)+1
    data_contour_w[which(longi_sorted==bray_curtis$Long[j]),which(lati==bray_curtis$Lat[j])]= round(dc_w,-1)+1
    data_contour_06[which(longi_sorted==bray_curtis$Long[j]),which(lati==bray_curtis$Lat[j])]= round(dc_06,-1)+1
    data_contour_90[which(longi_sorted==bray_curtis$Long[j]),which(lati==bray_curtis$Lat[j])]= round(dc_90,-1)+1
  }
  
  pdf(family="Helvetica",file=paste('delta_dom_',fractions[i],'_w.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour_w,  col=col_vec,levels=c(1:101))
  #points(x=bray_curtis$Long, y =bray_curtis$Lat, col=col_vec[round(delta_vec_all[,i])+1],pch=15,cex=0.229)
  SDMTools::legend.gradient(pnt,col_vec1,limits=c(0,1),title='Delta p',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  pdf(family="Helvetica",file=paste('delta_dom_',fractions[i],'_06.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour_06,  col=col_vec,levels=c(1:101))
  #points(x=bray_curtis$Long, y =bray_curtis$Lat, col=col_vec[round(delta_vec_all[,i])+1],pch=15,cex=0.229)
  SDMTools::legend.gradient(pnt,col_vec1,limits=c(0,1),title='Delta p',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  pdf(family="Helvetica",file=paste('delta_dom_',fractions[i],'_90.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour_90,  col=col_vec,levels=c(1:101))
  #points(x=bray_curtis$Long, y =bray_curtis$Lat, col=col_vec[round(delta_vec_all[,i])+1],pch=15,cex=0.229)
  SDMTools::legend.gradient(pnt,col_vec1,limits=c(0,1),title='Delta p',cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  pdf(family="Helvetica",file=paste('delta_dom_mean',fractions[i],'.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour_m,  col=col_vec,levels=c(1:101))
  #points(x=bray_curtis$Long, y =bray_curtis$Lat, col=col_vec[round(mean_delta_vec_all[,i])+1],pch=15,cex=0.229)
  SDMTools::legend.gradient(pnt,col_vec1,limits=c(0,1),title='Delta p',cex=0.4)
  hide_arctic()
  axis_map0() 
  dev.off()
}

code_letter <- c('F', 'E', 'D', 'C', 'B', 'A')
for (frac in 1:6){
  letter <- code_letter[frac]
  good <- best_models$Gen[best_models$Fraction==fractions0[frac]]
  covered_areas <- NULL
  for (k in good){
    cond0 <- dominant_communities_woa[,2*frac-1] == k
    cond <- dominant_communities_2006[,2*frac-1] == k
    cond1 <- dominant_communities_2090[,2*frac-1]== k
    area0 <- sum(cond0*111*111*weights_cos*dominant_communities_woa[,2*frac])
    area <- sum(cond*111*111*weights_cos*dominant_communities_2006[,2*frac])#area in km^2
    area1 <- sum(cond1*111*111*weights_cos*dominant_communities_2090[,2*frac])#area in km^2
    delta <- area1 - area     
    covered_areas <- rbind(covered_areas, c(k,  fractions[frac], model0, area0,area, area1, delta))
  }
  full_area <- sum(111*111*weights_cos)
  write(full_area, paste(model0,"_full_area.txt", sep=''))
  #colnames(covered_areas)<- c('Metacommunity','fraction', 'model','Area woa' ,'Area 2006', 'Area 2090', 'Delta area (2090-2006)')
  write.table(covered_areas, paste(model0,'_',fractions[frac],'_covered_areas.txt', sep=''), sep="\t", row.names = F, col.names = F)
  
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_woa.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  color_code0 <- scales::alpha(color_sets[[frac]][dominant_communities_woa[,2*frac-1]],
                               alpha = dominant_communities_woa[,2*frac])
  points(x=bray_curtis$Long, y =bray_curtis$Lat, col=color_code0,pch=15,cex=0.229)
  legend(x = 75, y = 60, legend = paste(letter,good, sep=''), col=color_sets[[frac]][good],pch=15,cex=8, ncol=2)
  text(x= 80, y = 65, labels = 'Metacommunity' ,cex=6)
  hide_arctic()
  axis_map0()
  dev.off()
  
  good1<-c(good,15)
  data_contour <- matrix(NA, ncol=150, nrow=360)
  alpha_contour <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(bray_curtis)[1]){
    dc=dominant_communities_woa[,2*frac-1][i]
    data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= dc
    alpha_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= round(dominant_communities_woa[,2*frac][i]*100)+1
  }
#   for (i in good){
#     data_contour[[i]][is.na(data_contour[[i]])]<-0
#   }
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_woa_contour.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
#   set_col=color_sets[[frac]]
#   ec=scales::alpha('white',0)
#   for (i in good){
#     .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
#   }
  .filled.contour(longi_sorted, lati_sorted,data_contour,  col=color_sets[[frac]][good],levels=good1)
  contour(x=longi_sorted,y=lati_sorted, z=data_contour,col='white',add=T, lwd=0.4,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
  
  png("contour.png", width=10,height=4.065, res=200, units = 'in', family ='Helvetica')
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()

  set_al <- NULL
  for (i in seq(0,1,0.01)){
    set_al <-append(set_al, rgb(i,0,0))
  }
  png('alpha.png', width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
  .filled.contour(longi_sorted, lati, alpha_contour,
                  col=set_al, levels = 1:101)
  hide_arctic()
  dev.off()
  
  im <- readPNG('contour.png')
  alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_woa_contour_alpha.pdf', sep=''),width=10,height=4.065)
  grid.raster(w)
  dev.off()

  png("contour.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)

  grid <-  expand.grid('lat'=seq(-89.5, 59.5,1), 'lon'=seq(-179.5, 179.5,1))
  grid$cell <- paste(grid$lat, grid$lon, sep='_')

  cond = bray_curtis$cell[delta_vec_all_w[,frac]<50]
  cond1 = grid$cell %in% cond
  
  co <- 1
  for (i in unique(grid$lon)){
    v <- sum(grid$lon==i)
    if (co%%4==1 ){
      cond1[grid$lon==i][c(1:v)[!(seq(1,v,1) %in% seq(1, v,4))]] <- F
    } else if(co%%4==3){
      cond1[grid$lon==i][c(1:v)[!(seq(1,v,1) %in% seq(3, v,4))]] <- F
    } else {
      cond1[grid$lon==i] <- F
    } 
    co <- co+1
  }
  points(x=grid$lon[cond1], y =grid$lat[cond1], col='black',pch=20,cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()

  im <- readPNG('contour.png')
  alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_woa_contour_alpha_1.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()

  png("contour0.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
  
  for (u in BGCP_shp@polygons){
    for (v in u@Polygons){
      polygon(v@coords, lwd = 8)
    }
  }
  hide_arctic()
  dev.off()
  
  im <- readPNG('contour0.png')
  # alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_woa_contour_alpha_longhurst.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()

  png("contour2.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
  biome_t <- biome
  for (q in unique(biome_t[biome_t!=0])){
    hi<-biome_t
    hi[biome!=q]=0
    contour(lo, lt, hi, col = 'black',add=T, drawlabels = F, cex=0.284)
  }
  plot(coastline,lwd=0.0475, col='black', add=T)
  hide_arctic()
  axis_map0()
  dev.off()
  
  im <- readPNG('contour2.png')
  # alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_woa_contour_alpha_McKingley.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()
  
  for (cou in 1:3){
    png("contour2.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
    par(mar=c(0,0,0,0))
    maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
    .filled.contour(longi_sorted, lati_sorted,data_contour,
                    col=color_sets[[frac]][good], levels=good1)
    contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
    biogeog <- t(to_plot[[cou]])
    for (q in unique(biogeog[biogeog !=0])){
      hi<-biogeog
      hi[biogeog!=q]=0
      contour(lo, lt, hi, col = 'black',add=T, drawlabels = F, cex=0.284)
    }
    plot(coastline,lwd=0.0475, col='black', add=T)
    hide_arctic()
    axis_map0()
    dev.off()
    
    im <- readPNG('contour2.png')
    w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
    pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],
                   '_metacommunities_woa_contour_alpha_',names[cou],'.pdf', sep=''),
        width=10,height=4.065)
    grid.raster(w)
    dev.off()
  }
  
  
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2006-15.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  color_code0 <- scales::alpha(color_sets[[frac]][dominant_communities_2006[,2*frac-1]],
                               alpha = dominant_communities_2006[,2*frac])
  points(x=bray_curtis$Long, y =bray_curtis$Lat, col=color_code0,pch=15,cex=0.229)
  legend(x = 75, y = 60, legend = good, col=color_sets[[frac]][good],pch=15,cex=8, ncol=2)
  text(x= 80, y = 65, labels = 'Metacommunity' ,cex=6)
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- matrix(NA, ncol=150, nrow=360)
  alpha_contour <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(bray_curtis)[1]){
    dc=dominant_communities_2006[,2*frac-1][i]
    data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= dc
    alpha_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= round(dominant_communities_2006[,2*frac][i]*100)+1
  }

  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2006-15_contour.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,  col=color_sets[[frac]][good],levels=good1)
  contour(x=longi_sorted,y=lati_sorted, z=data_contour,col='white',add=T, lwd=0.4,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
  
  png("contour.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
 
  
  set_al <- NULL
  for (i in seq(0,1,0.01)){
    set_al <-append(set_al, rgb(i,0,0))
  }
  png('alpha.png', width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
  .filled.contour(longi_sorted, lati, alpha_contour,
                  col=set_al, levels = 1:101)
  hide_arctic()
  dev.off()
  
  im <- readPNG('contour.png')
  alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2006-15_contour_alpha.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()

  png("contour.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)

  cond = bray_curtis$cell[delta_vec_all[,frac]<50]
  cond1 = grid$cell %in% cond
  
  co <- 1
  for (i in unique(grid$lon)){
    v <- sum(grid$lon==i)
    if (co%%4==1 ){
      cond1[grid$lon==i][c(1:v)[!(seq(1,v,1) %in% seq(1, v,4))]] <- F
    } else if(co%%4==3){
      cond1[grid$lon==i][c(1:v)[!(seq(1,v,1) %in% seq(3, v,4))]] <- F
    } else {
      cond1[grid$lon==i] <- F
    } 
    co <- co+1
  }
  points(x=grid$lon[cond1], y =grid$lat[cond1], col='black',pch=20,cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  im <- readPNG('contour.png')
  alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2006-15_contour_alpha_1.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()
  
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2090-99.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  color_code0 <- scales::alpha(color_sets[[frac]][dominant_communities_2090[,2*frac-1]],
                               alpha = dominant_communities_2090[,2*frac])
  points(x=bray_curtis$Long, y =bray_curtis$Lat, col=color_code0,pch=15,cex=0.229)
  legend(x = 75, y = 60, legend = good, col=color_sets[[frac]][good],pch=15,cex=8, ncol=2)
  text(x= 80, y = 65, labels = 'Metacommunity' ,cex=6)
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- matrix(NA, ncol=150, nrow=360)
  alpha_contour <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(bray_curtis)[1]){
    dc=dominant_communities_2090[,2*frac-1][i]
    data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= dc
    alpha_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= round(dominant_communities_2090[,2*frac][i]*100)+1
  }

  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2090-99_contour.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,  col=color_sets[[frac]][good],levels=good1)
  contour(x=longi_sorted,y=lati_sorted, z=data_contour,col='white',add=T, lwd=0.4,drawlabels=F)
  axis_map0()
  dev.off()
  
  png("contour.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()

  
  set_al <- NULL
  for (i in seq(0,1,0.01)){
    set_al <-append(set_al, rgb(i,0,0))
  }
  png('alpha.png', width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
  .filled.contour(longi_sorted, lati, alpha_contour,
                  col=set_al, levels = 1:101)
  dev.off()
  
  im <- readPNG('contour.png')
  alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2090-99_contour_alpha.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()

  png("contour.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)

  cond = bray_curtis$cell[delta_vec_all_f[,frac]<50]
  cond1 = grid$cell %in% cond
  
  co <- 1
  for (i in unique(grid$lon)){
    v <- sum(grid$lon==i)
    if (co%%4==1 ){
      cond1[grid$lon==i][c(1:v)[!(seq(1,v,1) %in% seq(1, v,4))]] <- F
    } else if(co%%4==3){
      cond1[grid$lon==i][c(1:v)[!(seq(1,v,1) %in% seq(3, v,4))]] <- F
    } else {
      cond1[grid$lon==i] <- F
    } 
    co <- co+1
  }
  points(x=grid$lon[cond1], y =grid$lat[cond1], col='black',pch=20,cex=0.4)
  hide_arctic()
  axis_map0()
  dev.off()
  
  im <- readPNG('contour.png')
  alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2090-99_contour_alpha_1.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()
  
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2090-99_noT.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  color_code0 <- scales::alpha(color_sets[[frac]][dominant_communities_2090_noT[,2*frac-1]],
                               alpha = dominant_communities_2090_noT[,2*frac])
  points(x=bray_curtis$Long, y =bray_curtis$Lat, col=color_code0,pch=15,cex=0.229)
  legend(x = 75, y = 60, legend = good, col=color_sets[[frac]][good],pch=15,cex=8, ncol=2)
  text(x= 80, y = 65, labels = 'Metacommunity' ,cex=6)
  hide_arctic()
  axis_map0()
  dev.off()
  
  data_contour <- matrix(NA, ncol=150, nrow=360)
  alpha_contour <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(bray_curtis)[1]){
    dc=dominant_communities_2090_noT[,2*frac-1][i]
    data_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= dc
    alpha_contour[which(longi_sorted==bray_curtis$Long[i]),which(lati==bray_curtis$Lat[i])]= round(dominant_communities_2090_noT[,2*frac][i]*100)+1
  }
#   for (i in good){
#     data_contour[[i]][is.na(data_contour[[i]])]<-0
#   }
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2090-99_noT_contour.pdf', sep=''),
      width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
#   set_col=color_sets[[frac]]
#   ec=scales::alpha('white',0)
#   for (i in good){
#     .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
#   }
  .filled.contour(x=longi_sorted,y=lati_sorted, z=data_contour,  col=color_sets[[frac]][good], levels = good1)
  contour(x=longi_sorted,y=lati_sorted, z=data_contour,col='white',add=T, lwd=0.4,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
  
  png("contour.png", width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
#   set_col=color_sets[[frac]]
#   ec=scales::alpha('white',0)
#   for (i in good){
#     .filled.contour(longi_sorted, lati_sorted,data_contour[[i]],  col=c(ec,set_col[i], set_col[i], set_col[i]), levels = c(0,0.1,i,15))
#   }
  .filled.contour(longi_sorted, lati_sorted,data_contour,
                  col=color_sets[[frac]][good], levels=good1)
  contour(longi_sorted, lati_sorted,data_contour,col='white',add=T, cex=0.284,drawlabels=F)
  hide_arctic()
  axis_map0()
  dev.off()
  
  set_al <- NULL
  for (i in seq(0,1,0.01)){
    set_al <-append(set_al, rgb(i,0,0))
  }
  png('alpha.png', width=10,height=4.065, res=200, units = 'in', family="Helvetica")
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="white",border="white",xpd=TRUE)
  .filled.contour(longi_sorted, lati, alpha_contour,
                  col=set_al, levels = 1:101)
  dev.off()
  
  im <- readPNG('contour.png')
  alph <- readPNG('alpha.png')
  w <- matrix(rgb(im[,,1], im[,,2], im[,,3], alph[,,1]), nrow = dim(im)[1], ncol = dim(im)[2])
  pdf(family="Helvetica",file=paste('projection_',model0,'_',fractions[frac],'_metacommunities_2090-99_noT_contour_alpha.pdf', sep=''),
      width=10,height=4.065)
  grid.raster(w)
  dev.off()
}

changes_metac <- t(matrix(changes_metac, ncol=6))
write.table(changes_metac, paste(model0, 'changes_metac.txt', sep='_'))
#write.table(areas, paste(model0, 'changes_metac_basins_frac.txt', sep='_'))
changes_metac_basins <- matrix(changes_metac_basins, ncol=4, byrow = T)
write.table(changes_metac_basins, paste(model0, 'changes_metac_basins_frac.txt', sep='_'))

changes_metac_basins<-read.table('model-mean_changes_metac_basins_frac.txt')
changes_metac <- read.table('model-mean_changes_metac.txt')

changes_metac_basins_0<-NULL
count=1
for (condition in basins_conds){
  area_no_change <- sum(111*111*weights_cos*as.numeric(condition)*as.numeric(number_of_changes==0))
  area_change <- sum(111*111*weights_cos*as.numeric(condition)*as.numeric(number_of_changes>0))
  tot <- area_no_change+area_change
  area_no_change1 <- area_no_change*100/tot
  area_change1<-area_change*100/tot
  changes_metac_basins_0<-append(changes_metac_basins_0, c(basins[count],area_change1, area_no_change1))
  count=count+1
}
changes_metac_basins_0 <- t(matrix(changes_metac_basins_0, ncol=5))
write.table(changes_metac_basins_0,paste(model0, 'changes_metac_basins.txt', sep='_'))
changes_metac_basins_0 <- read.table('model-mean_changes_metac_basins.txt')

area_change<-sum(as.numeric(number_of_changes>0)*111*111*weights_cos)
area_no_change<-sum(as.numeric(number_of_changes==0)*111*111*weights_cos)
tot=area_no_change+area_change
area_no_change1 <- area_no_change*100/tot
area_change1<-area_change*100/tot
write(c(model0, area_change1), 'changes_metac_area_total.txt', append=F)  

changes_metac_numb <- NULL
for (i in 0:5){
  area_change<-sum(as.numeric(number_of_changes>i)*111*111*weights_cos)*100/tot
  print(area_change)
}



