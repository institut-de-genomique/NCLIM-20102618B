#!bin/usr/bin/env Rscript
# written by paul frémont
library("gbm")
library("dismo")
library("FactoMineR")
# library("factoextra")
# library("readxl")
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
setwd("/env/export/cns_n02_scratch/scratch_TaraOcean/Niches_models_Neogen/final_model_5/")
source('vioplot.R')
source('hide_arctic.R')
source('axis_map0.R')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')

cambria <- readRDS('cambria.rds')
DS <- readRDS('shifts_dominant.rds')
colnames(DS)<- c("Fraction" ,"Metacommunity" , "Model","Location","Lat_2006", 'Long_2006',"Lat_2090", 'Long_2090',"Shift" ,"Lat_shift" ,"Long_shift")
DS1 <- read.table('Distances_shifts_model-mean.txt')
colnames(DS1)<- c("Fraction" ,"Metacommunity" , "Model","Location","Lat_2006", 'Long_2006',"Lat_2090", 'Long_2090',"Shift" ,"Lat_shift" ,"Long_shift")
CA <- read.table('Projections/covered_areas_all.txt')
colnames(CA) <- c("Metacommunity" , "Fraction",	"Model","Area woa","Area 2006"	,"Area 2090",	"Delta area (2090-2006)")
RL_w <- read.table('Projections/relative_influences_all.txt', header=F)
colnames(RL_w) <- c( 'Inf.', 'Metacommunity', 'Frac', 'Param')

d_long <- density(abs(DS$Long_shift))
d_lat <- density(abs(DS$Lat_shift))

x1=abs(DS$Long_shift)
x2=abs(DS$Lat_shift)
x3=abs(DS$Shift)
pdf(family="Helvetica",'shifts_distributions.pdf', width = 20, height = 10)
par(mar=c(5.1,5.1,4.1,2.1))
hist(x1, breaks=seq(0,5000, 50),col=scales::alpha("blue", 0.5),xlab='Shift (km)',
     main='Latitudinal vs Longitudinal shifts', 
     border=scales::alpha('blue', 0.5), ylab='Number', cex.lab=2, cex.axis=2,cex.main=2)
hist(x2, breaks=seq(0,5000, 50), add=T,col=scales::alpha("red", 0.5), 
     border=scales::alpha('red', 0.5))

#hist(x3[x1>x2], breaks=seq(0,5000, 50),col=scales::alpha("blue", 0.5),xlab='Shift (km)',ylab='Number', 
 #    border=scales::alpha('blue', 0.5), ylim=c(0,9), main='Total shifts: Latitudinal vs Longitudinal',
 #    cex.lab=2, cex.axis=2,cex.main=2)
#hist(x3[x1<x2], breaks=seq(0,5000, 50), add=T,col=scales::alpha("red", 0.5), 
 #    border=scales::alpha('red', 0.5))

#hist(x1[x1>x2], breaks=seq(0,5000, 50),col=scales::alpha("blue", 0.5),xlab='Shift (km)',ylab='Number', 
 #    border=scales::alpha('blue', 0.5), ylim=c(0,10), 
 #    main='Directional shifts: Latitudinal vs Longitudinal', cex.lab=2, cex.axis=2,cex.main=2)
#hist(x2[x1<x2], breaks=seq(0,5000, 50), add=T,col=scales::alpha("red", 0.5), 
 #    border=scales::alpha('red', 0.5))

#hist(abs(x1-x2)[x1>x2], breaks=seq(0,5000, 50),col=scales::alpha("blue", 0.5),xlab='Shift (km)',ylab='Number', 
 #    border=scales::alpha('blue', 0.5), ylim=c(0,15), 
 #   main='Directional shifts: Latitudinal - Longitudinal', cex.lab=2, cex.axis=2,cex.main=2)
#hist(abs(x1-x2)[x1<x2], breaks=seq(0,5000, 50), add=T,col=scales::alpha("red", 0.5), 
 #    border=scales::alpha('red', 0.5))
dev.off()

df <- readRDS('Genocenoses_env_parameters_woa.rds')

Fractions = c('180-2000', '20-180', '43952', '0.8-5', '0.22-3', '0-0.2')
Fractions1 = c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
CA$percentage <- CA$`Delta area (2090-2006)`/CA$`Area 2006`

shift_north_a_com_5 <- DS$Shift[DS$Fraction=='180-2000' & DS$Metacommunity==5 & DS$Location=='North_A']
write(shift_north_a_com_5, 'Shift_North_A_180-2000_5.txt')

trops <- c(8,5,6,11,5,3)
tropical_com <- NULL
temps <- c(5,6,3,8,7,7)
temp_com <- NULL
for (i in 1:6){
  cond = CA$Fraction==Fractions1[i] & CA$Metacommunity==trops[i]
  cond_t <- CA$Fraction==Fractions1[i] & CA$Metacommunity==temps[i]
  tropical_com<-rbind(tropical_com, c(Fractions1[i] ,CA$percentage[cond],
                                      CA$`Area 2006`[cond], CA$`Area 2090`[cond],
                                      CA$`Delta area (2090-2006)`[cond]))
  temp_com<-rbind(temp_com, c(Fractions1[i], CA$percentage[cond_t],
                              CA$`Area 2006`[cond_t], CA$`Area 2090`[cond_t],
                              CA$`Delta area (2090-2006)`[cond_t]))
}
me <- mean(as.numeric(temp_com[,2]))
me_06 <- mean(as.numeric(temp_com[,3]))
me_90 <- mean(as.numeric(temp_com[,4]))
me_d <- mean(as.numeric(temp_com[,5]))
me1 <- mean(as.numeric(tropical_com[,2]))
me1_06 <- mean(as.numeric(tropical_com[,3]))
me1_90 <- mean(as.numeric(tropical_com[,4]))
me1_d <- mean(as.numeric(tropical_com[,5]))
temp_com <- rbind(temp_com, c('average', me, me_06, me_90,me_d ))
tropical_com <- rbind(tropical_com, c('average', me1, me1_06, me1_90,me1_d ))
write.table(tropical_com, 'changes_trop_com.txt')
write.table(temp_com, 'changes_temp_com.txt')

for (frac in unique(RL_w$Frac)){
  print(frac)
  print(mean(RL_w$Inf.[RL_w$Param=='T' & RL_w$Frac==frac]))
  print(median(RL_w$Inf.[RL_w$Param=='T' & RL_w$Frac==frac]))
}
cols<-c('red', 'green', 'blue')
c<-1
plot(c(10,10), ylim=c(0,0.003), xlim=c(0,max(DS$Shift)))
for (frac in unique(DS$Fraction)[c(1,4,5)]){
  d <- density(DS$Shift[DS$Fraction==frac])
  polygon(d, col=scales::alpha(cols[c], 0.3))
  c <-c+1
}
dev.off()

shifts_stats<-NULL
for (frac in Fractions){
  shifts_stats <- rbind(shifts_stats, c(frac, mean(DS$Shift[DS$Fraction==frac]), median(DS$Shift[DS$Fraction==frac])))
}
write.table(shifts_stats, 'shifts_stats_fractions.txt', col.names = F)
perc_poleward=(sum(DS$Lat_shift[DS$Location=='North_A' | DS$Location=='North_P']>0)+sum(DS$Lat_shift[DS$Location=='South_A' | DS$Location=='South_P' | DS$Location=='Indian']<0 ))/dim(DS)[1]
polewards <- DS$Lat_shift[DS$Location=='South_A' | DS$Location=='South_P' | DS$Location=='Indian'][DS$Lat_shift[DS$Location=='South_A' | DS$Location=='South_P' | DS$Location=='Indian']<0]
polewards <- append(polewards, DS$Lat_shift[DS$Location=='North_A' | DS$Location=='North_P'][DS$Lat_shift[DS$Location=='North_A' | DS$Location=='North_P']>0])
data_to_write <- c('perc_poleward', 'mean_shift', 'median_shift', 'mean_poleward', 'median_poleward', 'max')
data_to_write <- rbind(data_to_write,c(perc_poleward, mean(DS$Shift),median(DS$Shift), mean(abs(polewards)), median(abs(polewards)), max(DS$Shift, na.rm = T)))
write.table(data_to_write, 'shifts_stats_pole_mean_med.txt', col.names = F, row.names = F)

shifts_stats<-NULL
for (frac in Fractions){
  shifts_stats <- rbind(shifts_stats, c(frac, mean(DS1$Shift[DS1$Fraction==frac]), median(DS1$Shift[DS1$Fraction==frac])))
}
write.table(shifts_stats, 'shifts_stats_fractions_all.txt', col.names = F)
perc_poleward=(sum(DS1$Lat_shift[DS1$Location=='North_A' | DS1$Location=='North_P']>0)+sum(DS1$Lat_shift[DS1$Location=='South_A' | DS1$Location=='South_P' | DS1$Location=='Indian']<0 ))/dim(DS1)[1]
polewards <- DS1$Lat_shift[DS1$Location=='South_A' | DS1$Location=='South_P' | DS1$Location=='Indian'][DS1$Lat_shift[DS1$Location=='South_A' | DS1$Location=='South_P' | DS1$Location=='Indian']<0]
polewards <- append(polewards, DS1$Lat_shift[DS1$Location=='North_A' | DS1$Location=='North_P'][DS1$Lat_shift[DS1$Location=='North_A' | DS1$Location=='North_P']>0])
data_to_write <- c('perc_poleward', 'mean_shift', 'median_shift', 'mean_poleward', 'median_poleward', 'max')
data_to_write <- rbind(data_to_write,c(perc_poleward, mean(DS1$Shift),median(DS1$Shift), mean(abs(polewards)), median(abs(polewards)), max(DS1$Shift, na.rm = T)))
write.table(data_to_write, 'shifts_stats_pole_mean_med_all.txt', col.names = F, row.names = F)


pdf(family="Helvetica",'Relative_influences.pdf', height = 10, width = 15 )
boxplot(Inf.~Param, data=RL_w, ylab='Relative influence')
dev.off()
for (frac in unique(RL_w$Frac)){
  pdf(family="Helvetica",paste('Relative_influences_',frac,'.pdf', sep=''), height = 10, width = 15 )
  boxplot(Inf.~Param, data=RL_w[RL_w$Frac==frac,], ylab='Relative influence')
  dev.off()
}
pdf(family="Helvetica",'distribution_shifts.pdf')
test<-density(DS$Shift)
plot(test, xlim=c(0, max(DS$Shift)), xaxs='i', xlab='Shift (kms)', main='Shifts distribution')
dev.off()

pairwise.wilcox.test(RL_w$Inf., g = RL_w$Param)

x1 <- RL_w$Inf.[RL_w$Param=='T']
x2 <- RL_w$Inf.[RL_w$Param=='Sal']
x3 <- RL_w$Inf.[RL_w$Param=='Si']
x4 <- RL_w$Inf.[RL_w$Param=='NO3']
x5 <- RL_w$Inf.[RL_w$Param=='Phos']
x6 <- RL_w$Inf.[RL_w$Param=='Fe']
x8 <- RL_w$Inf.[RL_w$Param=='SI_NO3']
colores = brewer.pal(8, "Dark2")
pdf(family="Helvetica",'relative_influence_violin.pdf')
vioplot(x1, x2, x3,x4,x5,x6,x8, names=c('T', 'Sal', 'Si', 'NO3', 'Phos', 'Fe',  'SI_NO3'), col =colores)
dev.off()
for (frac in unique(RL_w$Frac)){
  x1 <- RL_w$Inf.[RL_w$Param=='T' & RL_w$Frac==frac]
  x2 <- RL_w$Inf.[RL_w$Param=='Sal' & RL_w$Frac==frac]
  x3 <- RL_w$Inf.[RL_w$Param=='Si' & RL_w$Frac==frac]
  x4 <- RL_w$Inf.[RL_w$Param=='NO3' & RL_w$Frac==frac]
  x5 <- RL_w$Inf.[RL_w$Param=='Phos' & RL_w$Frac==frac]
  x6 <- RL_w$Inf.[RL_w$Param=='Fe' & RL_w$Frac==frac]
  x8 <- RL_w$Inf.[RL_w$Param=='SI_NO3' & RL_w$Frac==frac]
  colores = brewer.pal(8, "Dark2")
  pdf(family="Helvetica",paste('Relative_influences_violin_',frac,'.pdf', sep=''))
  vioplot(x1, x2, x3,x4,x5,x6,x8, names=c('T', 'Sal', 'Si', 'NO3', 'Phos', 'Fe',  'SI_NO3'), col =colores)
  dev.off()
}

Fractions = c('180-2000', '20-180', '43952', '0.8-5', '0.22-3', '0-0.2')
Fractions1 = c('180-2000', '20-180', '5-20', '0.8-5', '0.22-3', '0-0.2')
locations <- unique(DS$Location)
models = c('model-mean')
global_stats_area <- NULL
for (frac in Fractions1){
  pdf(family="Helvetica",file = paste('areas_', frac, '.pdf', sep =''), width=10,height=10)
  for (mod in models){
    if (frac != '0-0.2'){
      color = colorRampPalette(brewer.pal(11,"Spectral"))(max(CA$Metacommunity[CA$Fraction == frac & CA$Model == mod]))
    }
    else{
      color = colorRampPalette(brewer.pal(11,"Spectral"))(8)
    }
    full_area <- as.numeric(read.table(paste(mod,'_full_area.txt', sep='')))
    data <- as.matrix(CA[CA$Fraction == frac & CA$Model == mod,4:6])
    data <- rbind(data, c( (full_area-sum(CA[CA$Fraction == frac & CA$Model == mod,4])),(full_area-sum(CA[CA$Fraction == frac & CA$Model == mod,5])), (full_area-sum(CA[CA$Fraction == frac & CA$Model == mod,6]))))
    data <- rbind(data,c( 360700000-sum(data[,1]), 360700000-sum(data[,1]), 360700000-sum(data[,1])))
    barplot(data, col = c(color[CA$Metacommunity[CA$Fraction == frac & CA$Model == mod]], 'white', 'lightblue'), 
            legend = c(CA$Metacommunity[CA$Fraction == frac & CA$Model == mod], 'not covered' , 'arctic (no samples)'), ylab = 'Oceanic surface', main=mod,
            xlim=c(0, 7), args.legend=list(title="Metacommunities", bty= 'n', cex=1.5), cex.axis = 1.3, cex.lab =1.3, cex.names =1.3)
    data_percent <- data/3607000
    percent_covered_woa <- sum(data_percent[1:(dim(data)[1]-2),1])
    percent_covered_2006 <- sum(data_percent[1:(dim(data)[1]-2),2])
    percent_covered_2090 <- sum(data_percent[1:(dim(data)[1]-2),3])
    global_stats_area<-append(global_stats_area, c(mod,frac, percent_covered_woa,percent_covered_2006, percent_covered_2090))
  }
  dev.off()
}
global_stats_area<-t(matrix(global_stats_area, ncol=6))
write.table(global_stats_area, 'global_stats_area.txt', row.names = F)

# temperate communities
temps<-c(5, 6, 3, 8, 7,7)
trops<-c(8,5,6,11,5,3)
for (i in 1:6){
  print(Fractions1[i])
  print(CA$percentage[CA$Metacommunity== temps[i] & CA$Fraction==Fractions1[i]])
  print(CA$percentage[CA$Metacommunity== trops[i] & CA$Fraction==Fractions1[i]])
}

locs <- c('North Atlantic', 'South Atlantic', 'North Pacific', 'South Pacific', 'Indian')
positions =c(-5:5)
positions_0 = NULL
posits =seq(0.1, 1,0.1)
posits_0 =seq(1,0.1,-0.1)
for (i in c(-5:4)){
  if (i<0){
    positions_0 <- append(positions_0, i+posits_0)
  } else{
    positions_0 <- append(positions_0, i+posits)
  }
}


labels = as.character(c(-c(1:1 %o% 10^(4:0)),0, c(1:1 %o% 10^(0:4))))
labels_0 = c(-c(9:2 %o% 10^(4:0)), c(2:9 %o% 10^(0:4)))
plot(x=0, y=0, 
     xlim = c(-6,6), 
     ylim =c(-6,6),
     axes=F,xlab = 'Longitude (kms)', ylab ='Latitude (kms)') 
axis(side=1, at=positions, labels=labels)
axis(side=2, at=positions, labels=labels)
axis(side=1, at=positions_0, labels=rep(NA, length(positions_0)), tck=-0.005)
axis(side=2, at=positions_0, labels=rep(NA, length(positions_0)), tck=-0.005)


posits_1<- seq(1, 9,1)
posits_2 <- seq(9,1,-1)

posits_3<- seq(1, 10,0.1)
posits_4 <- seq(10,1,-0.1)

labels = as.character(c(-c(1:1 %o% 10^(4:1)),0, c(1:1 %o% 10^(1:4))))

positions_2=c(0:10)
positions_3 <- NULL
positions_4<- NULL
labels_3 <- NULL
values_3 <- NULL
for (j in c(1:9)){
    if (j>4){
      positions_3 <- append(positions_3, j+log10(posits_1)-1)
      positions_4 <- append(positions_4, j+log10(posits_3)-1)
      for (i in 1:9){
        if (i==1){
          print(j)
          labels_3 <- append(labels_3, labels[j])
          values_3 <- append(values_3, as.numeric(labels[j]))
        } else{
          labels_3 <- append(labels_3, NA)
          values_3 <- append(values_3, as.numeric(labels[j])*i*0.1)
        }
      }
    } else{
      positions_3 <- append(positions_3, j-log10(posits_2))
      positions_4 <- append(positions_4, j-log10(posits_4))
      for (i in 1:9){
        if (i==1){
          print(j)
          labels_3 <- append(labels_3, labels[j])
          values_3 <- append(values_3, as.numeric(labels[j]))
        } else{
          labels_3 <- append(labels_3, NA)
          values_3 <- append(values_3, as.numeric(labels[j])*(11-i)*0.1)
        }
      }
    }
    
}
values_4<- NULL
for (u in 4:1){
  values_4 <- append(values_4, seq(-10^u, -10^(u-1),  10^(u-2)))
}
for (u in 1:5){
  values_4 <- append(values_4, seq(10^(u-1), 10^u,  10^(u-2)))
}

plot(positions_3, rep(1, length(positions_3)))
plot(x=0, y=0, xlim = c(0,8), 
     ylim =c(0,8),
     axes=F,xlab = 'Longitude (kms)', ylab ='Latitude (kms)', col='white') 
axis(side=1, at=positions_3, labels=labels_3)
axis(side=2, at=positions_3, labels=labels_3)


col_code=c('red', 'blue', 'green')
colors1 <- c('saddlebrown', 'red', 'dodgerblue2', 'darkgreen', 'darkviolet', 'darkorange')
DS$col <- colors1[7-match(DS$Fraction,Fractions )]
for (loc in locations){
  pdf(family="Helvetica",file = paste('shifts_', loc, '_fraction.pdf', sep =''), width=15,height=15)
  for (mod in models){
    cond0=DS$Model == mod & DS$Location == loc & DS$Lat_shift>0 & DS$Long_shift>0
    cond1=DS$Model == mod & DS$Location == loc  & DS$Lat_shift>0 & DS$Long_shift<0
    cond2=DS$Model == mod & DS$Location == loc  & DS$Lat_shift<0 & DS$Long_shift>0
    cond3=DS$Model == mod & DS$Location == loc  & DS$Lat_shift<0 & DS$Long_shift<0
    conditions = list(cond0,cond1, cond2, cond3)
    
    cond4 = abs(DS$Lat_2006) < 30
    cond5 = abs(DS$Lat_2006) > 30 & abs(DS$Lat_2006) < 50
    cond6 = abs(DS$Lat_2006) > 50
    conditions_0 = list(cond4, cond5, cond6)
    
    cond7 = log10(DS$Lat_shift)<0 & log10(DS$Long_shift)>0
    cond7[is.na(cond7)]<-TRUE
    cond8 = log10(DS$Lat_shift)>0 & log10(DS$Long_shift)<0
    cond8[is.na(cond8)]<-TRUE
    cond9 = log10(-DS$Lat_shift)<0 & log10(DS$Long_shift)>0
    cond9[is.na(cond9)]<-TRUE
    cond10 = log10(-DS$Lat_shift)>0 & log10(DS$Long_shift)<0
    cond10[is.na(cond10)]<-TRUE
    cond11 = log10(DS$Lat_shift)<0 & log10(-DS$Long_shift)>0
    cond11[is.na(cond11)]<-TRUE
    cond12 = log10(DS$Lat_shift)>0 & log10(-DS$Long_shift)<0
    cond12[is.na(cond12)]<-TRUE
    cond13 = log10(-DS$Lat_shift)<0 & log10(-DS$Long_shift)>0
    cond13[is.na(cond13)]<-TRUE
    cond14 = log10(-DS$Lat_shift)>0 & log10(-DS$Long_shift)<0
    cond14[is.na(cond14)]<-TRUE
    t=sum(DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) < 30 , na.rm =T)
    t1=sum(DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 30 & abs(DS$Lat_2006) < 50 , na.rm = T)
    t2=sum(DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 50 , na.rm = T)
    ts = c(t,t1,t2)
#     plot(x=0, y=0, 
#          xlim = c(-6,6), 
#          ylim =c(-6,6),
#          axes=F,xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=mod) 
#     axis(side=1, at=positions, labels=labels)
#     axis(side=2, at=positions, labels=labels)
#     axis(side=1, at=positions_0, labels=rep(NA, length(positions_0)), tck=-0.005)
#     axis(side=2, at=positions_0, labels=rep(NA, length(positions_0)), tck=-0.005)
    plot(x=0, y=0, xlim = c(0,8), 
         ylim =c(0,8.5),
         axes=F,xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=mod, col='white', cex.lab=1.5) 
    axis(side=1, at=positions_3, labels=labels_3, cex.axis=2)
    axis(side=2, at=positions_3, labels=labels_3, cex.axis=2)

    
    x <- NULL
    y <- NULL
    x2 <- NULL
    y2 <- NULL
    color_vec <- NULL
    for (j in 1:3){
      for (u in 1:4){
        t = ts[j]
        if(!(is.na(t)) & t != 0){
          conda=conditions[[u]]
          condb=conditions_0[[j]]
          if (sum(conda*condb)>0){
            x <- append(x, rep(4,sum(conda*condb)))
            y <- append(y, rep(4,sum(conda*condb)))
            x2 <- append(x2,DS[conda & condb , 11])
            y2 <- append(y2, DS[conda & condb, 10])
            color_vec <- append(color_vec, DS[conda & condb, 12])
            #color_vec <- append(color_vec, col_code[j])
          }
        }
      }
    }
    x3<- x2
    y3<- y2
    for (count in 1:length(x2)){
      x2[count]=positions_4[which.min(abs(x3[count]-values_4))]
      y2[count]=positions_4[which.min(abs(y3[count]-values_4))]
    }
    arrows(x0 = x, y0 = y,x1 = x2, y1 = y2, col = color_vec, lwd=2.5)
    #legend(x = 'topleft', bty= 'n',
    #legend = c('Equatorial communities', 'Temperate communities', 'Polar communities' ), fill =c('red', 'blue', 'green'), cex = 1.5)
    legend(x='topright', bty='n', legend=Fractions[6:1], fill=colors1, cex=1.5, ncol = 3)
  }
  dev.off()
}


shift_map <- function(shifts, name_pdf){
  cases <- NULL
  for (u in 1:dim(shifts)[1]){
    if (shifts$Long_2006[u] < -100 & shifts$Long_2090[u] > 100 ){
      cases <-append(cases, 1)
      shifts$Long_shift[u] <- -shifts$Long_shift[u]
    } else if (shifts$Long_2090[u] < -20& shifts$Long_2006[u] > 100){
      cases <-append(cases, 2)
      shifts$Long_shift[u] <- -shifts$Long_shift[u] 
    } else{
      cases <-append(cases, 3)
    }
  }
  colors1 <- c('saddlebrown', 'red', 'dodgerblue2', 'darkgreen', 'darkviolet', 'darkorange')
  pdf(family="Helvetica",file=name_pdf,width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475,  col='black',add=T)
  cols <- colors1[7-match(shifts$Fraction, Fractions)]
  arrows(x0=shifts$Long_2006[cases==3], y0=shifts$Lat_2006[cases==3], x1=shifts$Long_2090[cases==3], 
         y1=shifts$Lat_2090[cases==3], lwd=1.7, col = cols[cases==3], length = 0.025)
  
  arrows(x0=shifts$Long_2006[cases==2], y0=shifts$Lat_2006[cases==2], x1=shifts$Long_2090[cases==2]+360, 
         y1=shifts$Lat_2090[cases==2], lwd=1.7, col = cols[cases==2], length = 0.025)
  arrows(x0=shifts$Long_2006[cases==2]-360, y0=shifts$Lat_2006[cases==2], x1=shifts$Long_2090[cases==2], 
         y1=shifts$Lat_2090[cases==2], lwd=1.7, col = cols[cases==2], length = 0.025)
  
  arrows(x0=shifts$Long_2006[cases==1], y0=shifts$Lat_2006[cases==1], x1=shifts$Long_2090[cases==1]-360, 
         y1=shifts$Lat_2090[cases==1], lwd=1.7, col = cols[cases==1], length = 0.025)
  arrows(x0=shifts$Long_2006[cases==1]+360, y0=shifts$Lat_2006[cases==1], x1=shifts$Long_2090[cases==1], 
         y1=shifts$Lat_2090[cases==1], lwd=1.7, col = cols[cases==1], length = 0.025)
  #legend(x = 75, y = 70, legend = fractions[6:1], col = colors1, cex = 8, pch=15)
  hide_arctic()
  axis_map0()
  dev.off()
}
shift_map(DS, 'shifts_map.pdf')
shift_map(DS1, 'shifts_map_all.pdf') 
# for (frac in Fractions){
#   colors = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(max(df$Genocenose[df$Fraction==frac], na.rm = T))
#   for (mod in models){
#     pdf(family="Helvetica",file = paste('shifts_', mod,'_',frac, '_com.pdf', sep =''), width=15,height=10)
#     par(mfrow=c(2,3))
#     for (loc in locations){
#       plot(x=0, y=0, 
#            xlim = c(-2000,2000), 
#            ylim =c(-2000,2000),
#            xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=locs[which(locations==loc)], cex.axis = 1.5, cex.lab = 1.5, cex.main =2)
#       t = sum(DS$Model == mod & DS$Location == loc & DS$Fraction==frac, na.rm =T)
#       if(!(is.na(t))){
#         arrows(x0 = rep(0, t+sum(is.na(t))),
#                y0 = rep(0,t+sum(is.na(t))),
#                x1 = DS[DS$Model == mod & DS$Location == loc & DS$Fraction==frac, 10], 
#                y1 = DS[DS$Model == mod & DS$Location == loc & DS$Fraction==frac, 11], col = colors[DS[DS$Model == mod & DS$Location == loc & DS$Fraction==frac, 2]])
#       }
#       legend(x = 'topleft', bty= 'n',
#              legend = unique(DS$Metacommunity[DS$Fraction==frac]), fill =colors[unique(DS$Metacommunity[DS$Fraction==frac])], cex = 1.5)
#     }
#     dev.off()
#   }
# }
# 
# locs <- c('North Atlantic', 'South Atlantic', 'North Pacific', 'South Pacific', 'Indian')
# for (mod in models){
#   pdf(family="Helvetica",file = paste('shifts_', mod, '.pdf', sep =''), width=15,height=10)
#   par(mfrow=c(2,3))
#   for (loc in locs){
#     plot(x=0, y=0, 
#          xlim = c(-2000, 2000), 
#          ylim =c(-2000, 2000),
#          xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=locs[which(locations==loc)], cex.axis = 1.5, cex.lab = 1.5, cex.main =2)
#     t = sum(DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) < 30, na.rm =T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS[DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) < 30, 10], 
#              y1 = DS[DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) < 30, 11], col = 'red')
#     }
#     t = sum(DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 30 & abs(DS$Lat_2006) < 50, na.rm = T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS[DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 30 & abs(DS$Lat_2006) < 50, 10], 
#              y1 = DS[DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 30 & abs(DS$Lat_2006) < 50, 11], col = 'blue')
#     }
#     t =  sum(DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 50, na.rm = T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0,t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS[DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 50, 10], 
#              y1 = DS[DS$Model == mod & DS$Location == loc & abs(DS$Lat_2006) > 50, 11], col = 'green')
#     }
#     legend(x = 'topleft', bty= 'n',
#            legend = c('Equatorial communities', 'Temperate communities', 'Polar communities' ), fill =c('red', 'blue', 'green'), cex = 1.5)
#   }
#   
#   dev.off()
# }
# 
# for (mod in models){
#   for (loc in locations){
#     print(loc)
#     for (frac in Fractions){
#     print(frac)
#     print(mean(DS$Shift[DS$Location == loc & DS$Model==mod & DS$Fraction ==frac], na.rm =T))
#     }
#   }
#   print('')
# }
# 
# colors0 <- brewer.pal(6, "Paired")
# #colors0[4] <- 'darkorange'
# for (mod in models){
#   pdf(family="Helvetica",file = paste('shifts_Fractions_', mod, '.pdf', sep =''), width=15,height=10)
#   par(mfrow=c(2,3))
#   for (loc in locations){
#     plot(x=0, y=0, 
#          xlim = c(-2000,2000), 
#          ylim =c(-2000,2000),
#          xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=locs[which(locations==loc)], cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
#     i = 1
#     for (frac in Fractions){
#       t = sum(DS$Model == mod & DS$Location == loc & DS$Fraction == frac, na.rm =T)
#       if(!(is.na(t))){
#         arrows(x0 = rep(0, t+sum(is.na(t))), 
#                y0 = rep(0,t+sum(is.na(t))),
#                x1 = DS[DS$Model == mod & DS$Location == loc & DS$Fraction == frac, 10], 
#                y1 = DS[DS$Model == mod & DS$Location == loc & DS$Fraction == frac, 11], col = colors0[i])
#         i = i + 1
#       }
#     }
#     legend(x = 'topleft', bty= 'n',
#            legend = Fractions1, fill = colors0, cex = 1.5)
#   }
#   dev.off()
# }
# 
# locations <- unique(DS0$Location)
# for (loc in locations){
#   pdf(family="Helvetica",file = paste('shifts_', loc, '_all.pdf', sep =''), width=18,height=5)
#   par(mfrow=c(1,4))
#   for (mod in models){
#     plot(x=0, y=0, 
#          xlim = c(min(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T)), 
#          ylim =c(min(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T)),
#          xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=mod, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
#     t = sum(DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) < 30, na.rm =T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) < 30, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) < 30, 11], col = 'red')
#     }
#     t = sum(DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 30 & abs(DS0$Lat_2006) < 50, na.rm = T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 30 & abs(DS0$Lat_2006) < 50, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 30 & abs(DS0$Lat_2006) < 50, 11], col = 'blue')
#     }
#     t =  sum(DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 50, na.rm = T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0,t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 50, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 50, 11], col = 'green')
#     }
#     legend(x = 'topleft', bty= 'n',
#            legend = c('Equatorial communities', 'Temperate communities', 'Polar communities' ), fill =c('red', 'blue', 'green'), cex = 1.5)
#   }
#   
#   dev.off()
# }
# 
# for (mod in models){
#   pdf(family="Helvetica",file = paste('shifts_', mod, '_all.pdf', sep =''), width=15,height=10)
#   par(mfrow=c(2,3))
#   for (loc in locations){
#     plot(x=0, y=0, 
#          xlim = c(min(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T)), 
#          ylim =c(min(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T)),
#          xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=locs[which(locations==loc)], cex.axis = 1.5, cex.lab = 1.5, cex.main =2)
#     t = sum(DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) < 30, na.rm =T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) < 30, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) < 30, 11], col = 'red')
#     }
#     t = sum(DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 30 & abs(DS0$Lat_2006) < 50, na.rm = T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 30 & abs(DS0$Lat_2006) < 50, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 30 & abs(DS0$Lat_2006) < 50, 11], col = 'blue')
#     }
#     t =  sum(DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 50, na.rm = T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0,t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 50, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) > 50, 11], col = 'green')
#     }
#     legend(x = 'topleft', bty= 'n',
#            legend = c('Equatorial communities', 'Temperate communities', 'Polar communities' ), fill =c('red', 'blue', 'green'), cex = 1.5)
#   }
#   
#   dev.off()
# }
# 
# 
# colors = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(max(DS0$Metacommunity, na.rm = T))
# for (mod in models){
#   pdf(family="Helvetica",file = paste('shifts_', mod, '_all_com.pdf', sep =''), width=15,height=10)
#   par(mfrow=c(2,3))
#   for (loc in locations){
#     plot(x=0, y=0, 
#          xlim = c(min(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T)), 
#          ylim =c(min(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T)),
#          xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=locs[which(locations==loc)], cex.axis = 1.5, cex.lab = 1.5, cex.main =2)
#     t = sum(DS0$Model == mod & DS0$Location == loc, na.rm =T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc, 11], col = colors[DS0[DS0$Model == mod & DS0$Location == loc, 2]])
#     }
#   }
#   dev.off()
# }
# 
# locations <- unique(DS0$Location)
# for (loc in locations){
#   pdf(family="Helvetica",file = paste('shifts_', loc, '_all_com.pdf', sep =''), width=18,height=5)
#   par(mfrow=c(1,4))
#   for (mod in models){
#     plot(x=0, y=0, 
#          xlim = c(min(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 10], na.rm = T)), 
#          ylim =c(min(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T),max(DS0[DS0$Model == mod & DS0$Location == loc, 11], na.rm = T)),
#          xlab = 'Longitude (kms)', ylab ='Latitude (kms)', main=mod, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
#     t = sum(DS0$Model == mod & DS0$Location == loc, na.rm =T)
#     if(!(is.na(t))){
#       arrows(x0 = rep(0, t+sum(is.na(t))), 
#              y0 = rep(0,t+sum(is.na(t))),
#              x1 = DS0[DS0$Model == mod & DS0$Location == loc & abs(DS0$Lat_2006) < 30, 10], 
#              y1 = DS0[DS0$Model == mod & DS0$Location == loc , 11], col = colors[DS0[DS0$Model == mod & DS0$Location == loc , 2]])
#     }
#   }
#   dev.off()
# }

# cosin_mat <- function(vec_1){
#   print(sqrt(t(vec_1) %*% vec_1))
#   cos <- apply(matrice,1,FUN = function(vec_2)
#     (t(vec_1) %*% vec_2)/ (sqrt(t(vec_1) %*% vec_1 * t(vec_2) %*% vec_2)) )
# }
# pdf(family="Helvetica",'shifts_correlation.pdf', width = 10, height = 10)
# par(mfrow=c(2,3))
# br <- c(seq(-1, 1, length.out = 100))
# for (loc in unique(DS$Location)){
#   shifts_loc <- DS[DS$Location==loc,]
#   levels(shifts_loc$Fraction)<- c(levels(shifts_loc$Fraction), '5-20')
#   shifts_loc$Fraction[shifts_loc$Fraction=='43952']<-'5-20'
#   shifts_loc$Long_2006[shifts_loc$Long_2006<0] <- shifts_loc$Long_2006[shifts_loc$Long_2006<0]
#   matrice <- as.matrix(shifts_loc[,5:8])
# #   row.names(matrice)<-NULL
# #   colnames(matrice)<-NULL
#   cosinus_mat <- apply(matrice, 1,FUN=cosin_mat)
# #   heatmap.2(cor(t(shifts_loc[,5:8])), trace="none",symm=TRUE, 
# #             labRow = paste(shifts_loc$Fraction, shifts_loc$Metacommunity, sep='_'),
# #             labCol = paste(shifts_loc$Fraction, shifts_loc$Metacommunity, sep='_'),
# #             dendrogram = "none", keysize=1,margins=c(10,22), col= greenred(99),breaks = br, 
# #             symkey = F, hclustfun = hclust, main=loc)
#   heatmap.2(cosinus_mat, trace="none",symm=TRUE, 
#             labRow = paste(shifts_loc$Fraction, shifts_loc$Metacommunity, sep='_'),
#             labCol = paste(shifts_loc$Fraction, shifts_loc$Metacommunity, sep='_'),
#             dendrogram = "none", keysize=1,margins=c(10,22), col= greenred(99),breaks = br, 
#             symkey = F, hclustfun = hclust, main=loc)
# }
# dev.off()
