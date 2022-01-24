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
#library("gplots")
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
library('grid')
library('matlab')
library('shape')

cambria <- readRDS('cambria.rds')
source('axis_map0.R')
source('hide_arctic.R')
model0 <-'model-mean'
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
bray_curtis <- readRDS('model-mean_bray_curtis.rds')
bc_new <-readRDS('model-mean_bray-curtis_dom.rds')
mapping <- readRDS('mapping_lo_lt.rds')
bray_curtis <- bray_curtis[mapping,]
lati<-seq(-89.5, 89.5, 1)
longi<-seq(0.5,359.5,1)
longi[181:360] =as.numeric(longi[181:360])-360
lati_sorted<-seq(-89.5, 60, 1)
longi_sorted =sort(longi)

coffs <- c(-0.1, 1/6)
for (cof in coffs){
  cond <- bc_new>cof
  cond1 <- (1-bray_curtis$b_c)>cof
  bc_cond <- bray_curtis[cond,]
  bc_cond1 <- bray_curtis[cond1,]
  bc_new1 <- bc_new[bc_new>cof]
  bc_prob <- 1-bray_curtis$b_c[(1-bray_curtis$b_c)>cof]
  values <- 100*bc_new1+1
  values_prob <- 100*bc_prob+1
  data_contour <- matrix(NA, ncol=150, nrow=360)
  data_contour1 <- matrix(NA, ncol=150, nrow=360)
  data_contour_prob <- matrix(NA, ncol=150, nrow=360)
  data_contour1_prob <- matrix(NA, ncol=150, nrow=360)
  for (i in 1:dim(bc_cond)[1]){
    data_contour[which(longi_sorted==bc_cond$Long[i]),which(lati==bc_cond$Lat[i])]= round(values[i])
  }
  for (i in 1:dim(bc_cond1)[1]){
    data_contour_prob[which(longi_sorted==bc_cond1$Long[i]),which(lati==bc_cond1$Lat[i])]= round(values_prob[i])
  }
  data_contour1 <- round(data_contour, -1)/10+1
  data_contour1_prob <- round(data_contour_prob, -1)/10+1
  OldMin <- 1
  OldMax <-11
  NewMin <- 1
  NewMax <-100
  data_contour1 <- (((data_contour1 - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  data_contour1_prob <- (((data_contour1_prob - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  set = colorRampPalette(c('blue','yellow', 'red'))(NewMax-1)
  set1 = colorRampPalette(c('blue','lightblue','green','yellow','orange' ,'red'))(NewMax-1)
  set2 = colorRampPalette(c('purple','blue','green','yellow' ,'red'))(NewMax-1)
  set3=jet.colors(NewMax-1)
  list_cols <- list(set1, set2,set3)
  colo <- c('set1', 'set2', 'set3')
  i=1
  for (set0 in list_cols){
    pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
    if (cof==-0.1){
      pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_dom_2006-2090_',colo[i] ,'.pdf', 
                   sep=''),width=10,height=4.065)
    } else{
      pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_dom_2006-2090_',colo[i] ,'_coff.pdf', 
                     sep=''),width=10,height=4.065)
    }
    par(mar=c(0,0,0,0))
    maps::map(database="world",fill=T,col="grey80",border="grey80",xpd=TRUE)
    plot(coastline,lwd=0.0475, col='black', add=T)
    .filled.contour(x=longi_sorted,y=lati_sorted, z=round(data_contour1),  col=set0, 
                    levels = c(1:NewMax))
    #contour(x=longi_sorted,y=lati_sorted, z=round(data_contour1),  col=scales::alpha('black',0.2), add=T)
    # points(x=bray_curtis$Long, y =bray_curtis$Lat, col=set[values],bg= set[values],pch=15,cex=0.229)
    #SDMTools::legend.gradient(pnt,set0,limits=c(0,1),title='Bray-curtis \ndissimilarity index',cex=0.4)
    Arrows(x0=-40,y0=15,x1=-40,y1=29.5, arr.type='triangle', col='red', lwd=6, arr.length=0.1, arr.width=0.2)
    Arrows(x0=-20,y0=-14.2,x1=-20,y1=-26, arr.type='triangle', col='red', lwd=6, arr.length=0.1, arr.width=0.2)
    hide_arctic()
    axis_map0()
    dev.off()
    
    if (cof==-0.1){
      pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_prob_2006-2090_',colo[i] ,'.pdf', 
                     sep=''),width=10,height=4.065)
    } else{
      pdf(family="Helvetica",file=paste('projection_',model0,'_metacommunities_bray-curtis_prob_2006-2090_',colo[i] ,'_coff.pdf', 
                     sep=''),width=10,height=4.065)
    }
    par(mar=c(0,0,0,0))
    maps::map(database="world",fill=T,col="grey80",border="grey80",xpd=TRUE)
    plot(coastline,lwd=0.0475, col='black', add=T)
    .filled.contour(x=longi_sorted,y=lati_sorted, z=round(data_contour1_prob),  col=set0, 
                    levels = c(1:NewMax))
    #contour(x=longi_sorted,y=lati_sorted, z=round(data_contour1),  col=scales::alpha('black',0.2), add=T)
    # points(x=bray_curtis$Long, y =bray_curtis$Lat, col=set[values],bg= set[values],pch=15,cex=0.229)
    #SDMTools::legend.gradient(pnt,set0,limits=c(0,1),title='Bray-curtis \ndissimilarity index',cex=0.4)
    hide_arctic()
    axis_map0()
    dev.off()
    i=i+1
  }
}
