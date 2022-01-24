#!bin/usr/bin/env Rscript
library('mapproj')
library('mapplots')
library('SDMTools')
library('RColorBrewer')
library('parallel')
library('StatMeasures')


cambria <- readRDS('cambria.rds')
# Ind_fishery <- read.csv('CatchInd_2010_2014.csv', header = T)
# saveRDS(Ind_fishery, "CatchInd_2010_2014.rds")
Ind_fishery <- readRDS("CatchInd_2010_2014.rds")
Ind_fishery$Reported[is.na(Ind_fishery$Reported)]<-0
Ind_fishery$IUU[is.na(Ind_fishery$IUU)]<-0
Ind_fishery$Discards[is.na(Ind_fishery$Discards)]<-0
Ind_fishery$Fish <- Ind_fishery$Reported + Ind_fishery$IUU + Ind_fishery$Discards
# Nind_fishery <- read.table('CatchNInd_2010_2014.csv', header = T, sep=',')
# saveRDS(Nind_fishery, "CatchNInd_2010_2014.rds")
Nind_fishery <- readRDS("CatchNInd_2010_2014.rds")
Nind_fishery$Reported[is.na(Nind_fishery$Reported)]<-0
Nind_fishery$IUU[is.na(Nind_fishery$IUU)]<-0
Nind_fishery$Discards[is.na(Nind_fishery$Discards)]<-0
Nind_fishery$Fish <- Nind_fishery$Reported + Nind_fishery$IUU + Nind_fishery$Discards
grid <- read.table('Cells.csv', header = T, sep=',')

sum_fish <- function(data){
  data$lat <- grid$Lat[match(data$Cell, grid$Cell)]
  decimals <- data$lat%%1 
  data$lat1 <- data$lat-decimals+0.5
  data$lon <- grid$Lon[match(data$Cell, grid$Cell)]
  decimals1 <- data$lon%%1 
  data$lon1 <- data$lon-decimals1+0.5
  data$cell_new <- as.character(paste(data$lat1, '_',data$lon1, sep=''))
  data_1 = data.frame(
    sum     = with(data, tapply(Fish, cell_new, sum))
  )
  data_1$cell <- rownames(data_1)
  data_2 <- data.frame(do.call('rbind', strsplit(as.character(data_1$cell),'_',fixed=TRUE)))
  data_2 <- cbind(data_2, data_1$sum)
  data_2 <- cbind(data_2,data_1$cell)
  rownames(data_2)<-NULL
  colnames(data_2)<- c('lat', 'lon', 'fish', 'cell')
  data_2$lat<-as.numeric(levels(data_2$lat))[data_2$lat]
  data_2$lon<-as.numeric(levels(data_2$lon))[data_2$lon]
  sum_fish<- data_2
}

Ind_fishery_new <- sum_fish(Ind_fishery)
Nind_fishery_new <- sum_fish(Nind_fishery)
values <- Nind_fishery_new$fish[match(Ind_fishery_new$cell,Nind_fishery_new$cell)]
values[is.na(values)]<-0
data_fisheries<- Ind_fishery_new
data_fisheries$fish <- data_fisheries$fish+values
col_code <- colorRampPalette(c('blue','green','yellow', 'red'))(10)
data_fisheries <- as.data.frame(data_fisheries)
colnames(data_fisheries)<- c('Lat', 'Long', 'Fish', 'Cell')
colors <- col_code[decile(as.numeric(data_fisheries$Fish))]
data_fisheries$deciles <- decile(as.numeric(data_fisheries$Fish))

saveRDS(data_fisheries, 'data_fisheries.rds')
data_fisheries<- readRDS('data_fisheries.rds')

pdf(family="Helvetica",'fisheries_map.pdf', width=10,height=4.065)
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
pnt=cbind(x =c(100,105,105,100), y =c(25,25,45,45))
points(x=data_fisheries$Long, y=data_fisheries$Lat, 
       col=colors,
       pch=15,cex=2.72)
SDMTools::legend.gradient(pnt,col_code,limits=c(round(min(data_fisheries$Fish, na.rm=T ), digits = 3),round(max(data_fisheries$Fish, na.rm=T ), digits = 3)),title='Fishing rate\n(kg/km^2/year)',cex=0.4)
dev.off()
  
