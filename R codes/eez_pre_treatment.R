library('rgdal')

cambria <- readRDS('cambria.rds')
shp='World_EEZ_v10_20180221/eez_boundaries_v10.shp'
myshp = readOGR(shp, layer = basename(strsplit(shp, "\\.")[[1]])[1])
data_eez= spTransform(myshp, CRS("+proj=longlat +datum=WGS84"))
test<- data_eez@data
# 
shp1='World_EEZ_v10_20180221/eez_v10.shp'
myshp1 = readOGR(shp1)
data_eez1= spTransform(myshp1, CRS("+proj=longlat +datum=WGS84"))
test1 <- data_eez1@data
test2 <- data_eez1@polygons
saveRDS(test2, 'eez_data.rds')
test2<-readRDS('eez_data.rds')
model0 = commandArgs(trailingOnly = T)
#bray_curtis<- readRDS(paste(model0,'_bray_curtis.rds', sep=''))

#if (model0=='model-mean'){
#  pdf(family="Helvetica",'EEZ_map.pdf', width=10,height=4.065)
#  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
#  for (data in test2){
#    test<-data@Polygons
#    for (data1 in data@Polygons){
#      test4<-data1@coords
#      polygon(test4, col=scales::alpha('red', 0.5))
#    }
#  }
#  dev.off()
#}

p_in_polygon<-function(nvert, vertx, verty, point){
  c=F
  j=1
  testx <- point[1]
  testy <- point[2]
  for (i in 1:(nvert-2)){
    if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ){
      c = !c
    }
    j=i
  }
  return(c)
}

eez_points <- NULL
eez_points1 <- NULL
passed <- NULL
for (data in test2){
  test<-data@Polygons
  for (data1 in data@Polygons){
    test4<-data1@coords
    for (long in seq( floor(min(test4[,1]))-0.5,floor(max(test4[,1]))+0.5, 1)){
      for (lat in seq( floor(min(test4[,2]))-0.5,floor(max(test4[,2]))+0.5, 1) ){
        test_func <- p_in_polygon(dim(test4)[1], test4[,1], test4[,2],c(long, lat))
        if ( !(paste(lat,long,sep='_') %in% passed) & paste(lat,long,sep='_') %in% bray_curtis$cell & test_func==T){
          eez_points<-append(eez_points, paste(lat,long,sep='_') )
          eez_points1 <- rbind(eez_points1, c(lat,long))
        }
        passed <- append(passed, paste(lat,long,sep='_'))
      }
    }
  }
}
saveRDS(eez_points, paste(model0,'eez_points.rds', sep='_'))
saveRDS(eez_points1, paste(model0,'eez_points1.rds', sep='_'))


pdf(family="Helvetica",paste(model0,'EEZ_map.pdf', sep='_'), width=10,height=4.065)
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
points(x=eez_points1[,2], y=eez_points1[,1], col='red',bg= 'red',pch=15,cex=0.4)
dev.off()
# pdf(family="Helvetica",'test.pdf')
# plot(-168.7, -14, col='red', cex= 5, xlim=c(-175,-160), ylim=c(-20,-10), pch=15)
# polygon(test4, col='blue')
# points(-168.7, -14, col='red')
# dev.off()
