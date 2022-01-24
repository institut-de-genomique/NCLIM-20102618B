library('vegan')
library('parallel')
library('maptools')
library('ncdf4')

cambria <- readRDS('cambria.rds')
BGCP_shp <- readShapeSpatial('longhurst_v4_2010/Longhurst_world_v4_2010.shp')
png("longhurst.png", width=10,height=4.065, res=200, units = 'in')
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
for (u in BGCP_shp@polygons){
  for (v in u@Polygons){
    polygon(v@coords)
  }
}
dev.off()

test <- ncdf4::nc_open('Time_Varying_Biomes.nc')
u <- ncvar_get(test, 'MeanBiomes')
u[is.nan(u)]<-NA
biomes <- as.vector(u)
lt <- ncvar_get(test, 'lat')
lo <- ncvar_get(test, 'lon')

v <-t(u)
v[is.na(v)]<--1
source('axis_map.R')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
pdf(family="Helvetica",file='McKingley.pdf',width=10,height=4.065)
par(mar=c(25,10,25,10))
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T)
for (j in unique(v[v!=-1])){
  h<-v
  h[v!=j]=0
  contour(lo, lt, h, col = 'black',add=T, drawlabels = F, cex=0.284)
}
axis_map()
dev.off()


grd_bio <- expand.grid(lat=lt, lon=lo)
grd_bio$biome <- biomes
grd_bio$cell <- paste(grd_bio$lat, grd_bio$lon, sep='_')

reyg_bgcp19 <- read.csv2('Reygondeau/BGCP_2019_REYGONDEAU.csv', header = T, sep=',')
reyg_bgcp19$cell <- paste(reyg_bgcp19$Latitude, reyg_bgcp19$Longitude, sep='_')
reyg_bgcp19$BGCP[is.nan(reyg_bgcp19$BGCP)]<-NA

reyg_bgcp <- read.csv2('Reygondeau/PROVINCE_REYGONDEAU.csv', header = F, sep = ',')
reyg_bgcp_vec <- as.vector(as.matrix(reyg_bgcp))
reyg_bgcp_vec[is.nan(reyg_bgcp_vec)]<-NA
reyg_bio <- read.csv2('Reygondeau/BIOME_REYGONDEAU.csv', header = F, sep = ',')
reyg_bio_vec <- as.vector(as.matrix(reyg_bio))
reyg_bio_vec[is.nan(reyg_bio_vec)]<-NA
grd_bio$reyg_bio <- reyg_bio_vec

reyg_bgcp19_1 <- matrix(NA, nrow=180, ncol = 360)
c_lt <-1
for (lat in lt){
  c_lo<-1
  for (lon in lo){
    lts <- c(lat-0.25, lat+0.25)
    los <- c(lon-0.25, lon+0.25)
    gr <- expand.grid(lati=lts, loni=los)
    cells_gr <- paste(gr$lati, gr$loni, sep='_')
    vec_19 <- reyg_bgcp19$BGCP[reyg_bgcp19$cell %in% cells_gr]
    x <- unique(vec_19)
    if (length(x)==1){
      reyg_bgcp19_1[c_lt, c_lo]<-x
    } else{
      y <- data.frame(table(vec_19, useNA = 'ifany'))
      reyg_bgcp19_1[c_lt, c_lo]<-x[which.max(y$Freq)]
    }
    c_lo<-c_lo+1
  }
  c_lt<-c_lt+1
}
saveRDS(reyg_bgcp19_1, 'BGCP_2019_REYGONDEAU_mat_1.rds')
reyg_bgcp19_1_vec <- as.vector(reyg_bgcp19_1)
reyg_bgcp19_1_vec[is.nan(reyg_bgcp19_1_vec)]<-NA
grd_bio$reyg_bgcp19 <- reyg_bgcp19_1

saveRDS(grd_bio, 'biogeographies.rds')

reyg_bgcp19_1 <-readRDS('BGCP_2019_REYGONDEAU_mat_1.rds')
to_plot <-list(reyg_bgcp19_1, reyg_bgcp, reyg_bio)
names <- c('Reygondeau_BGCP19.pdf','Reygondeau_BGCP.pdf','Reygondeau_BIOME.pdf' )
for (i in 1:3){
  pdf(family="Helvetica",file=names[i],width=10,height=4.065)
  par(mar=c(25,10,25,10))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  v <- t(to_plot[[i]])
  for (j in unique(v[v!=-1])){
    h<-v
    h[v!=j]=0
    contour(lo, lt, h, col = 'black',add=T, drawlabels = F, cex=0.284)
  }
  axis_map()
  dev.off()
}
