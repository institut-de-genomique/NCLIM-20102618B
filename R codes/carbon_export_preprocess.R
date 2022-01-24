library('ncdf4')
source('axis_map0.R')
source('hide_arctic.R')
cexp <- nc_open('Henson2012_Figure4_ncs/export_Henson2012.nc')
lon <- ncvar_get(cexp, 'lon')
lon[lon< -180]<- lon[lon< -180]+360
lat <- ncvar_get(cexp, 'lat')
lat <- ncvar_get(cexp, 'lat')
dt_cexp <- ncvar_get(cexp, 'EPCHn') 
dt_cexpv <- as.vector(dt_cexp)
lons <- rep(lon, length(lat))
lats <- NULL
for (l in lat){
  lats <- append(lats, rep(l, length(lon)))
}
dt_cexps <- data.frame('Latitude'=lats, 'Longitude'=lons, 'Export.flux'=dt_cexpv)
dt_cexps$col <- round(99*(dt_cexps$Export.flux-min(dt_cexps$Export.flux, na.rm=T))/(max(dt_cexps$Export.flux, na.rm=T)
                                                                                     -min(dt_cexps$Export.flux, na.rm=T)))+1
saveRDS(dt_cexps, 'Henson_2012_export.rds')
col_scale <- colorRampPalette(c('darkorchid4', 'blue', 'green', 'orange', 'red'))(100)
pdf(family="Helvetica",file='Henson_2012_export.pdf',
    width=10,height=4.065)
maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
plot(coastline,lwd=0.0475, col='black', add=T, cex=0.284)
points(x = dt_cexps$Longitude,y =  dt_cexps$Latitude, cex=0.15, pch=15, col=col_scale[dt_cexps$col])
hide_arctic()
axis_map0()
dev.off()

names <- c('eppley', 'laws', 'schlitzer')
for (na in names){
  dat <- nc_open(paste(na, '_export.nc', sep=''))
  lat_s <- ncvar_get(dat, 'ETOPO60Y')
  lon_s <- ncvar_get(dat, 'ETOPO60X')
  lon_s[lon_s> 179.5]<-lon_s[lon_s> 179.5]-360 
  cexp_s <- ncvar_get(dat, names(dat$var))
  dt_cexpv_s <- as.vector(cexp_s)
  lons_s <- rep(lon_s, length(lat_s))
  lats_s <- NULL
  for (l in lat){
    lats_s <- append(lats_s, rep(l, length(lon_s)))
  }
  dt_cexps_s <- data.frame('Latitude'=lats_s, 'Longitude'=lons_s, 'Export.flux'=dt_cexpv_s)
  dt_cexps_s$col <- round(99*(dt_cexps_s$Export.flux-min(dt_cexps_s$Export.flux, na.rm=T))/(max(dt_cexps_s$Export.flux, na.rm=T)
                                                                                      -min(dt_cexps_s$Export.flux, na.rm=T)))+1
  saveRDS(dt_cexps_s, paste(na,'_export.rds', sep=''))
  col_scale <- colorRampPalette(c('darkorchid4', 'blue', 'green', 'orange', 'red'))(100)
  pdf(family="Helvetica",file=paste(na,'_export.pdf', sep=''),
      width=10,height=4.065)
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T, cex=0.284)
  points(x = dt_cexps_s$Longitude,y =  dt_cexps_s$Latitude, cex=0.15, pch=15, col=col_scale[dt_cexps_s$col])
  hide_arctic()
  axis_map0()
  dev.off()
}
