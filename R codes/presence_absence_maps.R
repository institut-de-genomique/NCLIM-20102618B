library('maptools')
cambria <- readRDS('cambria.rds')
source('hide_arctic.R')
source('axis_map0.R')
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')
upper.half.circle <- function(x,y,r,nsteps=100,...){
  rs <- seq(0,pi,len=nsteps)
  xc <- x+r*cos(rs)
  yc <- y+r*sin(rs)
  polygon(xc,yc, border = NA,...)
}

lower.half.circle <- function(x,y,r,nsteps=100,...){
  rs <- seq(0,pi,len=nsteps)
  xc <- x-r*cos(rs)
  yc <- y-r*sin(rs)
  polygon(xc,yc,border = NA,...)
}

df <- readRDS('Genocenoses_env_parameters_woa.rds')
frac <- '180-2000'
clusts <- c(5, 8)
for(clust in clusts){
  df1 <-df[df$Fraction==frac,]
  df1$Genocenose <- as.numeric(df1$Genocenose==clust)+1
  colors <- c('red','green' )
  dcm <- grepl('DCM', df1$Station)
  sur <- grepl('SUR', df1$Station)
  pdf(family="Helvetica",paste(frac,'_', clust,'.pdf', sep=''),width=10,height=4.065)
  par(mar=c(0,0,0,0))
  maps::map(database="world",fill=T,col="grey80",border="gray80",xpd=TRUE)
  plot(coastline,lwd=0.0475, col='black', add=T)
  for (i in 1:length(df1$Long[sur])){
    upper.half.circle(df1$Long[sur][i], df1$Lat[sur][i], col=colors[df1$Genocenose[sur][i]], r=3)
  }
  for (i in 1:length(df1$Long[dcm])){
    lower.half.circle(df1$Long[dcm][i], df1$Lat[dcm][i], col=colors[df1$Genocenose[dcm][i]], r=3)
  }
  
  #points(df1$Long, df1$Lat, col=colors[df1$Genocenose], pch=19, cex=20)
  hide_arctic()
  axis_map0()
  dev.off()
}
