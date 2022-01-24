library('maptools')
library('RColorBrewer')
axis_map0 <- function(){
  axis(1, pos=-84.5, at=c(-179.5,seq(-120,180,60)), cex.axis=1.5, tck=-0.03, las=1, font=2,
       labels=c('180°W',  '120°W', '60°W',  '0',  '60°E', '120°E',  '180°E'), lwd= 1.4, padj = 0)
  axis(2, pos=-179.5, at=c( -84.5, seq(-60,60,30)), cex.axis=1.5, tck=-0.03, las= 1,  font=2,
       labels=c(NA, '60°S', '30°S', '0', '30°N', '60°N'), lwd= 1.4, padj=0.45)
  axis(3, pos=60, at=c(-179.5,seq(-150,180,30)),labels = F, tck=0, lwd= 1.4)
  axis(4, pos=180, at=c( -84.5, seq(-60,60,30)),labels = F,tck=0, lwd= 1.4)
}
coastline <- readShapeSpatial('ne_10m_coastline/ne_10m_coastline.shp')

